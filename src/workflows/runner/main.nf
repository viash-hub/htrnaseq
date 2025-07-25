def date = new Date().format('yyyyMMdd_hhmmss')

def viash_config = java.nio.file.Paths.get("${moduleDir}/_viash.yaml")
def version = get_version(viash_config)

workflow run_wf {
  take:
    raw_ch

  main:
    input_ch = raw_ch
      // List the FASTQ files per input directory
      // Be careful: an event per lane is created!
      | map {id, state ->
        def new_state = state + ["run_id": id]
        return [id, new_state]
      }

    save_params_ch = input_ch
      | toSortedList()
      | map { states ->
        def new_id = "save_params"
        def all_states = states.collect{it[1]}
        def run_params_output_templates = all_states.collect{it.run_params}
        assert run_params_output_templates.unique().size() == 1: "The value for the 'run_params' parameter is not the same across runs."
        def new_state = ["run_params": run_params_output_templates[0], "all_states": all_states]
        return [new_id, new_state]
      }

      | save_params.run(
        key: "save_params_runner",
        fromState: {id, state ->

          def convertPaths
          convertPaths = { value ->
            if (value instanceof java.nio.file.Path)
              return value.toUriString()
            else if (value instanceof List)
              return value.collect { convertPaths(it) }
            else if (value instanceof Collection)
              throw new UnsupportedOperationException("Collections other than Lists are not supported")
            else
              return value
          }
          
          // Apply conversion to all state values
          def convertedState = state.all_states.collect{it.collectEntries { k, v -> [(k): convertPaths(v)] }}
          
          def yaml = new org.yaml.snakeyaml.Yaml()
          def yamlString = yaml.dump(convertedState)
          def encodedYaml = yamlString.bytes.encodeBase64().toString()
          
          return [
            "id": id,
            "params_yaml": encodedYaml,
            "output": state.run_params
          ]
        },
        toState: ["run_params": "output"]
      )
    
    htrnaseq_ch = input_ch
      | listInputDir.run(
        fromState: [
          "input": "input",
          "ignore": "ignore",
        ],
        toState: { id, state, result ->
          def clean_state = state.findAll{ it.key != "input" }
          clean_state + result
        }
      )
      // ListInputDir puts the sample_id as the event ID (slot 0 from the tuple).
      // The sample_id was inferred from the start of the file name,
      // and it can be used to group the FASTQ files, because an input folder 
      // can contain input files from multiple samples (pools). Additionally,
      // there might be multiple FASTQs for a single sample that correspond to the
      // lanes. So the fastq files must be gathered across lanes and input folders
      // in order to create an input lists for R1 and R2.
      // The ID of the event here is important! It determines the name of the output
      // folders for the FASTQ files and these folders are published as-is later.
      // The folder where the FASTQ files are stored in should be named after the run ID.
      | map {id, state -> ["${state.sample_id}/${state.run_id}".toString(), state]}
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def new_r1 = states.collect{it.r1_output}
        def new_r2 = states.collect{it.r2_output}
        // This assumes that, except for r1 and r2, 
        // the keys across the grouped states are the same.
        // TODO: this can be asserted.
        def new_state = states[0] + [
          "r1": new_r1,
          "r2": new_r2,
        ]
        return [id, new_state]
      }
      | htrnaseq.run(
        args: [
          f_data: 'fData/$id.txt',
          p_data: 'pData/$id.txt',
          star_output: 'star_output/$id/*',
          fastq_output: 'fastq/*',
          eset: 'esets/$id.rds',
          nrReadsNrGenesPerChrom: 'nrReadsNrGenesPerChrom/$id.txt',
          star_qc_metrics: 'starLogs/$id.txt',
          html_report: "report.html",
          run_params: null
        ],
        fromState: [
          input_r1: "r1",
          input_r2: "r2",
          barcodesFasta: "barcodesFasta",
          genomeDir: "genomeDir",
          annotation: "annotation",
          umi_length: "umi_length",
          sample_id: "sample_id",
        ],
        toState: { id, result, state -> state + result }
      )

      // The HT-RNAseq workflow outputs multiple events, one per 'pool' (usually a plate)
      // but for publishing the results, this is not handy because we want to use the $id
      // variable as a pointer to the target data.
      //
      // So, we should combine everything together
      //
      // project_id / experiment_id / "data_processed" / date_workflow
    grouped_ch = htrnaseq_ch
      | toSortedList
      | map{ vs ->
          def all_fastqs
          [
            vs[0][1].run_id, // The original ID
            [
              star_output: reduce_paths(vs.collect{ it[1].star_output }.flatten()),
              nrReadsNrGenesPerChrom: reduce_paths(vs.collect{ it[1].nrReadsNrGenesPerChrom }),
              star_qc_metrics: reduce_paths(vs.collect{ it[1].star_qc_metrics }),
              eset: reduce_paths(vs.collect{ it[1].eset }),
              f_data: reduce_paths(vs.collect{ it[1].f_data }),
              p_data: reduce_paths(vs.collect{ it[1].p_data }),
              fastq_output: vs.collect{ it[1].fastq_output }.flatten().unique(),
              html_report: vs.collect{ it[1].html_report }[0],  // The report is for all pools
              plain_output: vs.collect{ it[1].plain_output }[0],
              project_id: vs.collect{ it[1].project_id }[0],
              experiment_id: vs.collect{ it[1].experiment_id }[0]
            ]
          ]
        }

    grouped_with_params_list_ch = grouped_ch.combine(save_params_ch)
      | map {new_id, grouped_ch_state, save_params_id, save_params_state ->
        def new_state = grouped_ch_state + ["run_params": save_params_state.run_params]
        return [new_id, new_state]
      
      }

    results_publish_ch = grouped_with_params_list_ch
      | publish_results.run(
        fromState: { id, state ->
          def project = (state.plain_output) ? id : "${state.project_id}"
          def experiment = (state.plain_output) ? id : "${state.experiment_id}"
          def id0 = "${project}/${experiment}"
          def id1 = (state.plain_output) ? id : "${id0}/data_processed/${date}"
          def id2 = (state.plain_output) ? id : "${id1}_htrnaseq_${version}"

          if (id == id2) {
            println("Publising results to ${params.results_publish_dir}")
          } else {
            println("Publising results to ${params.results_publish_dir}/${id2}")
          }

          [
            star_output: state.star_output,
            star_output: state.star_output,
            nrReadsNrGenesPerChrom: state.nrReadsNrGenesPerChrom,
            star_qc_metrics: state.star_qc_metrics,
            eset: state.eset,
            f_data: state.f_data,
            p_data: state.p_data,
            html_report: state.html_report,
            run_params: state.run_params,
            output: "${id2}"
          ]
        },
        toState: { id, result, state -> state },
        directives: [
          publishDir: [
            path: "${params.results_publish_dir}", 
            overwrite: false,
            mode: "copy"
          ]
        ]
      )

    fastq_publish_ch = grouped_ch
      | flatMap{id, state ->
        def new_states = state.fastq_output.collect{fastq_dir -> 
          def run_id = fastq_dir.name // The folder name corresponds to the run
          def sample = fastq_dir.parent.name // The parent folder name should be the sample
          def new_id = "${run_id}/${sample}"
          def fastq_files = fastq_dir.listFiles()
          def new_state = [
            "fastq_output": fastq_files,
            "sample_id": sample,
            "run_id": run_id 
          ]
          return [new_id, new_state]
        }
        return new_states
      }
      | publish_fastqs.run(
        fromState: { id, state ->
          def id0 = state.run_id
          def id1 = (state.plain_output) ? id : "${id0}/${date}"
          def id2 = (state.plain_output) ? id : "${id1}_htrnaseq_${version}/${state.sample_id}"

          if (id == id2) {
            println("Publising fastqs to ${params.fastq_publish_dir}")
          } else {
            println("Publising fastqs to ${params.fastq_publish_dir}/${id2}")
          }

          [
            input: state.fastq_output,
            output: "${id2}",
          ]
        },
        toState: { id, result, state -> state },
        directives: [
          publishDir: [
            path: "${params.fastq_publish_dir}", 
            overwrite: false,
            mode: "copy"
          ]
        ]
      )

  emit:
    grouped_ch
      | map{ id, state -> [ id, [ _meta: [ join_id: state.run_id ] ] ] }
}

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("HT-RNAseq version to be used: ${version}")
  return version
}

/*
 * This function uses a heuristic to group a list of paths so that the level of nesting
 * of IDs is represented in the output.
 *
 * We iterative of the path sections (subfolders) from the last (file) the first (root node).
 * The first path segment that is common across all 'events' is the cutoff. We cutoff the paths
 * at this level, select the unique elements from the list and use that as input for the next step.
 *
 * An optional offset allows one to shift the cutoff left or right.
 */
def reduce_paths(paths, offset = 0) {
  def path_length = paths.collect{ it.getNameCount() }[0]

  def unique_list = (path_length-1..0).collectEntries { i ->
    [ (i): paths.collect{ it.getName(i) }.unique().size() ]
  }

  def cutoff = unique_list.find{ it.value == 1 }.key

  def grouped_paths = paths.collect{ f -> "/" + f.subpath(0, cutoff+1+offset) }.unique()

  println("")
  println("Detecting the common path section to pass to the next step:")
  print("  From: ")
  print paths
  println("")
  print("  To:   ")
  print grouped_paths
  println("")

  return grouped_paths
}
