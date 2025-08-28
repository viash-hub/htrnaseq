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
      | map { id, state -> 
        // The argument names for this workflow and the htrnaseq workflow may overlap
        // here, we store a copy in order to make sure to not accidentally overwrite the state.
        def new_state = state + [
          "star_output_dir_workflow": state.star_output_dir,
          "nrReadsNrGenesPerChrom_dir_workflow": state.nrReadsNrGenesPerChrom_dir,
          "star_qc_metrics_dir_workflow": state.star_qc_metrics_dir,
          "eset_dir_workflow": state.eset_dir,
          "f_data_dir_workflow": state.f_data_dir,
          "p_data_dir_workflow": state.p_data_dir
        ]
        return [id, new_state]
      }
      | listInputDir.run(
        fromState: [
          "input": "input",
          "pools": "pools",
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
      | view {"Pool inputs after listing directory: $it"}
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
    // So, we should combine everything together
    results_publish_ch = htrnaseq_ch
      | combine(save_params_ch)
      | map {new_id, grouped_ch_state, save_params_id, save_params_state ->
        def new_state = grouped_ch_state + ["run_params": save_params_state.run_params]
        return [new_id, new_state]
      }
      | toSortedList
      | map{ vs ->
          def states = vs.collect{it[1]}

          // The STAR output is a directory for each well in a plate (or pool of plates).
          // The wells are grouped into a directory per sample. The name of this directory should
          // match the sample_id.
          def star_output_samples = states.collectMany{state -> 
            state.star_output.collect{
              def star_sample_dir = it.parent
              assert star_sample_dir.name == state.sample_id: "Unexpected state: the parent directory of STAR output \
                path '${it}' should match with the sample ID ${sample_id}"
              star_sample_dir
            }
          }
          def new_state = [
            "star_output": star_output_samples,
          ]

          // Keys for which the values should be the same across samples
          def state_keys_unique = [
            "html_report",
            "project_id",
            "experiment_id",
            "star_output_dir_workflow",
            "nrReadsNrGenesPerChrom_dir_workflow",
            "star_qc_metrics_dir_workflow",
            "eset_dir_workflow",
            "f_data_dir_workflow",
            "p_data_dir_workflow",
            "f_data",
            "run_params"
          ]
          def state_unique_keys = state_keys_unique.inject([:]) { state_to_update, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "State error: values for argument $argument_name should be the same across states. \
                                                 Argument values: $argument_values"
            // take the unique value from the set (there is only one)
            def argument_value
            argument_values.each { argument_value = it }
            state_to_update + [(argument_name): argument_value]
          }

          // Keys that just require gathering of values across samples
          def state_keys_collect = [
            "nrReadsNrGenesPerChrom",
            "star_qc_metrics",
            "eset",
            "p_data",
          ]
          def state_collect = state_keys_collect.collectEntries{ key_ ->
            [key_, states.collect{it.get(key_)}]
          }

          new_state = new_state + state_unique_keys + state_collect
          [states[0].run_id, new_state]  
      }
      | publish_results.run(
        fromState: { id, state ->
          def prefix = "${state.project_id}/${state.experiment_id}/data_processed/${date}_htrnaseq_${version}"

          println("Publising results to ${params.results_publish_dir}/${prefix}")

          [ 
            // Inputs
            star_output: state.star_output,
            nrReadsNrGenesPerChrom: state.nrReadsNrGenesPerChrom,
            star_qc_metrics: state.star_qc_metrics,
            eset: state.eset,
            f_data: state.f_data,
            p_data: state.p_data,
            html_report: state.html_report,
            run_params: state.run_params,
            // Output locations
            run_params_output: "${prefix}/${state.run_params.name}",
            html_report_output: "${prefix}/${state.html_report.name}", 
            star_output_dir: "${prefix}/${state.star_output_dir_workflow}",
            nrReadsNrGenesPerChrom_dir: "${prefix}/${state.nrReadsNrGenesPerChrom_dir_workflow}",
            star_qc_metrics_dir: "${prefix}/${state.star_qc_metrics_dir_workflow}",
            eset_dir: "${prefix}/${state.eset_dir_workflow}",
            f_data_dir: "${prefix}/${state.f_data_dir_workflow}",
            p_data_dir: "${prefix}/${state.p_data_dir_workflow}"
          ]
        },
        toState: { id, result, state -> result },
        directives: [
          publishDir: [
            path: "${params.results_publish_dir}", 
            overwrite: false,
            mode: "copy"
          ]
        ]
      )
      | setState([
          "star_output_dir",
          "nrReadsNrGenesPerChrom_dir",
          "star_qc_metrics_dir",
          "eset_dir",
          "f_data_dir",
          "p_data_dir",
        ]
      )

    fastq_publish_ch = htrnaseq_ch
      // The output from the htrnaseq workflow is on sample (i.e. pool) level
      // Multiple sequencing runs may have contributes to the FASTQ files from this pool.
      // So the fastq_output is a list of directories, one for each run.
      // We assume that the names of the folders containing the FASTQ files are equal to the pool names.
      | flatMap {id, state ->
          state.fastq_output.collect{fastq_dir ->
            def run_id = fastq_dir.name
            def new_id = "${run_id}/${state.sample_id}"
            def new_state = [
              "fastq_output": fastq_dir.listFiles(),
              "sample_id": state.sample_id,
              "run_id": run_id,
              "output": "${run_id}/${date}_htrnaseq_${version}/${state.sample_id}".toString()
            ]
            [new_id, new_state]
          }
      }
      // A folder containing the FASTQ files from a certain pool may be present in the state from
      // multiple samples; if that pool contributed to the data from those samples.
      // Those FASTQ files will only be published once by filtering out the duplicate events here.
      | unique{it[0]}
      | publish_fastqs.run(
        fromState: [
          "input": "fastq_output",
          "output": "output",
        ],
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
    results_publish_ch

}

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("HT-RNAseq version to be used: ${version}")
  return version
}