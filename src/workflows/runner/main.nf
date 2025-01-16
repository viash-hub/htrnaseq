def date = new Date().format('yyyyMMdd_hhmmss')

def viash_config = java.nio.file.Paths.get("$projectDir/../../../../").toAbsolutePath().normalize().toString() + "/_viash.yaml"
def version = get_version(viash_config)

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

      // Pass the original ID in the state to pick it up later
      | map{ id, state -> [ id, state + [ run_id: id ] ] }

      | listInputDir.run(
        fromState: [ input: "input" ],
        toState: { id, state, result ->
          def clean_state = state.findAll{ it.key != "input" }
          clean_state + result
        }
      )

      | htrnaseq.run(
        args: [
          f_data: 'fData/$id.txt',
          p_data: 'pData/$id.txt',
          star_output: 'star_output/$id/*',
          fastq_output_r1: 'fastq/*_R1_001.fastq',
          fastq_output_r2: 'fastq/*_R1_001.fastq',
          eset: 'esets/$id.rds',
          nrReadsNrGenesPerChrom: 'nrReadsNrGenesPerChrom/$id.txt',
          star_qc_metrics: 'starLogs/$id.txt',
          html_report: "report.html"
        ],
        fromState: [
          input_r1: "r1",
          input_r2: "r2",
          barcodesFasta: "barcodesFasta",
          genomeDir: "genomeDir",
          annotation: "annotation",
        ],
        toState: { id, result, state -> state + result }
      )

      // The HT-RNAseq workflow outputs multiple events, one per 'pool' (usually a plate)
      // but for publishing the results, this is not handy because we want to use the $id
      // variable as a pointer to the target data.
      //
      // So, we should combine everything together
      //
      // project_id / experiment_id / date_workflow

      | toSortedList

      | map{ vs ->
          def all_fastqs
          [
            vs[0][1].run_id, // The original ID
            [
              star_output: reduce_paths(vs.collect{ it[1].star_output }.flatten()),
              fastq_output_r1: reduce_paths(vs.collect{ it[1].fastq_output_r1 }.flatten(), 1),
              fastq_output_r2: reduce_paths(vs.collect{ it[1].fastq_output_r2 }.flatten(), 1),
              nrReadsNrGenesPerChrom: reduce_paths(vs.collect{ it[1].nrReadsNrGenesPerChrom }),
              star_qc_metrics: reduce_paths(vs.collect{ it[1].star_qc_metrics }),
              eset: reduce_paths(vs.collect{ it[1].eset }),
              f_data: reduce_paths(vs.collect{ it[1].f_data }),
              p_data: reduce_paths(vs.collect{ it[1].p_data }),
              html_report: vs.collect{ it[1].html_report }[0],  // The report is for all pools
              plain_output: vs.collect{ it[1].plain_output }[0],
              project_id: vs.collect{ it[1].project_id }[0],
              experiment_id: vs.collect{ it[1].experiment_id }[0]
            ]
          ]
        }

      | publish_results.run(
        fromState: { id, state ->
          def project = (state.plain_output) ? id : "${state.project_id}"
          def experiment = (state.plain_output) ? id : "${state.experiment_id}"
          def id0 = "${project}/${experiment}"
          def id1 = (state.plain_output) ? id : "${id0}/${date}"
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

      | publish_fastqs.run(
        fromState: { id, state ->
          def id0 = "${id}"
          def id1 = (state.plain_output) ? id : "${id0}/${date}"
          def id2 = (state.plain_output) ? id : "${id1}_htrnaseq_${version}"

          println(state.plain_output)

          if (id == id2) {
            println("Publising fastqs to ${params.fastq_publish_dir}")
          } else {
            println("Publising fastqs to ${params.fastq_publish_dir}/${id2}")
          }

          [
            input_r1: state.fastq_output_r1,
            input_r2: state.fastq_output_r2,
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
    output_ch
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
