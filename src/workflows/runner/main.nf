def date = new Date().format('yyyyMMdd_hhmmss')

def viash_config = java.nio.file.Paths.get("$projectDir/../../../../").toAbsolutePath().normalize().toString() + "/_viash.yaml"
def version = get_version(viash_config)

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

      | listInputDir.run(
        fromState: [ input: "input" ],
        toState: { id, state, result ->
          def clean_state = state.findAll{ it.key != "input" }
          clean_state + result
        }
      )

      | htrnaseq.run(
        fromState: [
          input_r1: "r1",
          input_r2: "r2",
          barcodesFasta: "barcodesFasta",
          genomeDir: "genomeDir",
          annotation: "annotation"

        ],
        toState: { id, result, state -> result }
      )

      | niceView()

      | publish_results.run(
        fromState:
          [
            star_output: "star_output", 
            fastq_output_r1: "fastq_output_r1",
            fastq_output_r2: "fastq_output_r2",
            star_output: "star_output",
            nrReadsNrGenesPerChrom: "nrReadsNrGenesPerChrom",
            star_qc_metrics: "star_qc_metrics",
            eset: "eset",
            f_data: "f_data",
            p_data: "p_data",
            html_report: "html_report",
          ],
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
          [
            input_r1: state.fastq_output_r1,
            input_r2: state.fastq_output_r2,
            output: "${id}/",
          ]
        },
        toState: { id, result, state -> [:] },
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
      | map{ id, state -> [ id, state + [ _meta: [ join_id: "run" ] ] ] }
}

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("Version to be used: ${version}")
  return version
}
