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
