def date = new Date().format('yyyyMMdd_hhmmss')

def viash_config = java.nio.file.Paths.get("${moduleDir}/_viash.yaml")
def version = get_version(viash_config)

workflow run_wf {
  take:
    raw_ch

  main:
    input_ch = raw_ch
      // Use the ID of the event as the run ID
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
    
    /*
      The provided input is a single directory per input event. However, the well demultiplexing
      requires a list of FASTQ files to be demultiplexed. The `listInputDir` subworkflow lists the
      content of the input directory and deduces the pool (sample) IDs from the file names.
      When the FASTQ files are split per lane, they are gathered into the `r1_output` and `r2_output`
      lists of a single event so that the lanes can be demultiplexed in parallel and joined afterwards.
    */
    demultiplex_ch = input_ch
      | listInputDir.run(
        fromState: [
          "input": "input",
          "pools": "pools",
        ],
        toState: [
          "output_r1": "r1_output",
          "output_r2": "r2_output",
          "sample_id": "sample_id",
        ]
      )
      /*
        `listInputDir` puts the pool ID as the event ID (slot 0 from the tuple). The same pool
        can be present in multiple sequencing runs, so the run ID is added in order to keep the
        events unique; otherwise events with duplicate IDs are dropped.
      */
      | map {id, state -> ["${state.run_id}/${id}".toString(), state]}
      // Each pool is demultiplexed into its wells separately
      | well_demultiplex.run(
        fromState: [
          "input_r1": "output_r1",
          "input_r2": "output_r2",
          "barcodesFasta": "barcodesFasta",
        ],
        toState: { id, result, state ->
          def all_fastq = result.output_r1 + result.output_r2
          def output_dir = all_fastq.collect{it.parent}.unique()
          assert output_dir.size() == 1, "Expected output from well demultiplexing (id $id) to reside into one directory. Found: $output_dir"
          def new_state = state + [
            "output_r1": result.output_r1,
            "output_r2": result.output_r2,
            "fastq_output_directory": output_dir[0],
          ]
          return new_state
        }
      )

    // Publish the demultiplexed FASTQ files, one directory per pool.
    fastq_publish_ch = demultiplex_ch
      | map {id, state ->
        def fastq_prefix = "${state.project_id}/${state.experiment_id}/${state.run_id}/" +
          "${date}_well_demultiplex_${version}/${state.sample_id}"
        def new_state = state + ["fastq_prefix": fastq_prefix.toString()]
        return [id, new_state]
      }
      | publish_fastqs.run(
        fromState: { id, state ->
          println("Publishing fastqs to ${params.fastq_publish_dir}/${state.fastq_prefix}")

          [
            "input": state.output_r1 + state.output_r2,
            "output": state.fastq_prefix,
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
    /*
      The demultiplexing outputs an event per pool, while this runner outputs events on the
      level of the input events (i.e. per sequencing run). The events are grouped back together
      here and the YAML file with the run parameters is added to the output.
    */
    fastq_publish_ch
      | map {id, state -> [state.run_id, state]}
      | groupTuple(by: 0, sort: "hash")
      | combine(save_params_ch)
      | map {run_id, states, _save_params_id, save_params_state ->
        def new_state = [
          "run_params": save_params_state.run_params,
          "_meta": ["join_id": run_id],
        ]
        return [run_id, new_state]
      }
}

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("Version to be used: ${version}")
  return version
}