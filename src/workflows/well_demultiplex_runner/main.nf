def date = new Date().format('yyyyMMdd_hhmmss')

def viash_config = java.nio.file.Paths.get("$projectDir/../../../../").toAbsolutePath().normalize().toString() + "/_viash.yaml"
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
    
    demultiplex_ch = input_ch
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
      // Group FASTQ files by sample_id and run_id
      | map {id, state -> ["${state.sample_id}/${state.run_id}".toString(), state]}
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def new_r1 = states.collect{it.r1_output}
        def new_r2 = states.collect{it.r2_output}
        // This assumes that, except for r1 and r2, 
        // the keys across the grouped states are the same.
        def new_state = states[0] + [
          "r1": new_r1,
          "r2": new_r2,
        ]
        return [id, new_state]
      }

      | well_demultiplex.run(
        args: [
          output_r1: 'fastq/$id/*_R1_001.fastq',
          output_r2: 'fastq/$id/*_R2_001.fastq',
        ],
        fromState: [
          input_r1: "r1",
          input_r2: "r2",
          barcodesFasta: "barcodesFasta",
          sample_id: "sample_id",
        ],
        toState: { id, result, state -> state + result }
      )

    // Group all output by run_id for simplified publishing
    grouped_ch = demultiplex_ch
      | toSortedList()
      | map{ vs ->
          [
            vs[0][1].run_id, // The original ID
            [
              output_r1: vs.collect{ it[1].output_r1 }.flatten(),
              output_r2: vs.collect{ it[1].output_r2 }.flatten(),
              run_params: vs.collect{ it[1].run_params }[0],
              plain_output: vs.collect{ it[1].plain_output }[0],
              project_id: vs.collect{ it[1].project_id }[0],
              experiment_id: vs.collect{ it[1].experiment_id }[0]
            ]
          ]
        }

    grouped_with_params_ch = grouped_ch.combine(save_params_ch)
      | map {new_id, grouped_ch_state, save_params_id, save_params_state ->
        def new_state = grouped_ch_state + ["run_params": save_params_state.run_params]
        return [new_id, new_state]
      }

    // Publish demultiplexed FASTQ files
    fastq_publish_ch = grouped_ch
      | flatMap{id, state ->
      // Extract unique sample IDs from the output files
      def files = state.output_r1 + state.output_r2
      if (files.isEmpty()) return []
      
      // Get the sample directories (they contain the actual FASTQ files)
      def sampleDirs = files.collect { file -> 
        file.getParent() 
      }.unique()
      
      // For each sample directory, create a publishing event
      return sampleDirs.collect { dir ->
        def dirPath = dir.toString()
        def pathParts = dirPath.tokenize('/')
        def sampleIdx = pathParts.findIndexOf { it.contains("_") || it ==~ /[A-Za-z0-9]+/ }
        def sampleId = sampleIdx >= 0 ? pathParts[sampleIdx] : "unknown"
        
        def new_state = [
        "fastq_output": dir,
        "sample_id": sampleId,
        "run_id": id,
        "project_id": state.project_id,
        "experiment_id": state.experiment_id,
        "plain_output": state.plain_output
        ]
        
        ["${id}/${sampleId}", new_state]
      }
      }
      | publish_fastqs.run(
      fromState: { id, state ->
        def projectPath = state.project_id ?: "unknown_project"
        def experimentPath = state.experiment_id ?: "unknown_experiment"
        def samplePath = state.sample_id
        def runPath = state.run_id
        
        def publishPath = "${projectPath}/${experimentPath}/${runPath}/${date}_well_demultiplex_${version}/${samplePath}"
        
        println("Publishing fastqs to ${params.fastq_publish_dir}/${publishPath}")

        [
        "input": state.fastq_output,
        "output": publishPath
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

def get_version(input) {
  def inputFile = file(input)
  if (!inputFile.exists()) {
    // When executing tests
    return "unknown_version"
  }
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(inputFile)
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("Version to be used: ${version}")
  return version
}