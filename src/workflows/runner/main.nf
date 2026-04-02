import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.atomic.AtomicBoolean

def date = new Date().format('yyyyMMdd_hhmmss')

session = nextflow.Nextflow.getSession()
final service = session.publishDirExecutorService()



def viash_config = java.nio.file.Paths.get("${moduleDir}/_viash.yaml")
def version = get_version(viash_config)

if (!params.containsKey("results_publish_dir") || !params.containsKey("fastq_publish_dir")) {
   error "Please set both 'results_publish_dir' and 'fastq_publish_dir'."
}

// S3 paths containing double slashes might cause issues with empty objects being created
// Remove trailing slashes from the publish dir. The params map is immutable, so create a copy
def regex_trailing_slashes = ~/\/+$/
def results_publish_dir = params.results_publish_dir - regex_trailing_slashes
def fastq_publish_dir = params.fastq_publish_dir - regex_trailing_slashes


/* 
This is a utility workflow that gathers the input events and saves their state to a YAML file.
*/
workflow save_params_wf {
  take:
  input_ch

  main:
     
    output_ch = input_ch
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

  emit:
  output_ch
}


workflow run_wf {
  take:
    raw_ch

  main:
    input_ch = raw_ch
      // Use the ID of the event as a default for the `run_id`
      | map {id, state ->
        def new_state = state
        if (!state.run_id) {
          new_state = new_state + ["run_id": id]
        }
        return [id, new_state]
      }

    // Save the input parameters to a file
    // Adds a 'run_params' key to the state with the file.
    save_params_ch = input_ch
      | save_params_wf
  
    /*
      WELL DEMULTIPLEXING

      The provided input is a single directory per input event. However, the well demultiplexing
      requires a list of FASTQ files to be demultiplexed. The `listInputDir` workflow lists the 
      content of the input directory and deduces the sample IDs from the file name.

      Multiple experiments can be provided, meaning that a input events represent sequencing runs
      but provided per experiment (i.e. a sequencing run - experiment combination is a unique event).
      This is because a single sequencing run can contribute to multiple experiments. However, 
      the demultiplexing of the plate FASTQ files into well FASTQ files should ideally only be done once
      per sequencing run and not not per-experiment. The `run_id` is used as a way to group the sequencing
      runs and to determine the name of the output folder. 
    */
    demultiplex_ch = input_ch
      | map {id, state ->
        [state.run_id, state + ["event_id": id]]
      }
      | groupTuple(by: 0, sort: 'hash')
      // Here there might be multiple states if there are multiple experiments using the same sequencing run
      | map {run_id, states ->
        // The 'run_id' parameter _must_ map one to one to an input folder.
        def input_dirs = states.collect{it.input}
        assert input_dirs.collect{it.toUriString()}.unique().size() == 1, \
          "Found more than one input directory that corresponds with run ID '${run_id}'."

        // Its not possible to demultiplex the same sequencing run with different barcodes fasta
        def barcode_fastas = states.collect{it.barcodesFasta}
        assert barcode_fastas.collect{it.toUriString()}.unique().size() == 1, \
          "Demultiplexing the same run (${run_id}) with two different barcode layouts it not supported."

        // Its not allowed to use the same sequencing run twice or more for the same experiment
        def projects_experiments = states.collect{[it.project_id, it.experiment_id]}
        def project_experiments_str = projects_experiments.collect{project, experiment -> "${project}/${experiment}"}
        assert project_experiments_str.unique().size() == project_experiments_str.size(), \
          "It is not possible to use a sequencing run (${run_id}) for the same experiment multiple times. Found: ${project_experiments_str}"

        def new_state = [
          "run_id": run_id,
          "input": input_dirs[0],
          "barcodesFasta": barcode_fastas[0], // This is needed for the well demultiplexing
          "event_ids": states.collect{it.event_id},
          "project_experiment_ids": projects_experiments
        ]
        [run_id, new_state]
      } 
      /* ListInputDir puts the sample_id as the event ID (slot 0 from the tuple).
         If the FASTQ files are split per lane; the ListInputDir subworkflow will output the fastq files
         as a list in 'r1_output' and 'r2_output'. There are three things to take into account:
          (1) an input folder can contain input files from multiple samples (pools) - this is highly probable.
          (2) there might be multiple FASTQs for a single sample that correspond to the lanes - they originate from the same input directory.
          (3) there might be FASTQ files for the same sample across directories.
              This is a concatenation scenario: the sample was split across sequencing runs in order to increase sequencing depth.
                
         Item (2) -the lanes- are handled by providing a list of FASTQ files to the `well_demultiplexing` workflow.
         However, (3) is handled by the `htrnaseq` workflow. This means that each sequencing run is demultiplexed separately, 
         and the well FASTQs concatenated across sequencing runs instead of concatenating on plate level.
      */
      | listInputDir.run(
        fromState: [
          "input": "input"
        ],
        toState: [
            "output_r1": "r1_output",
            "output_r2": "r2_output",
            "sample_id": "sample_id"
          ]
      )
      /* Samples might be split across sequencing runs — this is item (3) from the comment above. 
         The listInputDir subworkflow uses just the sample ID for the 
         event. Even though the concatenation of samples is not handled by with well_demultiplexing workflow,
         there might be more than 1 event with the same ID after using the listInputDir subworkflow. 
         In order to make them unique again we just add the run ID, otherwise events with duplicate IDs
         are dropped.
      */
      | map {id, state -> ["${state.run_id}/${id}".toString(), state]}
      // Each plate (sample) is demultiplexed separately
      | well_demultiplex.run(
        fromState: [
            "input_r1": "output_r1",
            "input_r2": "output_r2",
            "barcodesFasta": "barcodesFasta"
        ],
        toState: {id, result, state ->
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
    
    /*
    Publish the results from the demultiplexing. The publishing is done for each sequencing run.
    Because the well demultiplexing is only done once per sequencing run anyway there is not
    much need for state manipulation here.
    */
    fastq_publish_ch = demultiplex_ch
      | map {id, state ->
        def new_state = state + [
          "fastq_prefix": "${state.run_id}/${date}_htrnaseq_${version}/${state.sample_id}".toString()
        ]
        [id, new_state]
      }
      | publish_fastqs.run(
        fromState: {id, state -> [
            "input": state.output_r1 + state.output_r2,
            "output": state.fastq_prefix,
          ]
        },
        toState: { id, result, state -> state },
        directives: [
          publishDir: [
            path: fastq_publish_dir, 
            overwrite: false,
            mode: "copy"
          ]
        ]
      )

    /*
    CREATE COUNT MATRICES

    The results from the well demultiplexing are the plates for each unique sequencing run.
    The information for the mapping (creating the count matrices) needs to be added from the input state.
    However, the input state provides events on sequencing runs per experiment, which means
    that the event for a single sequencing run is duplicated for how many experiments it will contribute to.
    In order to facilitate this join, the input event IDs were stored in the output state from
    the well demultiplexing. After the join the events are on the level of the input events, meaning
    that they are provided per experiment.


    After the join, there might be multiple events with the same sample ID: a single well plate 
    might have been put on different sequencing runs in order to increase sequencing depth.
    This means that they have the same 'sample_id' state key. If the sample ID would be used as event ID, 
    these duplicate events would be dropped. The event ID from the input channel is added as a prefix to
    the sample ID in order to make sure they are unique. In order to indicate to the subworkflow that 
    these samples need to be concatenated, value for the `sample_id` argument can be provided.
    This will cause the the `well_fastqs_to_esets` workflow to do the grouping and concatenation of the 
    events with the same sample ID.
    */
    demultiplex_with_input_ids_ch = demultiplex_ch
      // The IDs in the demultiplex_ch are of format '<run_id>/<sample_id>' (not split per experiment.)
      | flatMap {id, state ->
        state.event_ids.collect{ event_id ->
          def new_state = state.findAll{k, v -> !["event_ids", "project_experiment_ids"].contains(k)} + ["event_id": event_id]
          [event_id, new_state]
        }
      }
      // The event IDs at this point are the same IDs in the `input_ch` in order to do the join

    htrnaseq_ch = input_ch
      | cross(demultiplex_with_input_ids_ch)
      | map {input_event, demux_event  ->
        def (event_id, demux_state) = demux_event
        def input_state = input_event[1] 
        def keys_to_transfer = [
          "umi_length",
          "annotation",
          "genomeDir",
          "pools",
          "experiment_id",
          "project_id",
          "star_output_dir",
          "nrReadsNrGenesPerChrom_dir",
          "star_qc_metrics_dir",
          "eset_dir",
          "f_data_dir",
          "p_data_dir",
          "run_params"
        ]
        def new_state = demux_state + input_state.subMap(keys_to_transfer)

        def keys_to_rename = [
          "star_output_dir",
          "nrReadsNrGenesPerChrom_dir",
          "star_qc_metrics_dir",
          "eset_dir",
          "f_data_dir",
          "p_data_dir",
          "run_params"
        ]
        new_state = new_state.collectEntries{k, v ->
          def newKey = keys_to_rename.contains(k) ? "${k}_workflow" : k
          [newKey, v]
        }
        // Using 'event_id' ensures that the event IDs are unique.
        // The concatenation of samples is handled by providing the 'sample_id'
        // argument to the `well_fastqs_to_esets` workflow. 
        def new_id = "${event_id}/${demux_state.sample_id}"
        [new_id, new_state]
      }
      | filter {id, state ->
        def test_val = (state.pools == null || state.pools.isEmpty() || state.pools.contains(state.sample_id))
        if (!test_val) {
          log.info("Filtering out ${id} because it is not in the selected pools for this experiment (${state.pools}).")
        }
        test_val
      }
      | well_fastqs_to_esets.run(
          args: [
            f_data: 'fData/$id.txt',
            p_data: 'pData/$id.txt',
            star_output: 'star_output/$id/*',
            eset: 'esets/$id.rds',
            nrReadsNrGenesPerChrom: 'nrReadsNrGenesPerChrom/$id.txt',
            star_qc_metrics: 'starLogs/$id.txt',
            html_report: "report.html"
        ],
        fromState: {id, state -> [
            "input_r1": state.output_r1,
            "input_r2": state.output_r2,
            "annotation": state.annotation,
            "umi_length": state.umi_length,
            "barcodesFasta": state.barcodesFasta,
            "genomeDir": state.genomeDir,
            // Concatenate the samples with the same sample ID, but not across experiments.
            "sample_id": "${state.project_id}/${state.experiment_id}/${state.sample_id}",
          ]
        },
        toState: { id, result, state -> state + result.findAll{it.key != "run_params"} }
      )


    /* The well demultiplexing provides a list of FASTQ files as output 
       while this runner workflow needs a directory.
       Furthermore, the demultiplexing output is provided per unique sequencing run 
       (even when used by multiple experiments) but this runner workflow outputs
       the events per experiment. Because multiple sequencing runs can contribute to more
       than one experiment, a list of directories is output by the runner. Here the output
       from the demultiplexing is transformed in order to add it to the eset (results_publish_ch) output.

       Normally, one would expect to use fromState/toState logic for this with the
       `well_fastqs_to_esets` workflow. However, concatenation of samples is done within this workflow
       and it does not output the FASTQs that contributed to a certain sample.
    */
    fastq_output_directory_ch = demultiplex_ch
      | flatMap {id, state ->
        def base_state = state.findAll{k, v -> k != "project_experiment_ids"}
        state.project_experiment_ids.collect { project_id, experiment_id ->
          // Note: the format of the ID here must match with the format chosen for the 'sample_id' argument
          // in the well_fastqs_to_esets workflow in order to do the join with its output.
          ["${project_id}/${experiment_id}/${state.sample_id}".toString(), base_state + ["project_id": project_id, "experiment_id": experiment_id]]
        }
      }
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def fastq_output_dirs = states.collect{it.fastq_output_directory}
        def new_state = ["fastq_output_directory": fastq_output_dirs]
        def new_event = [id, new_state]
        return [id, new_state]
      }

    // While fastq_output_directory_ch contains all of the events,
    // the htrnaseq_ch may have been filtered to remove samples based on the 'pools' argument
    // This is why mismatches are allowed.
    results_publish_ch = htrnaseq_ch.join(fastq_output_directory_ch, failOnDuplicate: true, failOnMismatch: false)
      // Add the FASTQ directories from the demultiplexing output
      | map {id, htrnaseq_state, fastq_output_directory_state ->
        def new_state = htrnaseq_state + fastq_output_directory_state
        return [id, new_state]
      } 
      | combine(save_params_ch)
      // Add the run parameter YAML to the output
      | map {id, grouped_ch_state, _, save_params_state ->
        assert save_params_state.run_params.isFile()
        def new_state = grouped_ch_state + ["run_params": save_params_state.run_params]
        return [id, new_state]
      }
      // Group the events in order to publish the results per experiment
      | map {id, state -> 
        def new_id = "${state.project_id}/${state.experiment_id}".toString()
        [new_id, state]
      }
      | groupTuple(by: 0, sort: 'hash')
      // Join the results from the different plates in this experiment
      | map{ id, states ->
          // The STAR output is a directory for each well in a plate (or pool of plates).
          // The wells are grouped into a directory per sample. The name of this directory should
          // match the sample_id.
          def star_output_samples = states.collectMany{state -> 
            state.star_output.collect{
              def star_sample_dir = it.parent
              assert star_sample_dir.name == state.sample_id: "Unexpected state: the parent directory of STAR output \
                path '${it}' should match with the sample ID ${sample_id}"
              star_sample_dir
            }.unique()
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
            "run_params_workflow",
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
            "event_id"
          ]
          def state_collect = state_keys_collect.collectEntries{ key_ ->
            [key_, states.collect{it.get(key_)}]
          }

          new_state["results_prefix"] = "${state_unique_keys.project_id}/${state_unique_keys.experiment_id}/data_processed/${date}_htrnaseq_${version}"
          new_state = new_state + state_unique_keys + state_collect
          [id, new_state]  
      }
      | publish_results.run(
        fromState: { id, state ->

          println("Publising results to ${results_publish_dir}/${state.results_prefix}")

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
            html_report_output: "${state.results_prefix}/${state.html_report.name}", 
            run_params_output: "${state.results_prefix}/${state.run_params_workflow}",
            star_output_dir: "${state.results_prefix}/${state.star_output_dir_workflow}",
            nrReadsNrGenesPerChrom_dir: "${state.results_prefix}/${state.nrReadsNrGenesPerChrom_dir_workflow}",
            star_qc_metrics_dir: "${state.results_prefix}/${state.star_qc_metrics_dir_workflow}",
            eset_dir: "${state.results_prefix}/${state.eset_dir_workflow}",
            f_data_dir: "${state.results_prefix}/${state.f_data_dir_workflow}",
            p_data_dir: "${state.results_prefix}/${state.p_data_dir_workflow}"
          ]
        },
        toState: { id, result, state -> 
          result + [
            "run_params": state.run_params,
            "_meta": ["join_id": state.event_id[0]],
            "results_prefix": state.results_prefix
          ] 
        },
        directives: [
          publishDir: [
            path: results_publish_dir, 
            overwrite: false,
            mode: "copy"
          ]
        ]
      )



  has_published = new AtomicBoolean(false)

  interval_ch = channel.interval('30s'){ i ->
    // Allow this channel to stop generating events based on a later signal
    if (has_published.get()) {
      return channel.STOP
    }
    i
  }
  
  // Make sure that there is at least one event in the interval channel
  interval_at_least_one_event_ch = Channel.fromList([0]).concat(interval_ch)

  awaited_events_ch = results_publish_ch.mix(fastq_publish_ch)
    | toSortedList()
    | map {states ->
      if (states.size() == 0) {
        has_published.compareAndSet(false, true)
        error("There seems to be nothing to publish!")
      }
      states
    }

  await_ch = awaited_events_ch
    // Wait for processing events to be done
    // Create periodic events in order to check for the publishing to be done
    | combine(interval_at_least_one_event_ch)
    | until { event ->
      // Prevent until to output nothing by stopping on the first item of the channel.
      // It will output 'null' when its the first iteration.
      // This happens when there is not a lot of data to publish and/or the transfer is fast.
      if (event[-1] == 0) {
        return false
      }
      println("Checking if publishing has finished in service ${service}")
      def running_tasks = null
      if(service instanceof ThreadPoolExecutor) {
        def completed_tasks = service.getCompletedTaskCount()
        def task_count = service.getTaskCount()
        if (completed_tasks > 0) {
          running_tasks = completed_tasks - task_count
        }
      }
      else if( System.getenv('NXF_ENABLE_VIRTUAL_THREADS') ) {
        running_tasks = service.threadCount()
      }
      else {
        error("Publishing service of class ${service.getClass()} is not supported.")
      }
      
      if (running_tasks == 0) {
        println("Publishing has finished all current tasks. Continuing execution.")
        return true
      }
      println("Workflow is publishing. Waiting...")
      return false
    }
    | last()
    | map{ events ->
        // Signal to interval channel to stop generating events.
        has_published.compareAndSet(false, true)
        // Remove the value that was added by the interval channel.
        return events.dropRight(1)
    }
    | flatMap { events ->
        println("Creating transfer_complete.txt files.")
        result_events = []
        events.each {id, state -> 
          if (state.containsKey("fastq_prefix")) {
            def complete_file_fastqs = file("${fastq_publish_dir}/${state.fastq_prefix}/transfer_completed.txt")
            complete_file_fastqs.text = "" // This will create a file when it does not exist.
          } else if (state.containsKey("results_prefix")) {
            def complete_file_results = file("${results_publish_dir}/${state.results_prefix}/transfer_completed.txt") 
            complete_file_results.text = ""
            result_events.add([id, state])
          } else {
            error "State should contain either 'fastq_prefix' or 'results_prefix'"
          }
        }
        result_events
    }
    | setState([
        "star_output_dir",
        "nrReadsNrGenesPerChrom_dir",
        "star_qc_metrics_dir",
        "eset_dir",
        "f_data_dir",
        "p_data_dir",
        "run_params",
        "_meta"
      ]
    )

  emit:
    await_ch

}

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("HT-RNAseq version to be used: ${version}")
  return version
}