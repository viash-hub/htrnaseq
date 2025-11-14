workflow run_wf {
  take:
    raw_ch

  main:
    demultiplex_ch = raw_ch
      // Use the event ID as the default for the sample ID
      | map {id, state ->
        def sample_id = state.sample_id ?: id 
        def newState = state + ["sample_id": sample_id]
        return [id, newState]
      }
      /* 
         Save the input parameters to a YAML file.
        `Path` objects need to be cast to strings using `toUriString`,
         otherwise prefixes like `s3://` might be missing.
         A string representation for the YAML is generated here, 
         encoded as a base64 string and provided to the save_params component
         in order to save it to a file. 
      */
      | save_params.run(
        runIf: { id, state ->
          state.run_params != null
        },
        fromState: {id, state ->
          // Define the function before using it
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
          def convertedState = state.collectEntries { k, v -> [(k): convertPaths(v)] }
          
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
        The provided FASTQ files correspond to the contents of a well plate.
        Each well was tagged with barcode which resides at the beginning of the sequences.
        Well demultiplexing will split a FASTQ file for a plate into multiple FASTQ files,
        one for each well. `input_r1` and `input_r2` are lists: corresponding pairs from the list
        will be demultiplexed separately before joining the results. This is usefull when dealing 
        with multiple lanes. However, this is not yet joining the results based on identical sample name
        when concatenation FASTQ files across sequencing runs (e.g. to increase sequencing depth).

        When concatenating lanes, a joining processes is started for each FASTQ file. Each of these
        processes has their unique working directory. When necessary the `well_demultiplexing` workflow
        will move those files into one directory again. This is required because the `htrnaseq` workflow
        outputs directories instead of the individual files. The fact that all FASTQ files are in the 
        same directories is also checked by this workflow as a fail-safe.

        The input events are provided on plate level and the output is also provided on plate level.
      */
      | well_demultiplex.run(
        fromState: [
            "input_r1": "input_r1",
            "input_r2": "input_r2",
            "barcodesFasta": "barcodesFasta",
        ],
        toState: {id, result, state ->
          def all_fastq = result.output_r1 + result.output_r2
          def output_dir = all_fastq.collect{it.parent}.unique()
          assert output_dir.size() == 1, "Expected output from well demultiplexing (id $id) to reside into one directory. Found: $output_dir"
          def new_state = state + [
            "input_r1": result.output_r1,
            "input_r2": result.output_r2,
            "fastq_output_directory": output_dir[0],
          ]
          return new_state
        }
      )
    
    /* 
    The output FASTQ directories are added to the state by the `well_demultiplex` workflow.
    However, the well_demultiplex workflow outputs events on plate level per sequencing run, while the 
    `htrnaseq` workflow might output concatenated events: results across sequencing runs are joined.
    Here, we ascertain which FASTQ files have contributed to which result: when concatenation 
    was performed the list will contain multiple directories of FASTQ files. As a consequence,
    a FASTQ directory could be provided as output in more than 1 output event from the `htrnaseq` workflow.
    */ 
    fastq_output_directory_ch = demultiplex_ch
      | map {id, state ->
        def new_event = [state.sample_id, state]
        return new_event
      }
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def fastq_output_dirs = states.collect{it.fastq_output_directory}
        def new_state = ["fastq_output_directory": fastq_output_dirs]
        def new_event = [id, new_state]
        return [id, new_state]
      }

    /*
      Map a list of FASTQ files (one for each well) to a reference genome and generate count matrices.
      Sometimes counts from different FASTQ files need to be concatenated. This is done bases on the sample_id:
      if the sample ID of the two plates are identical, the FASTQ files will we joined _before_ mapping.
    */
    htrnaseq_ch = demultiplex_ch
      /* 
        Because the well_fastqs_to_esets workflow does concatenation of results
        across sequencing runs, the events are different from the input events.
        Here we tell Viash how to handle the join between the old events and the new events.
      */
      | map {id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        return [id, new_state]
      }
      | well_fastqs_to_esets.run(
        fromState: [
          // Actual input
          "input_r1": "input_r1",
          "input_r2": "input_r2",
          "annotation": "annotation",
          "umi_length": "umi_length",
          "barcodesFasta": "barcodesFasta",
          "genomeDir": "genomeDir",
          "sample_id": "sample_id",
          // Arguments which determine the output file names
          "star_output": "star_output",
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom",
          "star_qc_metrics": "star_qc_metric",
          "eset": "eset",
          "f_data": "f_data",
          "p_data": "p_data",
          "html_report": "html_report"
        ],
        toState: {id, result, state -> state + result}
      )

    output_ch = htrnaseq_ch.join(fastq_output_directory_ch, failOnMismatch: true, failOnDuplicate: true)
      | map {id, htrnaseq_state, fastq_output_directory_state ->
        def new_state = htrnaseq_state + fastq_output_directory_state
        return [id, new_state]
      } 
      | setState([
        "star_output": "star_output",
        "fastq_output": "fastq_output_directory",
        "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom",
        "star_qc_metrics": "star_qc_metrics",
        "eset": "eset",
        "f_data": "f_data",
        "p_data": "p_data",
        "html_report": "html_report",
        "run_params": "run_params",
        "_meta": "_meta",
      ])

  emit:
    output_ch
}
