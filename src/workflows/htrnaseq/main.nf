workflow run_wf {
  take:
    raw_ch

  main:
    input_ch = raw_ch
      // Use the event ID as the default for the sample ID
      | map {id, state ->
        def sample_id = state.sample_id ?: id 
        def newState = state + ["sample_id": sample_id, "run_id": id]
        return [id, newState]
      }

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

    // The featureData only has one requirement: the genome annotation.
    // It can be generated straight away. Most of the time, there is one shared 
    // annotation for all of the inputs and the fData should only be calculated once.
    // The state is manpulated in such a way that there is one event created per unique
    // input annotation file. In turn, the featureData file can joined into the original input
    // channel which allows it to be shared across events if required.
    f_data_ch = input_ch
      | toSortedList()
      | flatMap {ids_and_states ->
        def annotation_files = ids_and_states.inject([:]){ old_state, id_and_state ->
          def (id, state) = id_and_state
          def annotation_file = state.annotation
          def new_state = old_state + [(annotation_file): (old_state.getOrDefault(annotation_file, []) + [id])]
          return new_state
        }
        def file_names = annotation_files.keySet().collect{it.name}
        assert (file_names.toSet().size() == file_names.size()), 
          "Please make sure that the annotation files have unique file names."
        def new_states = annotation_files.collect{annotation_file, value ->
          def new_state = [annotation_file.name , ["annotation": annotation_file, "event_ids": value]]
          return new_state
        }
        return new_states 
      }
      | create_fdata.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: [
          "gtf": "annotation",
          "output": "f_data"
        ],
        toState: ["f_data": "output"]
      )
      | flatMap {_, state -> 
          def new_states = state.event_ids.collect{event_id ->
            [event_id, ["f_data": state.f_data]]
          }
          return new_states
      }

    // Perform mapping of each well.
    demultiplex_ch = input_ch
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


    concat_samples_ch = demultiplex_ch.join(f_data_ch, failOnMismatch: true, failOnDuplicate: true)
      | map {id, demultiplex_state, f_data_state ->
        def newState = demultiplex_state + ["f_data": f_data_state["f_data"]]
        [id, newState]
      }
      | concatRuns.run(
        fromState: [
          "input_r1": "input_r1",
          "input_r2": "input_r2",
          "sample_id": "sample_id",
        ],
        toState: {id, result, state ->
          def state_overwite = [
            "input_r1": result.output_r1,
            "input_r2": result.output_r2,
            "_meta": ["join_id": state.run_id]
          ]
          return state + state_overwite
        }
      )

    pool_ch = concat_samples_ch.join(fastq_output_directory_ch, failOnMismatch: true, failOnDuplicate: true)
      | map {id, concat_state, fastq_output_directory_state ->
        def new_state = concat_state + fastq_output_directory_state
        return [id, new_state]
      } 
      | parallel_map.run(
        directives: ["label": ["highmem", "highcpu"]],
        fromState: {id, state ->
          [
            "input_r1": state.input_r1,
            "input_r2": state.input_r2,
            "barcodesFasta": state.barcodesFasta,
            "umiLength": state.umi_length,
            "output": state.star_output[0],
            "genomeDir": state.genomeDir,
          ]
        },
        toState: [
          "star_output": "output",
        ]
      )
      // Split the events from 1 event per pool into events per well
      // and add extra metadata about the wells to the state.
      | well_metadata.run(
        fromState: [
          "barcodesFasta": "barcodesFasta",
          "input_r1": "input_r1",
          "input_r2": "input_r2",
          "star_mapping": "star_output"
        ],
        toState: [
          "input_r1": "output_r1",
          "input_r2": "output_r2",
          "pool": "pool",
          "well_id": "well_id",
          "barcode": "barcode",
          "lane": "lane",
          "n_wells": "n_wells",
          "star_mapping": "well_star_mapping",
        ]
      )
      // Use the bam file to generate statistics
      | generate_well_statistics.run(
        directives: [label: ["verylowmem", "verylowcpu"]],
        fromState: { id, state ->
          [
            "input": state.star_mapping.resolve('Aligned.sortedByCoord.out.bam'),
            "barcode": state.barcode,
            "well_id": state.well_id,
          ]
        },
        toState: [
          "nrReadsNrGenesPerChromWell": "nrReadsNrGenesPerChrom",
        ]
      )
      // Join the events back to pool-level
      | map {id, state ->
        // Create a special groupKey, such that groupTuple
        // knows when all the barcodes have been grouped into 1 event.
        // This way the processing is as distributed as possible.
        def key = groupKey(state.pool, state.n_wells)
        def newEvent = [key, state]
        return newEvent
      }
      // Use a custom sorting function because sort: 'hash'
      // requires a hash to be calculated on every entry of the state
      // This is inefficient when the number of events is large 
      // (i.e large number or barcodes).
      // Sorting on lexographical order of the barcode is sufficient here.
      | groupTuple(sort: {a, b -> a.barcode <=> b.barcode})
      | map {id, states ->
        // Gather the keys from all states. for some state items,
        // we need gather all the different items from across the states
        def barcodes = states.collect{it.barcode}
        assert barcodes.clone().unique().size() == barcodes.size(), \
          "Error when gathering information for pool ${id}, barcodes are not unique!"
        def well_ids = states.collect{it.well_id}
        assert well_ids.clone().unique().size() == well_ids.size(), \
          "Error when gathering information for pool ${id}, well IDs are not unique!"
        def custom_state = [
          "input_r1": states.collect{it.input_r1},
          "input_r2": states.collect{it.input_r2},
          "barcode": barcodes,
          "well_id": well_ids,
          "star_mapping": states.collect{it.star_mapping},
          // Well and pool stats should be carefully kept separate.
          // The workflow argument points to the name for the pool statistics:
          "nrReadsNrGenesPerChromWell": states.collect{it.nrReadsNrGenesPerChromWell},
          "nrReadsNrGenesPerChromPool": states[0].nrReadsNrGenesPerChrom
        ]
        //For many state items, the value is the same across states.
        def other_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
          }.minus(custom_state.keySet())
        // All other state should have a unique value
        def old_state_items = other_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across modalities. Please report this \
                                                 as a bug. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
          }

        def new_state = custom_state + old_state_items
        [id.getGroupTarget(), new_state]
      }

    pool_statistics_ch = pool_ch
      | generate_pool_statistics.run(
        directives: ["label": ["lowmem", "verylowcpu"]],
        fromState: [
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChromWell",
          "nrReadsNrGenesPerChromPool": "nrReadsNrGenesPerChromPool"
        ],
        toState: [
          "nrReadsNrGenesPerChromPool": "nrReadsNrGenesPerChromPool"
        ]
      )

    // The statistics from the STAR logs of different wells are joined
    // on pool level 
    star_logs_ch = pool_ch
      | combine_star_logs.run(
        directives: ["label": ["lowmem", "verylowcpu"]],
        fromState: {id, state -> [
            "star_logs": state.star_output.collect{it.resolve("Log.final.out")},
            "gene_summary_logs": state.star_output.collect{it.resolve("Solo.out/Gene/Summary.csv")},
            "reads_per_gene_logs": state.star_output.collect{it.resolve("ReadsPerGene.out.tab")},
            "barcodes": state.barcode,
            "output": state.star_qc_metrics
          ]
        },
        toState: [
          "star_qc_metrics": "output",
        ]
      )
    
    eset_ch = star_logs_ch.join(pool_statistics_ch, failOnMismatch: true, failOnDuplicate: true)
      | map {id, star_logs_state, pool_statistics_state ->
        def newState = star_logs_state + ["nrReadsNrGenesPerChromPool": pool_statistics_state.nrReadsNrGenesPerChromPool]
        return [id, newState]
      }
      | create_pdata.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: [
          "star_stats_file": "star_qc_metrics",
          "nrReadsNrGenesPerChromPool": "nrReadsNrGenesPerChromPool",
          "output": "p_data"
        ],
        toState: ["p_data": "output"],
      )
      | create_eset.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: [
          "pDataFile": "p_data",
          "fDataFile": "f_data",
          "mappingDir": "star_output",
          "output": "eset",
          "barcodes": "barcode",
          "poolName": "pool",
        ],
        toState: [
          "eset": "output",
        ]
      )

    report_channel = eset_ch
      | toSortedList()
      | filter{it}
      | map {ids_and_states ->
        assert ids_and_states.size() > 0 
        def states = ids_and_states.collect{it[1]}
        def html_report = states[0].html_report
        def ids = ids_and_states.collect{it[0]}
        def esets = states.collect{it.eset}
        ["report", ["esets": esets, "html_report": html_report, "original_ids": ids]]
      }
      | create_report.run(
        fromState: [
          "eset": "esets",
          "output_report": "html_report",
        ],
        toState: [
          "html_report": "output_report"
        ]
      )
      | flatMap {id, state ->
        state.original_ids.collect{original_id ->
          [original_id, ["html_report": state.html_report]]
        }
      }

    output_ch = eset_ch.join(report_channel, failOnMismatch: true, failOnDuplicate: true)
      | map {id, state_eset, state_report ->
        def new_state = state_eset + [
          "html_report": state_report.html_report,
        ]
        [id, new_state]
      }
      | setState([
        "star_output": "star_output",
        "fastq_output": "fastq_output_directory",
        "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChromPool",
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
