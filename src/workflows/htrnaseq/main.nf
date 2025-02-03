workflow run_wf {
  take:
    input_ch

  main:
    // The featureData only has one requirement: the genome annotation.
    // It can be generated straight away.
    f_data_ch = input_ch
      | create_fdata.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: [
          "gtf": "annotation",
          "output": "f_data"
        ],
        toState: {id, result, state -> ["f_data": result.output]}
      )

    // Perform mapping of each well.
    mapping_ch = input_ch
      | well_demultiplex.run(
        fromState: [
            "input_r1": "input_r1",
            "input_r2": "input_r2",
            "barcodesFasta": "barcodesFasta",
        ],
        toState: [
          "input_r1": "output_r1",
          "input_r2": "output_r2",
        ]
      )
      | parallel_map_wf.run(
        fromState: {id, state ->
          [
            "input_r1": state.input_r1,
            "input_r2": state.input_r2,
            "barcodesFasta": state.barcodesFasta,
            "umi_length": state.umi_length,
            "output": state.star_output[0],
            "genomeDir": state.genomeDir,
          ]
        },
        toState: [
          "star_output": "output",
        ]
      )

    // From the mapped wells, create statistics based on the BAM file
    // and join the events back to pool level.
    pool_ch = mapping_ch
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

    // The well statistics are merged on pool level. 
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
    
    p_data_ch = star_logs_ch.join(pool_statistics_ch, remainder: true)
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

    eset_ch = p_data_ch.join(f_data_ch, remainder: true)
      | map {id, p_data_state, f_data_state ->
        def newState = p_data_state + ["f_data": f_data_state["f_data"]]
        [id, newState]
      }
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
      | map {ids_and_states ->
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

    output_ch = eset_ch.join(report_channel)
      | map {id, state_eset, state_report ->
        def new_state = state_eset + ["html_report": state_report.html_report]
        [id, new_state]
      }
      | setState([
        "star_output": "star_output", 
        "fastq_output_r1": "input_r1",
        "fastq_output_r2": "input_r2",
        "star_output": "star_output",
        "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChromPool",
        "star_qc_metrics": "star_qc_metrics",
        "eset": "eset",
        "f_data": "f_data",
        "p_data": "p_data",
        "html_report": "html_report",
      ])


  emit:
    output_ch
}
