workflow run_wf {
  take:
    input_ch

  main:
    mapping_ch = input_ch
      | map {id, state ->
        def newState = state + ["n_barcodes": state.barcodesFasta.countFasta() as int]
        [id, newState]
      }
      | well_demultiplex.run(
        fromState: [
            "input_r1": "input_r1",
            "input_r2": "input_r2",
            "barcodesFasta": "barcodesFasta",
        ],
        toState: { id, result, state ->
          def newState = state + result + [
              "fastq_output_r1": result.output_r1, 
              "fastq_output_r2": result.output_r2, 
              "input_r1": result.output_r1,
              "input_r2": result.output_r2,
            ]
          return newState
        }
      )
      | parallel_map_wf.run(
        fromState: {id, state ->
          [
            "input_r1": state.input_r1[0],
            "input_r2": state.input_r2[0],
            "barcode": state.barcode,
            "pool": state.pool,
            "output": state.star_output[0],
            "genomeDir": state.genomeDir,
          ]
        },
        toState: ["star_output": "output"]
      )
      | generate_well_statistics.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: { id, state ->
          [
            "input": state.star_output.resolve('Aligned.sortedByCoord.out.bam'),
            "barcode": state.barcode,
          ]
        },
        toState: [
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom",
          "nrReadsNrUMIsPerCB": "nrReadsNrUMIsPerCB",
        ]
      )
      | map {id, state ->
        // Create a special groupKey, such that groupTuple
        // knows when all the barcodes have been grouped into 1 event.
        // This way the processing is as distributed as possible.
        def key = groupKey(state.pool, state.n_barcodes)
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
        def collected_state = [
          "fastq_output_r1": states.collect{it.fastq_output_r1[0]},
          "fastq_output_r2": states.collect{it.fastq_output_r2[0]},
          "annotation": states.collect{it.annotation}[0],
          "barcodes": states.collect{it.barcode},
          "star_output": states.collect{it.star_output},
          "nrReadsNrGenesPerChrom": states.collect{it.nrReadsNrGenesPerChrom},
          "star_qc_metrics": states.collect{it.star_qc_metrics}[0],
          "eset": states.collect{it.eset}[0],
          "pool": states.collect{it.pool}[0],
        ]
        [id.getGroupTarget(), collected_state]
      }
      | generate_pool_statistics.run(
        directives: ["label": ["lowmem", "lowcpu"]],
        fromState: [
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom",
        ],
        toState: [
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChromPool"
        ]
      )
      | map {id, state ->
        def extraState = [
            "star_logs": state.star_output.collect{it.resolve("Log.final.out")},
            "summary_logs": state.star_output.collect{it.resolve("Solo.out/Gene/Summary.csv")},
            "reads_per_gene": state.star_output.collect{it.resolve("ReadsPerGene.out.tab")}
        ]
        def newState = state + extraState
        return [id, newState]
      }
      | combine_star_logs.run(
        directives: ["label": ["lowmem", "lowcpu"]],
        fromState: [
          "star_logs": "star_logs",
          "gene_summary_logs": "summary_logs",
          "reads_per_gene_logs": "reads_per_gene",
          "barcodes": "barcodes",
          "output": "star_qc_metrics",
        ],
        toState: [
          "star_qc_metrics": "output",
        ]
      )


    f_data_ch = mapping_ch
      | create_fdata.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: ["gtf": "annotation"],
        toState: ["f_data": "output"],
      )

    p_data_ch = mapping_ch
      | create_pdata.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: [
          "star_stats_file": "star_qc_metrics",
          "nrReadsNrGenesPerChromPool": "nrReadsNrGenesPerChrom",
        ],
        toState: ["p_data": "output"],
      )

    output_ch = f_data_ch.join(p_data_ch, remainder: true)
      | map {id, f_data_state, p_data_state ->
        def newState = f_data_state + ["p_data": p_data_state["p_data"]]
        [id, newState]
      }
      | create_eset.run(
        directives: [label: ["lowmem", "lowcpu"]],
        fromState: [
          "pDataFile": "p_data",
          "fDataFile": "f_data",
          "mappingDir": "star_output",
          "output": "eset",
          "barcodes": "barcodes",
          "poolName": "pool",
        ],
        toState: [
          "eset": "output",
        ]
      )
      | setState([
        "star_output", 
        "fastq_output_r1",
        "fastq_output_r2",
        "star_output",
        "nrReadsNrGenesPerChrom",
        "star_qc_metrics",
        "eset",
      ])
      | niceView()


  emit:
    output_ch
}
