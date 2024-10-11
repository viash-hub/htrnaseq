workflow run_wf {
  take:
    input_ch

  main:
    mapping_ch = input_ch
      | well_demultiplex.run(
        fromState: { id, state ->
          [
            input_r1: state.input_r1,
            input_r2: state.input_r2,
            barcodesFasta: state.barcodesFasta,
          ]
        },
        toState: { id, result, state ->
          state + result + [
              fastq_output_r1: result.output_r1, 
              fastq_output_r2: result.output_r2, 
              input_r1: result.output_r1,
              input_r2: result.output_r2,
            ]
        },
        directives: [label: ["midmem", "midcpu"]]
      )
      | parallel_map_wf.run(
        fromState: {id, state ->
          def star_output = state.star_output[0]
          [
            "input_r1": state.input_r1[0],
            "input_r2": state.input_r2[0],
            "barcode": state.barcode,
            "pool": state.pool,
            "output": state.star_output[0],
            "genomeDir": state.genomeDir,
          ]
        },
        toState: {id, result, state -> 
          state + ["star_output": result.output,]
        },
      )
      | generate_well_statistics.run(
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
        [state.pool, id, state]
      }
      | groupTuple(by: 0, sort: "hash")
      | map {id, well_ids, states ->
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
        [id, collected_state]
      }
      | generate_pool_statistics.run(
        fromState: [
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom",
        ],
        toState: {id, result, state -> 
          state + ["nrReadsNrGenesPerChrom": result.nrReadsNrGenesPerChromPool]
        }
      )
      | map {id, state ->
        def extraState = [
            "star_logs": state.star_output.collect{it.resolve("Log.final.out")},
            "summary_logs": state.star_output.collect{it.resolve("Solo.out/Gene/Summary.csv")},
            "reads_per_gene": state.star_output.collect{it.resolve("ReadsPerGene.out.tab")}
        ]
        return [id, state + extraState]
      }
      | combine_star_logs.run(
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
        fromState: ["gtf": "annotation"],
        toState: ["f_data": "output"],
      )

    p_data_ch = mapping_ch
      | create_pdata.run(
        fromState: [
          "star_stats_file": "star_qc_metrics",
          "nrReadsNrGenesPerChromPool": "nrReadsNrGenesPerChrom",
        ],
        toState: ["p_data": "output"],
      )

    output_ch = f_data_ch.join(p_data_ch, remainder: true)
      | map {id, f_data_state, p_data_state ->
        [id, f_data_state + ["p_data": p_data_state["p_data"]]]
      }
      | create_eset.run(
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
        "f_data",
        "p_data"
      ])
      | niceView()


  emit:
    output_ch
}
