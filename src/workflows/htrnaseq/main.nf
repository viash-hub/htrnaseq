workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
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

      // TODO: Expand this into matching a whitelist/blacklist of barcodes
      // ... and turn into separate component
      | filter{ id, state -> state.barcode != "unknown" }
      | concat_text.run(
        key: "concat_txt_r1",
        runIf: {id, state -> state.input_r1.size() > 1},
        fromState: { id, state ->
          [
            input: state.input_r1,
            gzip_output: true,
          ]
        },
        toState: { id, result, state ->
          state + [ input_r1: [ result.output ] ]
        }
      )
      | concat_text.run(
        key: "concat_text_r2",
        runIf: {id, state -> state.input_r2.size() > 1},
        fromState: { id, state ->
          [
            input: state.input_r2,
            gzip_output: true
          ]
        },
        toState: { id, result, state ->
          state + [ input_r2: [ result.output ] ]
        }
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
          state + ["star_output": result.output]
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
          "nrReadsNrGenesPerChrom": states.collect{it.nrReadsNrGenesPerChrom},
        ]
        def newState = states[0] + collected_state
        [id, newState]
      }
      | generate_pool_statistics.run(
        fromState: [
          "nrReadsNrGenesPerChrom": "nrReadsNrGenesPerChrom",
        ],
        toState: {id, result, state -> [
          state + ["nrReadsNrGenesPerChrom": result.nrReadsNrGenesPerChrom], 
        ]},
      )
      | setState([
        "star_output", 
        "fastq_output_r1",
        "fastq_output_r2",
        "star_output",
        "nrReadsNrGenesPerChrom",
      ])

  emit:
    output_ch
}
