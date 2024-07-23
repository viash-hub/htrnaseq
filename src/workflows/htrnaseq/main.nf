workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      //| niceView()
      | well_demultiplex.run(
        fromState: { id, state ->
          [
            input_r1: state.input_r1,
            input_r2: state.input_r2,
            barcodesFasta: state.barcodesFasta,
          ]
        },
        toState: { id, result, state ->
          state + [
              output: result.output,
            ]
        },
        directives: [label: ["midmem", "midcpu"]]
      )

      | niceView()

      // TODO: turn this into a 'component': splitBarcodes
      | flatMap{ id, state ->
        state.output.collect{ p ->
          barcode = (p =~ /.*\\/([ACTG]*|unknown)_R?.*/)[0][1]
          pair_end = (p =~ /.*_(R[12])_.*/)[0][1]
          lane = (id =~ /.*_(L\d+).*/) ? (id =~ /.*_(L\d+).*/)[0][1] : "no_lanes"
          [ id + "__" + barcode + "__" + pair_end, state + [ sample_id: id, barcode: barcode, barcode_path: p, lane: lane, pair_end: pair_end ] ]
        }
      }

      | groupPairs.run(
        fromState: { id, state ->
          [
            input: state.barcode_path
          ]
        },
        toState: { id, result, state ->
          state + result
        }
      )

      // Does the sequencing platform use lanes?
      // Should those lanes be discriminated over?
      //
      | groupLanes.run(
        fromState: { id, state ->
          [
            input_r1: state.r1,
            input_r2: state.r2
          ]
        },
        toState: { id, result, state ->
          state + result + [ multiple_lanes: result.output_r1.size() > 1]
        }
      )

      | setState(
        [
          "input_r1": "r1",
          "input_r2": "r2",
          "barcodesFasta": "barcodesFasta",
          "genomeDir": "genomeDir",
          "barcode": "barcode",
        ]
      )

      | concat_text.run(
        key: "concat_text_r1",
        fromState: { id, state ->
          [
            input: state.input_r1,
            gzip_output: true
          ]
        },
        toState: { id, result, state ->
          state + [ input_r1: result.output ]
        }
      )
      | concat_text.run(
        key: "concat_text_r2",
        fromState: { id, state ->
          [
            input: state.input_r2,
            gzip_output: true
          ]
        },
        toState: { id, result, state ->
          state + [ input_r2: result.output ]
        }
      )

      | niceView()

      // Now we need a `groupBarcodes` again in order to run the mapping
      // without the need for an explicit list of barcodes

      // Prep the handover to the mapping step
      //| map { id, state ->
      //  // Group by the 'sample' id (ie barcode after demultiplexing)
      //  // Omit the 'unknown' group
      //  newState = state
      //    .output
      //    .collect{ p -> [ (p =~ /.*\\/([ACTG]*|unknown)_R?.*/)[0][1], p ] }
      //    .groupBy{ k, v -> k }
      //    .collect{ k, v -> [ k, v.collect{it[1]} ]}
      //    .findAll{ k, v -> k != "unknown"}
      //    .collect{ k, v -> [ barcode: k, fastq: v ] }
      //  [
      //    id, 
      //    state + [
      //      output: newState,
      //      barcodes: newState.collect{it.barcode},
      //      all: newState.collect{it.fastq}.flatten(),
      //      out: newState.collect{it.fastq}.flatten().flatten() // to keep Viash IO happy
      //    ]
      //  ]
      //}

      //| parallel_map.run(
      //  fromState: { id, state ->
      //    [
      //      poolName: id,
      //      input: state.all,
      //      genomeDir: state.genomeDir,
      //      barcodes: state.barcodes.join(",").toString(),
      //      wellBarcodesLength: 10,
      //      umiLength: 10,
      //    ]
      //  },
      //  toState: { id, result, state ->
      //    [
      //      output: result.output,
      //    ]
      //  }
      //)
      //
      //| niceView()
      //
      //| setState( [ "output": "out" ] )

  emit:
    output_ch
}
