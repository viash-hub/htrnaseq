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

      | splitWells.run(
        fromState: { id, state ->
          [
            input: state.output,
          ]
        },
        toState: { id, result, state ->
          state + result
        }
      )

      | setState(
        [
          "input": "barcode_path",
          "barcode": "barcode",
          "barcodesFasta": "barcodesFasta",
          "genomeDir": "genomeDir",
        ]
      )

      // TODO: Expand this into matching a whitelist/blacklist of barcodes
      // ... and turn into separate component
      | filter{ id, state -> state.barcode != "unknown" }

      | groupPairs.run(
        fromState: { id, state ->
          [
            input: state.input
          ]
        },
        toState: { id, result, state ->
          state + result
        }
      )

      // Does the sequencing platform use lanes?
      // Should those lanes be discriminated over?
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

      | groupWells.run(
        fromState: { id, state ->
          [
            input_r1: state.input_r1,
            input_r2: state.input_r2,
            well: state.barcode
          ]
        },
        toState: { id, result, state ->
          state + result
        }
      )

      | setState(
        [
          "input_r1": "output_r1",
          "input_r2": "output_r2",
          "wells": "wells",
          "barcodesFasta": "barcodesFasta",
          "genomeDir": "genomeDir",
        ]
      )

      | niceView()

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
