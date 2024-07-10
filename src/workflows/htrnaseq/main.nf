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
          [
            output: result.output,
          ]
        }
      )

      | niceView()

      // Prep the handover to the mapping step
      | map { id, state ->
        // Group by the 'sample' id (ie barcode after demultiplexing)
        // Omit the 'unknown' group
        newState = state
          .output
          .collect{ p -> [ (p =~ /.*\\/([ACTG]*|unknown)_R?.*/)[0][1], p ] }
          .groupBy{ k, v -> k }
          .collect{ k, v -> [ k, v.collect{it[1]} ]}
          .findAll{ k, v -> k != "unknown"}
          .collect{ k, v -> [ barcode: k, fastq: v ] }
        [
          id, 
          [
            output: newState,
            out: newState.collect{it.fastq}.flatten().flatten() // to keep Viash IO happy
          ]
        ]
      }

      | niceView()

      | setState( [ "output": "out" ] )

  emit:
    output_ch
}
