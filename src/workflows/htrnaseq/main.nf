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

  emit:
    output_ch
}
