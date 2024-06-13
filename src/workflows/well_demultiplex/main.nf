workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | niceView()
      | cutadapt.run(
          // TODO: Remove hard-coded directives and replace with profiles
          directives: [
              cpus: 4
            ],
          fromState: { id, state ->
              [
                input: state.input_r1,
                input_r2: state.input_r2,
                no_indels: true,
                action: "none",
                front_fasta: state.barcodesFasta,
                output: state.output
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
