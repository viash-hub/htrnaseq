workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | niceView()
      | cutadapt.run(
          debug: true,
          directives: [
              cpus: 4
            ],
          fromState: { id, state ->
              [
                no_indels: true,
                action: "none",
                cores: 0,
                front_fasta: state.barcodesFasta,
                outputDir: "fastq",
                output: "*_001.fastq",
                input: state.input_r1,
                input_r2: state.input_r2,
              ]
            },
          toState: { id, state, result ->
              [
                output: result.output
              ]
            }
        ) | niceView()

  emit:
    output_ch
}

