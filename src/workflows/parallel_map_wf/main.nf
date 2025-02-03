workflow run_wf {
    take:
    input_ch

    main:
    output_ch = input_ch
      | parallel_map.run(
        fromState: { id, state ->
         [
           "input_r1": state.input_r1,
           "input_r2": state.input_r2,
           "genomeDir": state.genomeDir,
           "barcodesFasta": state.barcodesFasta,
           "umiLength": state.umi_length,
           "output": state.output,
         ]
        },
        toState: ["output": "output"],
        directives: ["label": ["highmem", "lowcpu"]],
      )
      | setState(["output"])

    emit:
    output_ch
}