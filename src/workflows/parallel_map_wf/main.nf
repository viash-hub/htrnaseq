workflow run_wf {
    take:
    input_ch

    main:
    output_ch = input_ch
      | map {id, state -> [id, state + ["orig_id": id]]}
      | groupWells.run(
        fromState: { id, state ->
          [
            "input_r1": state.input_r1,
            "input_r2": state.input_r2,
            "well": state.barcode,
            "pool": state.pool,
          ]
        },
        toState: { id, result, state ->
          state + [ 
            "wells": result.wells,
            "input_r1": result.output_r1,
            "input_r2": result.output_r2,
            "_meta": ["join_id": state.orig_id]
          ]
        }
      )
      | parallel_map.run(
        fromState: { id, state ->
         [
           "input_r1": state.input_r1,
           "input_r2": state.input_r2,
           "genomeDir": state.genomeDir,
           "barcodes": state.wells,
           "pool": state.pool,
           "wellBarcodesLength": 10,
           "umiLength": 10,
           "output": state.output[0],
         ]
        },
        toState: { id, result, state ->
          state + [
            output: result.output,
          ]
        },
        directives: [label: ["midmem", "midcpu"]]
      )
      | setState(["output", "_meta"])
      
    emit:
    output_ch
}