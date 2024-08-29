workflow run_wf {
    take:
    input_ch

    main:
    pool_ch = input_ch
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
           "output": state.output,
         ]
        },
        toState: { id, result, state ->
          state + [
            output: result.output,
          ]
        },
        directives: [label: ["midmem", "midcpu"]]
      )
      | setState(["output"])

    input_join_ch = input_ch
      | map {id, state ->
        [state.pool, id, state]
      }
    output_ch = input_join_ch.combine(pool_ch, by: 0)
      | map {pool, well_id, state_well, state_pool ->
        well_output = state_pool.output.findAll{star_output_dir ->
          def barcodes_list = []
          def barcode_file_regex = ~/.*\/raw\/barcodes\.tsv$/
          star_output_dir.eachFileRecurse{barcode_file ->
            if (barcode_file =~ barcode_file_regex) {
              assert barcode_file.countLines() == 1, "Expected only one barcode in a single STAR output."
              barcodes_list.add(barcode_file.text.trim())
            }
          }
          assert barcodes_list.size() == 1, "Exactly one file should have matched the barcodes file regex (found: $barcodes_list)."
          def barcode
          barcodes_list.each{ it -> barcode = it }
          return barcode == state_well.barcode
        }
        assert well_output.size() == 1, "Two or more outputs from the mapping seemed to have processed barcode '$barcode'."
        [well_id, ["output": well_output[0]]]
      }


    emit:
    output_ch
}