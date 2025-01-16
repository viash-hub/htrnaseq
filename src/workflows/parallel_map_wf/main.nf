workflow run_wf {
    take:
    input_ch

    main:
    pool_ch = input_ch
      | groupWells.run(
        fromState: [
          "input_r1": "input_r1",
          "input_r2": "input_r2",
          "well": "barcode",
          "pool": "pool",
        ],
        toState: [
          "wells": "wells",
          "input_r1": "output_r1",
          "input_r2": "output_r2",
        ]
      )
      | parallel_map.run(
        fromState: { id, state ->
         [
           "input_r1": state.input_r1,
           "input_r2": state.input_r2,
           "genomeDir": state.genomeDir,
           "barcodes": state.wells,
           "pool": state.pool,
           "umiLength": state.umi_length,
           "output": state.output,
         ]
        },
        toState: ["output": "output"],
        directives: ["label": ["highmem", "lowcpu"]],
      )
      | setState(["output", "pool"])

    // input_ch is on pool level, while parallel_map
    // outputs multiple events per pool. 
    // Join the results back to pool level
    input_join_ch = input_ch
      | map {id, state ->
        def newEvent = [state.pool, id, state]
        return newEvent
      }

    output_ch = input_join_ch.combine(pool_ch, by: 0)
      | map {pool, well_id, state_well, state_pool ->
        def well_output = state_pool.output.findAll{star_output_dir ->
          def barcodes_list = []
          // Get the barcode from the STAR file. 
          // One STAR output contains the results for one
          // well barcode. We can look for the barcode in
          // the 'Solo.out/Gene/raw/barcode.tsv' file.
          def barcodes_files = files("${star_output_dir}/Solo.out/Gene/raw/barcodes.tsv")
          assert barcodes_files.size() == 1, \
            "Exactly one file should have matched the barcodes files (found: $barcodes_files)."
          def barcode
          barcodes_files.each{ it ->
            assert it.countLines() == 1,
              "Expected only one barcode in a single STAR output."
            barcode = it.text.trim()
          }
          return barcode == state_well.barcode
        }
        assert well_output.size() == 1, \
          "Two or more outputs from the mapping seemed to have processed barcode '$barcode'."
        [well_id, ["output": well_output[0]]]
      }


    emit:
    output_ch
}