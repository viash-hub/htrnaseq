workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
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
            output: "fastq/*_001.fastq"
          ]
        },
        toState: { id, result, state ->
          [
            output: result.output,
          ]
        }
      )
      // Parse the file names to obtain metadata about the output
      | flatMap{ id, state ->
        state.output.collect{ p ->
          def barcode = (p =~ /.*\\/([ACTG]*|unknown)_R?.*/)[0][1]
          def pair_end = (p =~ /.*_(R[12])_.*/)[0][1]
          def lane = (p =~ /.*_(L\d+).*/) ? (p =~ /.*_(L\d+).*/)[0][1] : "NA"
          def new_id = id + "__" + barcode
          [
            new_id,
            [
              pool: id,
              barcode: barcode,
              output: p,
              lane: lane,
              pair_end: pair_end,
              _meta: [ join_id: id ]
            ]
          ]
        }
      }
      // Group the outputs from across lanes
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def r1_output = states.findAll{ it.pair_end == "R1" }.collect{it.output}
        def r2_output = states.findAll{ it.pair_end == "R2" }.collect{it.output}
        assert r1_output.size() == r2_output.size()
        // Here we pick the state from the first item in the list of states
        // and overwrite the keys which are different across states
        // TODO: we can assert that these keys are the same
        def first_state = states[0]
        def new_id = first_state.pool + "__" + first_state.barcode
        def new_state = first_state + ["output_r1": r1_output, "output_r2": r2_output]
        [new_id, new_state]
      }
      | setState(["pool", "barcode", "lane", "_meta", "output_r1", "output_r2"])

  emit:
    output_ch
}
