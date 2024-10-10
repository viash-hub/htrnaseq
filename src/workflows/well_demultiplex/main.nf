workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | flatMap {id, state ->
         [state.input_r1, state.input_r2].transpose().withIndex().collect{ input_pair, index ->
          def single_input_r1 = input_pair[0]
          def single_input_r2 = input_pair[1]
          def newState = state + ["input_r1": single_input_r1, "input_r2": single_input_r2, "pool": id, "lane_sorting": index]
          def newId = id + "_" + index 
          [newId, newState]
        }
      }
      | cutadapt.run(
        // TODO: Remove hard-coded directives and replace with profiles
        directives: [
            cpus: 4
          ],
        fromState: { id, state ->
          def new_output = ("fastq_${state.lane_sorting}/*_001.fastq")
          [
            input: state.input_r1,
            input_r2: state.input_r2,
            no_indels: true,
            action: "none",
            front_fasta: state.barcodesFasta,
            output: new_output,
            error_rate: 0.10,
            demultiplex_mode: "single",
          ]
        },
        toState: { id, result, state ->
          [
            pool: state.pool,
            output: result.output,
            lane_sorting: state.lane_sorting,
          ]
        }
      )
      // Parse the file names to obtain metadata about the output
      | flatMap{ id, state ->
        def pool = state.pool
        state.output.collect{ p ->
          def barcode = (p =~ /.*\\/([ACTG]*|unknown)_R?.*/)[0][1]
          def pair_end = (p =~ /.*_(R[12])_.*/)[0][1]
          def lane = (p =~ /.*_(L\d+).*/) ? (p =~ /.*_(L\d+).*/)[0][1] : "NA"
          def new_id = pool + "__" + barcode
          [
            new_id,
            [
              pool: pool,
              barcode: barcode,
              output: p,
              lane: lane,
              pair_end: pair_end,
              lane_sorting: state.lane_sorting,
              _meta: [ join_id: pool ]
            ]
          ]
        }
      }
      // Group the outputs from across lanes
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def r1_output = states.findAll{ it.pair_end == "R1" }.collect{it.output}
        def r2_output = states.findAll{ it.pair_end == "R2" }.collect{it.output}
        def lane_sorting_r1 = states.findAll{ it.pair_end == "R1" }.collect{it.lane_sorting}
        def lane_sorting_r2 = states.findAll{ it.pair_end == "R2" }.collect{it.lane_sorting}

        // At this point, the lane_sorting hold the positios the items in r1_output and r2_output
        // should become in a new list.
        def r1_output_sorted = new ArrayList(r1_output.size())
        def r2_output_sorted = new ArrayList(r2_output.size())

        lane_sorting_r1.eachWithIndex { pos, index ->
          r1_output_sorted[pos] = r1_output[index]
        }

        lane_sorting_r2.eachWithIndex { pos, index ->
          r2_output_sorted[pos] = r2_output[index]
        }

        assert r1_output.size() == r2_output.size()
        // Here we pick the state from the first item in the list of states
        // and overwrite the keys which are different across states
        // TODO: we can assert that these keys are the same
        def first_state = states[0]
        def new_id = first_state.pool + "__" + first_state.barcode
        def new_state = first_state + ["output_r1": r1_output_sorted, "output_r2": r2_output_sorted]
        [new_id, new_state]
      }
      // TODO: Expand this into matching a whitelist/blacklist of barcodes
      // ... and turn into separate component
      | filter{ id, state -> state.barcode != "unknown" }
      | concat_text.run(
        key: "concat_txt_r1",
        runIf: {id, state -> state.output_r1.size() > 1},
        fromState: { id, state ->
          [
            input: state.output_r1,
            gzip_output: false,
            output: "${id}_R1.fastq"
          ]
        },
        toState: { id, result, state ->
          state + [ output_r1: [ result.output ] ]
        }
      )
      | concat_text.run(
        key: "concat_text_r2",
        runIf: {id, state -> state.output_r2.size() > 1},
        fromState: { id, state ->
          [
            input: state.output_r2,
            gzip_output: false,
            output: "${id}_R2.fastq",
          ]
        },
        toState: { id, result, state ->
          state + [ output_r2: [ result.output ] ]
        }
      )
      | setState(["pool", "barcode", "lane", "_meta", "output_r1", "output_r2"])

  emit:
    output_ch
}
