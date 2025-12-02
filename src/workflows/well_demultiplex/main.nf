workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      /*
      For each pool (i.e. event) in the channel, a list of R1 and R2 input
      reads is provided which correspond to the lanes. If there are multiple lanes,
      we can demultiplex into the wells for each lane in parallel. Therefore, cutadapt
      must be started multiple times and we need an event per lane. The events are
      created by taking the R1 and R2 pairs from the input lists. The index of the elements
      in these lists are added to the ID in order to make them unique.

      The same pools may be present in multiple sequencing runs. Here, the events must be unique
      across boths runs and samples. When called from the htrnaseq workflow; the events have the
      format '{pool_id}/{run_id}'
      */
      | flatMap {id, state ->
        assert state.input_r1.size() == state.input_r2.size(), \
          "Expected equal number of inputs for R1 and R2"
        if (state.input_r1.size() == 1) {
          // special case where we do not want to adjust the ID to add an index.
          // If we do add an index, the file paths will contain "_0", which
          // will not be removed. For the scenarios where we do have multiple lanes,
          // the files will be concatenated later and a new file path without the index
          // is created at that point.
          def newState = state + [
            "input_r1": state.input_r1[0],
            "input_r2": state.input_r2[0],
            "pool_and_run_id": id,
            "n_lanes": 1,
            "lane_sorting": 1,
          ]
          return [[id, newState]]
        }
        // Store the number of lanes that were encountered here in order to
        // group them together in an asynchronous manner later by providing
        // the expected number of events to be grouped to groupTuple.
        // see https://www.nextflow.io/docs/latest/reference/operator.html#grouptuple
        [state.input_r1, state.input_r2].transpose().withIndex().collect{ input_pair, index ->
          def single_input_r1 = input_pair[0]
          def single_input_r2 = input_pair[1]
          def newState = state + ["input_r1": single_input_r1,
                                  "input_r2": single_input_r2,
                                  "pool_and_run_id": id,
                                  "n_lanes": state.input_r1.size(),
                                  "lane_sorting": index]
          def newId = id + "_" + index
          [newId, newState]
        }
      }
      | cutadapt.run(
        directives: [label: ["midmem", "midcpu"]],
        fromState: { id, state ->
          [
            input: state.input_r1,
            input_r2: state.input_r2,
            no_indels: true,
            action: "none",
            front_fasta: state.barcodesFasta,
            output: "*_001.fastq.gz",
            error_rate: 0.10,
            demultiplex_mode: "single",
            output_r1: state.output_r1,
            output_r2: state.output_r2,
          ]
        },
        toState: { id, result, state ->
          def newState = [
            "pool_and_run_id": state.pool_and_run_id,
            "n_lanes": state.n_lanes,
            "output": result.output,
            "lane_sorting": state.lane_sorting,
          ]
          return newState
        }
      )
      | flatMap{ id, state ->
        if (!state.output) {
          log.error("Well demultiplexing seems to have yielded no FASTQ files! (ID $id, found ${state.output})")
        }
        // The output from cutadapt should be in the format {name}_R(1|2)_001.fastq.gz
        // See https://github.com/viash-hub/biobox/blob/952ff0843093b538cbfd6fefdecf2e7a0bc9e70b/src/cutadapt/script.sh#L226
        // Here, {name} is the name of the sequence in the barcode fasta: https://cutadapt.readthedocs.io/en/v5.0/guide.html#named-adapters
        state.output.collect{ p ->
          def path_as_string = p.name
          // Check for correct output file name format
          assert (path_as_string.endsWith("_R1_001.fastq.gz") || path_as_string.endsWith("_R2_001.fastq.gz")), \
            "Expected cutadapt output to contain files ending in '_R1_001.fastq.gz' or _R2_001.fastq.gz' only. Found: ${p}."
          // Detect read orientation from file name
          def pair_end = path_as_string.endsWith("_R1_001.fastq.gz") ? "R1" : "R2"
          // Use the start of the file
          def barcode_id = p.name - ~/_R(1|2)_001\.fastq\.gz$/
          def new_id = state.pool_and_run_id + "__" + barcode_id
          [
            new_id,
            [
              "pool_and_run_id": state.pool_and_run_id,
              "barcode_id": barcode_id,
              "output": p,
              "pair_end": pair_end,
              "n_lanes": state.n_lanes,
              "lane_sorting": state.lane_sorting,
            ]
          ]
        }
      }
      /*
      At this point, the events are provided on the smallest possible level,
      as each event represents the reads for a certain orientation from a
      particular lane and a single well. Here, we join these events back together
      on well level, gathering FASTQS across the lanes and read orientations.
      In order to make this joining as efficient as possible, the number of
      lanes which are expected to be gathered were stored in the state earlier.
      This way, the processing of a well can continue as as soon as all of
      the lanes have been gathered. The number of lanes times 2 (forward
      and reverse orientation) represents the total number of FASTQS (events)
      to be included for a certain well.
      */
      | map {id, state ->
          def group_key = groupKey(id, state.n_lanes * 2)
          return [group_key, state]
      }
      | groupTuple(by: 0, remainder: true, sort: {a, b ->
        // Make sure that the grouped states are in order,
        // meaning forward and reverse FASTQs are paired and the FASTQ
        // for the forward reads comes before the reverse reads FASTQ.
        if (a.lane_sorting == b.lane_sorting) {
          return a.pair_end <=> b.pair_end
        }
        return a.lane_sorting <=> b.lane_sorting
      })
      | map {group_key, states ->
        // The states are in one long flat list, group them into pairs
        // This assumes that the FASTQ files are already in order!
        // (See the 'sort' argument of groupTuple above)
        def output_pairs = states.collate(2)

        // Sanity check the state
        output_pairs.each{ pair ->
          assert pair.size() == 2, \
            "State error: expected FASTQ pairs as output from cutadapt, " +
            "found output state: $pair"
          def (first, second) = pair
          def should_be_the_same = [
            "barcode_id",
            "pool_and_run_id",
            "lane_sorting",
          ]
          should_be_the_same.each { attr_to_check ->
            first_attr = first.get(attr_to_check)
            second_attr = second.get(attr_to_check)
            assert first_attr == second_attr, \
              "State error: expected FASTQ pairs from cutadapt to have " +
              "the same detected ${attr_to_check}. Found: " +
              "$first_attr and $second_attr"
          }
        }
        // Forward and reverse reads should be designated
        // by 'R1' and 'R2', and sorted lexographically.
        def r1_output = output_pairs.collect{
          def forward_output = it[0].output
          assert forward_output.name.endsWith("R1_001.fastq.gz"), \
             "State error: expected first item from FASTQ pair to have " +
             "orientation 'R1', found ${forward_output.name}."
          return it[0].output
        }
        def r2_output = output_pairs.collect{
          def forward_output = it[1].output
          assert forward_output.name.endsWith("R2_001.fastq.gz"), \
             "State error: expected first item from FASTQ pair to have " +
             "orientation 'R2', found ${forward_output.name}."
          return it[1].output
        }
        assert r1_output.size() == r2_output.size()

        /* The lane sorting represents the order of the FASTQ files
           as provided by the input. The order of the FASTQ files should
           remain the same in the well output. This is because the result of STAR
           can differ based on the order of the reads in the FASTQ file.
           Even when the same reads are provided, the order of them matters.
        */
        def lane_sorting = output_pairs.it[0].lane_sorting
        def sorting_is_monotonically_increasing = lane_sorting.withIndex().every { i, idx ->
          idx == 0 || lane_sorting[idx - 1] <= i
        }
        assert sorting_is_monotonically_increasing, \
          "State error: expected the order of the FASTQ files after grouping " +
          "the cutadapt output to be the same as the order in the input. " +
          "Found sorting ${lane_sorting}, R1 output: ${r1_output}, R2 output: ${r2_output}."

        // Here we pick the state from the first item in the list of states
        // and overwrite the keys which are different across states
        def first_state = states[0]
        // The id is the sequence name for the barcode (from the FASTA file).
        def new_state = first_state + ["output_r1": r1_output, "output_r2": r2_output]
        // group_key.target is an attribute from an object created with nextflow's groupKey()
        // It is the Id by which the events were joined using groupTuple
        return [group_key.target, new_state]
      }
      | view {"State after running cutadapt: $it"}
      // TODO: Expand this into matching a whitelist/blacklist of barcodes
      // ... and turn into separate component

      // This is contatenation of the FASTQ files from different lanes
      // Concatenation of FASTQ files from the different runs is done later.
      | concat_text.run(
        directives: [label: ["verylowmem", "verylowcpu"]],
        key: "concat_txt_r1",
        runIf: {id, state -> state.output_r1.size() > 1},
        fromState: { id, state ->
          [
            input: state.output_r1,
            gzip_output: true,
            output: "${state.barcode_id}_R1_001.fastq.gz"
          ]
        },
        toState: { id, result, state ->
          def newState = state + [ output_r1: [ result.output ] ]
          return newState
        }
      )
      | concat_text.run(
        directives: [label: ["verylowmem", "verylowcpu"]],
        key: "concat_text_r2",
        runIf: {id, state -> state.output_r2.size() > 1},
        fromState: { id, state ->
          [
            input: state.output_r2,
            gzip_output: true,
            output: "${state.barcode_id}_R2_001.fastq.gz",
          ]
        },
        toState: { id, result, state ->
          def newState = state + [ output_r2: [ result.output ] ]
          return newState
        }
      )
      // Group the concatenated files back on pool level
      | map {id, state ->
        def new_event = [state.pool_and_run_id, state]
        return new_event
      }
      | groupTuple(by: 0, sort: {a, b -> a.barcode_id <=> b.barcode_id})
      | map {id, states ->
        def output_r1 = states.collect{it.output_r1}.flatten()
        def output_r2 = states.collect{it.output_r2}.flatten()
        def pools = states.collect{it.pool_and_run_id}
        assert pools.toSet().size() == 1, "Unexpected state: pool ID to be unique. Found: ${pools}."
        def output_state = ["output_r1": output_r1, "output_r2": output_r2, "pool_and_run_id": pools[0]]
        return [id, output_state]
      }
      // The concatenation of lanes happens in different work directories (each well is processed a different 
      // concat_text process). Here we make sure that the FASTQ files are gathered in a single directory. 
      // This could be skipped when no concatenation was done since cutadapt will output in a directory already.
      // But since we are copying symlinks most of the time there is almost no performance penalty here.
      | move_files_to_directory.run(
        fromState: { id, state ->
          [
            "input": state.output_r1 + state.output_r2,
            // Remark: the fastq path part may seem superfluous but is necessary for publising later
            "output": "fastq/${state.pool_and_run_id}/",
            "keep_symbolic_links": true
          ]
        },
        toState: {id, result, state ->
          def new_state = [
            "output_r1": state.output_r1.collect{result.output.resolve(it.name)},
            "output_r2": state.output_r2.collect{result.output.resolve(it.name)},
          ]
          new_state
        }
      )

  emit:
    output_ch
}
