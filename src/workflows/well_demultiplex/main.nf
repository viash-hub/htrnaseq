workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      /*
      Parse the fasta file containing the barcodes and do the following:
        - The sequence headers must not contain any whitespaces
        - The headers (Well IDs) must be unique
        - The barcodes must be unique
        - Store the number of barcodes in the state
        - Add a barcode to well ID (header) mapping to the state,
          in order to be able to retreive the well ID based on the FASTQ name after well demultiplexing
      */
      | map {id, state ->
        def n_wells = state.barcodesFasta.countFasta() as int
        // The header is the full header, the id is the part header up to the first whitespace character
        // We do not allow whitespace in the header of the fasta file, so assert this.
        def fasta_entries = state.barcodesFasta.splitFasta(
          record: ["id": true, "header": true, "seqString": true]
        )
        assert fasta_entries.every{it.id == it.header}, \
          "The barcodes FASTA headers must not contain any whitespace!"
        // Check if the fasta headers are unique
        def fasta_ids = fasta_entries.collect{it.id}
        assert fasta_ids.clone().unique() == fasta_ids, \
          "The barcodes FASTA entries must have a unique name!"
        // Check if the sequences are unique
        def fasta_sequences = fasta_entries.collect{it.seqString}
        assert fasta_sequences.clone().unique() == fasta_sequences, \
          "The barcodes FASTA sequences must be unique!"
        def well_id_matcher = /^([A-Za-z]+)0*([0-9]+)$/
        def entries_corrected_id = fasta_entries.collectEntries { it ->
          def unformatted_id = it.header
          def id_matched_to_format = unformatted_id =~ well_id_matcher
          assert (id_matched_to_format && id_matched_to_format.getCount() == 1), \
            "The FASTA headers must match the coordinate system of a well plate (e.g. A01, B01, ... or AA1, AB1, ...). Found: ${unformatted_id}"
          def id_letters = id_matched_to_format[0][1].toUpperCase()
          def id_numbers = id_matched_to_format[0][2]
          ["${id_letters}${id_numbers}", it.seqString]
        }
        def newState = state + [
          "n_wells": n_wells,
          "well_id_barcode_mapping": entries_corrected_id,
        ]
        [id, newState]
      }
      /*
      For each pool (i.e. event) in the channel, a list of R1 and R2 input
      reads is provided which correspond to the lanes. If there are multiple lanes,
      we can demultiplex into the wells for each lane in parallel. Therefore, cutadapt
      must be started multiple times and we need an event per lane. The events are
      created by taking the R1 and R2 pairs from the input lists. The index of the elements
      in these lists are added to the ID in order to make them unique.
      */
      | flatMap {id, state ->
        assert state.input_r1.size() == state.input_r2.size(), \
          "Expected equal number of inputs for R1 and R2"
        // Store the number of lanes that were encountered here in order to
        // group them together in an asynchronous manner later by providing
        // the expected number of events to be grouped to groupTuple.
        // see https://www.nextflow.io/docs/latest/reference/operator.html#grouptuple
        def n_lanes = state.input_r1.size()
        [state.input_r1, state.input_r2].transpose().withIndex().collect{ input_pair, index ->
          def single_input_r1 = input_pair[0]
          def single_input_r2 = input_pair[1]
          def newState = state + ["input_r1": single_input_r1,
                                  "input_r2": single_input_r2,
                                  "pool": id,
                                  "lane_sorting": index,
                                  "n_lanes": n_lanes]
          def newId = id + "_" + index
          [newId, newState]
        }
      }
      | cutadapt.run(
        directives: [label: ["highmem", "midcpu"]],
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
          def newState = [
            "pool": state.pool,
            "n_lanes": state.n_lanes,
            "output": result.output,
            "lane_sorting": state.lane_sorting,
            "n_wells": state.n_wells,
            "well_id_barcode_mapping": state.well_id_barcode_mapping,
          ]
          return newState
        }
      )
      // Parse the file names to obtain metadata about the output
      | flatMap{ id, state ->
        def pool = state.pool
        state.output.collect{ p ->
          def well_id_matcher = p =~ /.*\\/([A-Za-z0-9]*|unknown)_R?.*/
          assert well_id_matcher, \
            "Could not find Well ID in the name of FASTQ file ($p) output from cutadapt."
          def well_id = well_id_matcher[0][1]
          // Note: set the barcode to 'null' for reads that were put into 'unknown'
          def barcode = (well_id != "unknown") ? state.well_id_barcode_mapping[well_id].replaceAll("[^ACGTacgt]", "") : null
          assert (well_id == "unknown") || (barcode != null), \
            "After demultiplexing, no Well ID could be retreived for barcode ${barcode}."
          def pair_end_matcher = p =~ /.*_(R[12])_.*/
          assert pair_end_matcher, \
            "Could not find read orientation information in the name of the FASTQ file ($p) output from cutadapt."
          def pair_end = pair_end_matcher[0][1]
          def lane_matcher = p =~ /.*_(L\d+).*/
          def lane = lane_matcher ? lane_matcher[0][1] : "NA"
          def new_id = pool + "__" + well_id
          [
            new_id,
            [
              "pool": pool,
              "barcode": barcode,
              "well_id": well_id,
              "output": p,
              "lane": lane,
              "n_wells": state.n_wells,
              "pair_end": pair_end,
              "n_lanes": state.n_lanes,
              "lane_sorting": state.lane_sorting,
              "_meta": [ "join_id": pool ]
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
      | groupTuple(sort: {a, b ->
        // Make sure that the grouped states are in order,
        // meaning forward and reverse FASTQs are paired and the FASTQ
        // for the forward reads comes before the reverse reads FASTQ.
        if (a.lane_sorting == b.lane_sorting) {
          return a.pair_end <=> b.pair_end
        }
        return a.lane_sorting <=> b.lane_sorting
      })
      | map {_, states ->
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
            "barcode",
            "well_id",
            "lane",
            "pool",
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
          // Forward and reverse reads should be designated
          // by 'R1' and 'R2', and sorted lexographically.
          assert first.pair_end == "R1", \
            "State error: expected first item from FASTQ pair to have " +
             "orientation 'R1', found $first.pair_end"
          assert second.pair_end == "R2", \
            "State error: expected second item from FASTQ pair to have " +
             "orientation 'R2', found $second.pair_end"
        }

        def r1_output = output_pairs.collect{it[0].output}
        def r2_output = output_pairs.collect{it[1].output}
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
          "Found sorting $lane_sorting, R1 output: $r1_output, R2 output: $r2_output."

        // Here we pick the state from the first item in the list of states
        // and overwrite the keys which are different across states
        def first_state = states[0]
        def new_id = first_state.pool + "__" + first_state.well_id
        def new_state = first_state + ["output_r1": r1_output, "output_r2": r2_output]
        [new_id, new_state]
      }
      // TODO: Expand this into matching a whitelist/blacklist of barcodes
      // ... and turn into separate component
      | filter{ id, state -> state.well_id != "unknown" }
      | concat_text.run(
        directives: [label: ["lowmem", "lowcpu"]],
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
          def newState = state + [ output_r1: [ result.output ] ]
          return newState
        }
      )
      | concat_text.run(
        directives: [label: ["lowmem", "lowcpu"]],
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
          def newState = state + [ output_r2: [ result.output ] ]
          return newState
        }
      )
      | setState(["pool", "well_id", "n_wells", "barcode", "lane", "_meta", "output_r1", "output_r2"])

  emit:
    output_ch
}
