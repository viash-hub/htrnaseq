workflow run_wf {

  take:
    input_ch

  main:
    // Count the number of input events per sample
    // Results from events with the same sample ID need to be concatenated.
    event_counts_ch = input_ch
      | map {id, state ->
        def new_state = state + ["event_id": id]
        def new_event = [state.sample_id, new_state]
        return new_event
      }
      | groupTuple(by: 0)
      | flatMap { id, states ->
          def orig_event_ids = states.collect{it.event_id}
          def new_events = orig_event_ids.collect{ orig_event_id ->
            [orig_event_id, ["n_events": states.size()]]
           }
          return new_events
      }


    // The number of events per sample needs is passed number to `groupTuple()`
    // so that it can emit the sample as soon as it is ready. This makes sure
    // that the samples are processed asynchronously. 
    output_ch = input_ch.join(event_counts_ch, failOnMismatch: true)
      | flatMap {id, state_demultiplex, state_event_counts ->
          assert state_demultiplex.input_r1.size() == state_demultiplex.input_r2.size(),
            "Expected output from well demultiplexing to contain equal amount or forward and reverse FASTQ files."
          def new_states = [state_demultiplex.input_r1, state_demultiplex.input_r2].transpose().collect{ fastq_files ->
            def (r1_file, r2_file) = fastq_files
            def regex = ~/^(\w+)_R[12]{1}_001\.fastq(\.gz)?$/
            def parsed_file_name = r1_file.name =~ regex
            def parsed_file_name_r2 = r2_file.name =~ regex
            def well_id = parsed_file_name[0][1]
            def well_id_r2 = parsed_file_name_r2[0][1]
  
            assert (well_id.length() != 0) && (well_id == well_id_r2)
            def new_state = state_demultiplex + [
              "input_r1": r1_file,
              "input_r2": r2_file,
              "event_id": id,
            ]
            def group_settings = groupKey("${state_demultiplex.sample_id}_${well_id}", state_event_counts.n_events)
            return [group_settings, new_state]

          }
        return new_states 
      }
      | groupTuple(by: 0, sort: "hash", remainder: true)
      | map {group_settings, sample_states -> 
        def input_r1 = sample_states.collect{it.input_r1}.flatten()
        def input_r2 = sample_states.collect{it.input_r2}.flatten()
        def event_ids = sample_states.collect{it.event_id}
        def sample_id_list = sample_states.collect{it.sample_id}.unique()
        assert sample_id_list.size() == 1
        def sample_id = sample_id_list[0]
        assert input_r1.size() == input_r2.size()  
        
        def new_state = [
          "input_r1": input_r1, 
          "input_r2": input_r2,
          "event_id": event_ids,
          "sample_id": sample_id,
        ]
        return [group_settings.target, new_state]
      } 
      | concat_text.run(
        directives: [label: ["verylowmem", "verylowcpu"]],
        key: "concat_samples_r1",
        runIf: {id, state -> state.input_r1.size() > 1},
        fromState: { id, state ->
          def output_file_name = state.input_r1[0].name
          [
            input: state.input_r1,
            gzip_output: false,
            output: output_file_name
          ]
        },
        toState: { id, result, state ->
          def newState = state + [ input_r1: [ result.output ] ]
          return newState
        }
      )
      | concat_text.run(
        directives: [label: ["verylowmem", "verylowcpu"]],
        key: "concat_samples_r2",
        runIf: {id, state -> state.input_r2.size() > 1},
        fromState: { id, state ->
          def output_file_name = state.input_r2[0].name
          [
            input: state.input_r2,
            gzip_output: false,
            output: output_file_name
          ]
        },
        toState: { id, result, state ->
          def newState = state + [ input_r2: [ result.output ] ]
          return newState
        }
      )
      | map {id, state ->
          def new_state = [state.sample_id, state]
          return new_state
      }
      | groupTuple(by: 0, sort: 'hash')
      | map {id, states ->
        def new_state = [
          "input_r1": states.collect{it.input_r1}.flatten(),
          "input_r2": states.collect{it.input_r2}.flatten(),
          "_meta": ["join_id": states[0].event_id[0]]
        ]
        return [id, new_state]
      }
      | setState(
        [
          "output_r1": "input_r1",
          "output_r2": "input_r2",
          "_meta": "_meta"
        ]
      )

  emit: 
    output_ch

}
