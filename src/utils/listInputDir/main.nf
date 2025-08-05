workflow run_wf {

  take: in_

  main:

    out_ = in_
      | flatMap{ id, state ->
        println "Looking for fastq files in ${state.input}"
        def allFastqs = state.input
          .listFiles()
          .findAll{
            it.isFile() &&
            it.name ==~ /^.+\.fastq.gz$|^.+\.fastq$|^.+\.fasta$/
          }
        println "Found ${allFastqs.size()} fastq/fasta files in ${state.input}"
        assert allFastqs.size() > 0: "No fastq/fasta files found"

        println("Extracting information from fastq/fasta filenames")
        def processed_fastqs = allFastqs.collect { f ->
          def regex = ~/^(\S+)_S(\d+)_(L(\d+)_)?R(\d)_(\d+)\.fast[qa](\.gz)?$/
          def validFastq = f.name ==~ regex

          assert validFastq: "${f} does not match the regex ${regex}"

          def parsedFastq = f.name =~ regex
          def lane = parsedFastq[0][3]
          // Remove the trailing '_'
          def lane_remove_trailing = lane == null ? "" : lane.replaceAll('_$', "")
          def sample_id = parsedFastq[0][1]
          if (sample_id in ["Undetermined"] || (state.pools && !state.pools.isEmpty() && !state.pools.contains(sample_id))) {
            return null
          }
          return [
            "fastq": f,
            "sample_id": sample_id,
            "sample": parsedFastq[0][2],
            "lane": lane_remove_trailing,
            "read": parsedFastq[0][5],
          ]
        }

        println("Group paired fastq/fasta files")
        def grouped = processed_fastqs
          .findAll{it != null}
          .groupBy({it.sample_id}, {it.lane})
          .collectMany{ sample_id, states_per_lane ->
            def result = states_per_lane.collect{lane, lane_states ->
              assert lane_states.size() == 2, "Expected to find two fastq files per lane! " +
                "Found ${lane_states.size()}. State: ${states_per_lane}"
              def r1_state = lane_states.find({it.read == "1"})
              def r2_state = lane_states.find({it.read == "2"})
              def fastq_state = [
                "r1_output": r1_state.fastq,
                "r2_output": r2_state.fastq
              ]
              def new_state = fastq_state +
                r1_state.findAll{it.key in ["sample_id", "sample", "lane"]} + 
                ["_meta": ["join_id": id]]
              def new_id = lane?.trim() ? sample_id : "${sample_id}_${lane}".toString()
              return [new_id, new_state]
            }
            return result

          }
          return grouped

      }

  emit: out_

}
