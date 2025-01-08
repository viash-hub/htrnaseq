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
          def regex = ~/^(\w+)_S(\d+)_(L(\d+)_)?R(\d)_(\d+)\.fast[qa](\.gz)?$/
          def validFastq = f.name ==~ regex

          assert validFastq: "${f} does not match the regex ${regex}"

          def parsedFastq = f.name =~ regex

          return [
            fastq: f,
            sample_id: parsedFastq[0][1],
            sample: parsedFastq[0][2],
            lane: parsedFastq[0][3],
            read: parsedFastq[0][5],
          ]
        }

        println("Group paired fastq/fasta files")
        def grouped = processed_fastqs
          .groupBy({it.sample_id})
          .collect{ k, vs ->
            def r1 = vs.find{it.read == "1"}
            def r2 = vs.find{it.read == "2"}
            def newState = [
              r1: r1.fastq,
              r2: r2.fastq
            ]
            [ k, newState + vs[0].findAll{ it.key in ["sample", "lane"] } ]
          }

        grouped
          .collect{ [ it[0], it[1] + [ _meta: [ join_id: id ] ] ] }
      }

  emit: out_

}
