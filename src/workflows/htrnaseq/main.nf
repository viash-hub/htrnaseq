workflow run_wf {
  take:
    input_ch

  main:
    // Untar input genome if needed
    untar_ch = input_ch
      | toSortedList()
      | map {ids_and_states ->
        def ids = ids_and_states.collect{it[0]}
        def genomeDirs = ids_and_states.collect{it[1].genomeDir}.unique()
        assert genomeDirs.size() == 1, "Only one value for 'genomeDirs' should have been defined across all inputs."
        def genomeDir
        genomeDirs.each { genomeDir = it }
        ["unpack_genome", ["genomeDir": genomeDir, "original_ids": ids]]
      }
      | untar.run(
        runIf: {id, state ->
          def genomeDirStr = state.genomeDir.toString()
          genomeDirStr.endsWith(".tar.gz") || \
          genomeDirStr.endsWith(".tar") || \
          genomeDirStr.endsWith(".tgz") ? true : false
        },
        fromState: [
          "input": "genomeDir",
        ],
        toState: { id, result, state ->
          state + ["genomeDir": result.output]
        },
      )
      | flatMap{ id, state ->
        state.original_ids.collect{[it, state.genomeDir]}
      }

    output_ch = input_ch.join(untar_ch)
      | niceView()
      | map {id, state, genomeDir ->
        def newState = state + ["genomeDir": genomeDir]
        [id, newState]
      }
      | well_demultiplex.run(
        fromState: { id, state ->
          [
            input_r1: state.input_r1,
            input_r2: state.input_r2,
            barcodesFasta: state.barcodesFasta,
          ]
        },
        toState: { id, result, state ->
          state + [
              fastq_output: result.output,
              output: result.output,
            ]
        },
        directives: [label: ["midmem", "midcpu"]]
      )

      | splitWells.run(
        fromState: { id, state ->
          [
            input: state.output,
          ]
        },
        toState: { id, result, state ->
          state + result
        }
      )
      | setState(
        [
          "input": "barcode_path",
          "barcode": "barcode",
          "barcodesFasta": "barcodesFasta",
          "genomeDir": "genomeDir",
          "star_output": "star_output",
          "fastq_output": "fastq_output",
        ]
      )

      // TODO: Expand this into matching a whitelist/blacklist of barcodes
      // ... and turn into separate component
      | filter{ id, state -> state.barcode != "unknown" }

      | groupPairs.run(
        fromState: { id, state ->
          [
            input: state.input
          ]
        },
        toState: { id, result, state ->
          state + result
        }
      )

      // Does the sequencing platform use lanes?
      // Should those lanes be discriminated over?
      | groupLanes.run(
        fromState: { id, state ->
          [
            input_r1: state.r1,
            input_r2: state.r2
          ]
        },
        toState: { id, result, state ->
          state + result + [ multiple_lanes: result.output_r1.size() > 1]
        }
      )

      | setState(
        [
          "input_r1": "r1",
          "input_r2": "r2",
          "barcodesFasta": "barcodesFasta",
          "genomeDir": "genomeDir",
          "barcode": "barcode",
          "star_output": "star_output",
          "fastq_output": "fastq_output",
        ]
      )
      | concat_text.run(
        key: "concat_text_r1",
        fromState: { id, state ->
          [
            input: state.input_r1,
            gzip_output: true
          ]
        },
        toState: { id, result, state ->
          state + [ input_r1: result.output ]
        }
      )
      | concat_text.run(
        key: "concat_text_r2",
        fromState: { id, state ->
          [
            input: state.input_r2,
            gzip_output: true
          ]
        },
        toState: { id, result, state ->
          state + [ input_r2: result.output ]
        }
      )

      | groupWells.run(
        fromState: { id, state ->
          [
            input_r1: state.input_r1,
            input_r2: state.input_r2,
            well: state.barcode
          ]
        },
        toState: { id, result, state ->
          state + [ "wells": result.wells, "input_r1": result.output_r1, "input_r2": result.output_r2]
        }
      )
      | parallel_map.run(
        fromState: { id, state ->
         [
           input_r1: state.input_r1,
           input_r2: state.input_r2,
           genomeDir: state.genomeDir,
           barcodes: state.wells,
           wellBarcodesLength: 10,
           umiLength: 10,
           output: state.star_output[0],
         ]
        },
        toState: { id, result, state ->
          state + [
            star_output: result.output,
          ]
        },
        directives: [label: ["midmem", "midcpu"]]
      )
      | niceView()
      | setState(["star_output", "fastq_output"])
      
      //| niceView()
      //
      //| setState( [ "output": "out" ] )

  emit:
    output_ch
}
