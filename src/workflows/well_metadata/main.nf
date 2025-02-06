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
        def well_id_matcher = /^([A-Za-z]+)0*([1-9]?[0-9]+)$/
        def entries_corrected_id = fasta_entries.collectEntries { it ->
          def unformatted_id = it.header
          def id_matched_to_format = unformatted_id =~ well_id_matcher
          assert (id_matched_to_format && id_matched_to_format.getCount() == 1), \
            "The FASTA headers must match the coordinate system of a well plate (e.g. A01, B01, ... or AA1, AB1, ...). Found: ${unformatted_id}"
          def id_letters = id_matched_to_format[0][1].toUpperCase()
          def id_numbers = id_matched_to_format[0][2]
          [it.seqString.replaceAll("[^ACGTacgt]", ""), "${id_letters}${id_numbers}".toString()]
        }
        def newState = state + [
          "n_wells": n_wells,
          "barcode_well_id_mapping": entries_corrected_id,
        ]
        [id, newState]
      }
      | flatMap{ id, state ->
        def new_events = state.star_mapping.collect{ star_output_dir ->
          def pool = id
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
          def well_id = state.barcode_well_id_mapping[barcode]
          assert well_id, "Could not find Well ID in FASTA file for barcode ${barcode}."
          def return_state = [
              "${pool}__${well_id}".toString(),
              [
                "barcode": barcode,
                "well_id": well_id,
                "pool": pool,
                "n_wells": state.n_wells,
                "output_r1": state.input_r1,
                "output_r2": state.input_r2,
                "well_star_mapping": star_output_dir,
                "_meta": ["join_id": pool]
              ]
          ]
        }
        return new_events
      }
      // Parse the file names to obtain metadata about the output
      | map{ id, state ->
        // Populate the new state
        def fastq_files = [state.output_r1, state.output_r2].transpose().findResult{ fastq_pair ->
          def (forward_fastq, reverse_fastq) = fastq_pair
          def fastq_r1_name = forward_fastq.name
          def fastq_r2_name = reverse_fastq.name
          // Get the well ID, and also check if it matches between the forward and reverse FASTQ
          def well_id = null
          [fastq_r1_name, fastq_r2_name].each { file_name ->
            def well_id_matcher = file_name =~ /^([A-Za-z0-9]*|unknown)_R?.*/
            assert well_id_matcher, \
              "Could not find Well ID in the name of FASTQ file ($file_name) output from cutadapt."
            def current_well_id = well_id_matcher[0][1]
            if (!well_id) {
              well_id = current_well_id
            } else {
              assert well_id == current_well_id,
                "Well ID for forward and reverse fastq file did not match! File names: ${fastq_r1_name} and ${fastq_r2_name}"
            }
          }
          assert (well_id != null), \
            "No Well ID could be deduced from files ${fastq_r1_name} and ${fastq_r2_name}."

          if (well_id == "unknown" || well_id != state.well_id) {
            return null
          }
          return fastq_pair
        }
        
        def new_state = state + [
          "output_r1": fastq_files[0],
          "output_r2": fastq_files[1]
        ]
        return [id, new_state]
      }
      | setState(["output_r1", "output_r2", "pool", "well_id", "n_wells", "barcode", "well_star_mapping", "_meta"])

  emit:
    output_ch
}