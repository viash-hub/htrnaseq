include { well_demultiplex } from params.rootDir + "/target/nextflow/workflows/well_demultiplex/main.nf"

params.resources_test =  "gs://viash-hub-test-data/htrnaseq/v1/"

workflow test_wf {
  resources_test_file = file(params.resources_test)
  output_ch = Channel.fromList([
      [
        id: "SRR14730301",
        input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
        barcodesFasta: resources_test_file.resolve("2-wells.fasta"),
      ],
      [
        id: "SRR14730302",
        input_r1: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R2_001.fastq"),
        barcodesFasta:  resources_test_file.resolve("2-wells.fasta"),
      ],
    ])
    | map { state -> [ state.id, state ] }
    | well_demultiplex.run(
      fromState: { id, state ->
        [
          input_r1: state.input_r1,
          input_r2: state.input_r2,
          barcodesFasta: state.barcodesFasta,
        ]
      },
      toState: { id, output, state ->
        output }
    )
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      "Output: $output"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 4 : "2 samples, each with 2 barcodes"
    }
}


workflow test_wf2 {
  resources_test_file = file(params.resources_test)
  output_ch = Channel.fromList([
      [
        id: "SRR14730301",
        input_r1:
          [
            resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
            resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R1_001.fastq"),
          ],
        input_r2:
          [
            resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
            resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R2_001.fastq"),
          ],
        barcodesFasta: resources_test_file.resolve("2-wells.fasta"),
      ],
    ])
    | map { state -> [ state.id, state ] }
    | well_demultiplex.run(
      fromState: { id, state ->
        [
          input_r1: state.input_r1,
          input_r2: state.input_r2,
          barcodesFasta: state.barcodesFasta,
        ]
      },
      toState: { id, output, state ->
        output }
    )
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      "Output: $output"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 2 : "1 samples, and two barcodes"
    }
}
