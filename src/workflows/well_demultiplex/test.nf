include { well_demultiplex } from params.rootDir + "/target/nextflow/workflows/well_demultiplex/main.nf"

base = "gs://viash-hub-test-data/htrnaseq/v1/"

workflow test_wf {
  output_ch = Channel.fromList([
      [
        id: "SRR14730301",
        input_r1: base + "100k/SRR14730301/VH02001612_S9_R1_001.fastq",
        input_r2: base + "100k/SRR14730301/VH02001612_S9_R2_001.fastq",
        barcodesFasta: base + "2-wells.fasta",
      ],
      [
        id: "SRR14730302",
        input_r1: base + "100k/SRR14730302/VH02001614_S8_R1_001.fastq",
        input_r2: base + "100k/SRR14730302/VH02001614_S8_R2_001.fastq",
        barcodesFasta: base + "2-wells.fasta",
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
      assert output[1].output.size == 6 : "6 fastq files should be generated: pairs for 2 wells and 1 unkonwn"
      "Output: $output"
    }
    | view { id, state ->
      assert state.output.size == 6 : "6 fastq files should be generated: pairs for 2 wells and 1 unkonwn"
      "State: $state"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 2 : "2 samples in should result in 2 results out"
    }
}

