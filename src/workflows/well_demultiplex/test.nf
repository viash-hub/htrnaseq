include { well_demultiplex } from params.rootDir + "/target/nextflow/workflows/well_demultiplex/main.nf"
include { check_cutadapt_output } from params.rootDir + "/target/nextflow/integration_test_components/well_demultiplexing/check_cutadapt_output/main.nf"

params.resources_test =  "gs://viash-hub-test-data/htrnaseq/v1/"

workflow test_wf {
  resources_test_file = file(params.resources_test)
  output_ch = Channel.fromList([
      [
        id: "SRR14730301",
        input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
        barcodesFasta: resources_test_file.resolve("2-wells-with-ids.fasta"),
      ],
      [
        id: "SRR14730302",
        input_r1: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R2_001.fastq"),
        barcodesFasta:  resources_test_file.resolve("2-wells-with-ids.fasta"),
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
      assert output[1].output_r1.size() == 3: "Expected 3 forward fastq files: 2 wells and 1 unknown"
      assert output[1].output_r2.size() == 3: "Expected 3 reverse fastq files: 2 wells and 1 unknown"
      "Output: $output"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 2 : "Should have found two pools!"
    }
    | map {output ->
      def ids = output.collect{it[0]}
      def states = output.collect{it[1]}
      def output_r1 = states.collect{it.output_r1}.flatten()
      def output_r2 = states.collect{it.output_r2}.flatten()
      def ids_pool_1 = states[0].output_r1.collect{ids[0] + "__" + (it.name - ~/_R(1|2)_001.fastq$/) } 
      def ids_pool_2 = states[1].output_r2.collect{ids[1] + "__" + (it.name - ~/_R(1|2)_001.fastq$/) } 

      def new_state = [
        "ids": ids_pool_1 + ids_pool_2,
        "fastq_r1": output_r1,
        "fastq_r2": output_r2
      ]
      ["integration_test_check", new_state]
    }
    | check_cutadapt_output.run(
      fromState: {id, state -> state}
    )
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
        barcodesFasta: resources_test_file.resolve("2-wells-with-ids.fasta"),
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
      assert output[1].output_r1.size() == 3: "Expected 3 forward fastq files: 2 wells and 1 unknown"
      assert output[1].output_r2.size() == 3: "Expected 3 reverse fastq files: 2 wells and 1 unknown"
      "Output: $output"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 1 : "Should have found 1 pool"
    }
}
