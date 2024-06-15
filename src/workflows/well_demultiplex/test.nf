include { well_demultiplex } from params.rootDir + "/target/nextflow/workflows/well_demultiplex/main.nf"

workflow test_wf {
  output_ch = Channel.fromList([
        [
          // sample_sheet: resources_test.resolve("bcl_convert_samplesheet.csv"),
          // input: resources_test.resolve("iseq-DI/"),
          //sample_sheet: "https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/NovaSeq6000/SampleSheet.csv",
          input_r1: "test_resources/VH02001612_S9_R1_001-subsample.fastq",
          input_r2: "test_resources/VH02001612_S9_R2_001-subsample.fastq",
          barcodesFasta: "test_resources/barcodes.fasta",
          publish_dir: "output_dir/",
        ]
      ])
    | map { state -> [ "run", state ] }
    | well_demultiplex.run(
        toState: { id, output, state ->
          output + [ orig_input: state.input ] }
      )
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        "Output: $output"
      }
    // TODO: add assertions
}

