import java.nio.file.Files
import nextflow.exception.WorkflowScriptErrorException

def viash_config = java.nio.file.Paths.get("${params.rootDir}/target/nextflow/workflows/well_demultiplex_runner/_viash.yaml")

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("HT-RNAseq version to be used: ${version}")
  return version
}

// Create temporary directory for the publish_dir if it is not defined
if (!params.containsKey("fastq_publish_dir") && params.containsKey("publishDir")) {
    params.fastq_publish_dir = params.publishDir
}

if (!params.containsKey("fastq_publish_dir")) {
    def tempDir = Files.createTempDirectory("well_demultiplex_runner_integration_test")
    println "Created temp directory: $tempDir"
    // Register shutdown hook to delete it on JVM exit
    Runtime.runtime.addShutdownHook(new Thread({
        try {
            // Delete directory recursively
            Files.walk(tempDir)
                .sorted(Comparator.reverseOrder())
                .forEach { Files.delete(it) }
            println "Deleted temp directory: $tempDir"
        } catch (Exception e) {
            println "Failed to delete temp directory: $e"
        }
    }))
    params.fastq_publish_dir = tempDir
}
// The module inherits the parameters defined before the include statement,
// therefore any parameters set afterwards will not be used by the module.

include { well_demultiplex_runner } from params.rootDir + "/target/nextflow/workflows/well_demultiplex_runner/main.nf"
params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
    pipeline_version = get_version(viash_config)
    resources_test = file(params.resources_test)

    // fastq_publish_dir is inherited using params but it must be defined in the
    // state as well because viash will check if all arguments are present in the hashmap
    output_ch = Channel.fromList([
        [
          id: "test",
          input: resources_test.resolve("10k/SRR14730301"),
          barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
          project_id: "test_project",
          experiment_id: "test_experiment",
          fastq_publish_dir: params.fastq_publish_dir,
        ]
    ])
    | map {event -> [event.id, event] }
    | well_demultiplex_runner.run(
        fromState: {id, state -> state }
    )
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, state]"
        "Output: $output"
    }

    output_ch
        | toSortedList()
        | map { events ->
            assert events.size() == 1, "Expected 1 event to be output, found ${events.size()}"
            assert events[0][1].run_params.isFile()
            events
        }

    workflow.onComplete = {
        try {
            // Nexflow only allows exceptions generated using the 'error' function (which throws WorkflowScriptErrorException).
            // So in order for the assert statement to work (or allow other errors to let the tests to fail)
            // We need to wrap these in WorkflowScriptErrorException. See https://github.com/nextflow-io/nextflow/pull/4458/files
            // The error message will show up in .nextflow.log

            def run_dir = file("${params.fastq_publish_dir}/test_project/test_experiment/test")
            assert run_dir.isDirectory()
            // The name of the output directory contains a timestamp, so it is matched with a glob.
            def output_dirs = files("${run_dir.toUriString()}/*_well_demultiplex_${pipeline_version}", type: 'any')
            assert output_dirs.size() == 1, "Expected a single output directory, found ${output_dirs}"

            def sample_dir = output_dirs[0].resolve("VH02001612")
            assert sample_dir.isDirectory(), "Expected 'VH02001612' to be present in ${output_dirs[0]}"
            def fastq_files = files("${sample_dir}/*", type: 'any')
            assert fastq_files.every{it.isFile()}
            assert fastq_files.collect{it.name}.toSet() == [
                "A1_R1_001.fastq.gz",
                "A1_R2_001.fastq.gz",
                "B1_R1_001.fastq.gz",
                "B1_R2_001.fastq.gz",
                "unknown_R1_001.fastq.gz",
                "unknown_R2_001.fastq.gz"
            ].toSet()
        } catch (Exception e) {
            throw new WorkflowScriptErrorException("Integration test failed!", e)
        }
    }

}
