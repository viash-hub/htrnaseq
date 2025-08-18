import java.nio.file.Files
import nextflow.exception.WorkflowScriptErrorException

// Create temporary directory for the publish_dir if it is not defined
if (!params.fastq_publish_dir && params.publishDir) {
    params.fastq_publish_dir = params.publishDir
}

if (!params.fastq_publish_dir) {
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
params.resources_test = "gs://viash-hub-test-data/htrnaseq/v1/"

def date = new Date().format('yyyyMMdd_hhmmss')

workflow test_wf {
    output_ch = Channel.fromList([
        [
          id: "test",
          input: params.resources_test + "100k/SRR14730301",
          barcodesFasta: params.resources_test + "2-wells-with-ids.fasta",
          project_id: "test_project",
          experiment_id: "test_experiment",
          fastq_publish_dir: params.fastq_publish_dir,
        ]
    ])
    | map {event -> [event.id, event] }
    | well_demultiplex_runner.run(
        fromState: {id, state -> state }
    )
    workflow.onComplete = {
        try {
            // Nexflow only allows exceptions generated using the 'error' function (which throws WorkflowScriptErrorException).
            // So in order for the assert statement to work (or allow other errors to let the tests to fail)
            // We need to wrap these in WorkflowScriptErrorException. See https://github.com/nextflow-io/nextflow/pull/4458/files
            // The error message will show up in .nextflow.log

            def fastq_publish_dir = file("${params.fastq_publish_dir}/test_project/test_experiment/test/${date}_well_demultiplex_unknown_version/home/test/")
            assert fastq_publish_dir.isDirectory()
            def fastq_files = fastq_publish_dir.listFiles()
            assert fastq_files.collect{it.name}.toSet() == [
                "A1_R1_001.fastq",
                "A1_R2_001.fastq", 
                "B1_R1_001.fastq", 
                "B1_R2_001.fastq", 
                "unknown_R1_001.fastq", 
                "unknown_R2_001.fastq"
            ].toSet()
        } catch (Exception e) {
            throw new WorkflowScriptErrorException("Integration test failed!", e)
        }
    }

}