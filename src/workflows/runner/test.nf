import java.nio.file.Files
import nextflow.exception.WorkflowScriptErrorException

def viash_config = java.nio.file.Paths.get("${params.rootDir}/target/nextflow/workflows/runner/_viash.yaml")

def get_version(inputFile) {
  def yamlSlurper = new groovy.yaml.YamlSlurper()
  def loaded_viash_config = yamlSlurper.parse(file(inputFile))
  def version = (loaded_viash_config.version) ? loaded_viash_config.version : "unknown_version"
  println("HT-RNAseq version to be used: ${version}")
  return version
}

// Create temporary directory for the publish_dir if it is not defined
if (!params.containsKey("publish_dir") && params.containsKey("publishDir")) {
    params.publish_dir = params.publishDir
}

if (!params.containsKey("publish_dir")) {
    def tempDir = Files.createTempDirectory("demultiplex_runner_integration_test")
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
    params.publish_dir = tempDir
}

params.fastq_publish_dir = (file(params.publish_dir) / "fastq").toUriString()
params.results_publish_dir = (file(params.publish_dir) / "results").toUriString()
assert file(params.fastq_publish_dir).isEmpty()
assert file(params.results_publish_dir).isEmpty()

// The module inherits the parameters defined before the include statement, 
// therefore any parameters set afterwards will not be used by the module.

include { runner } from params.rootDir + "/target/nextflow/workflows/runner/main.nf"
params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  pipeline_version = get_version(viash_config)
  resources_test = file(params.resources_test)

  // results_publish_dir and results_publish_dir are inherited using params
  // but they must be defined in the state as well because viash will check
  // if all arguments are present in the hashmap
  output_ch = Channel.fromList([
    [
        id: "run_1",
        input: resources_test.resolve("10k/SRR14730301"),
        genomeDir: resources_test.resolve("genomeDir/subset/Homo_sapiens/v0.0.3"),
        barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
        annotation: resources_test.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        project_id: "foo",
        experiment_id: "bar",
        fastq_publish_dir: params.fastq_publish_dir,
        results_publish_dir: params.results_publish_dir,
    ],
    [
        id: "run_2",
        input:  resources_test.resolve("10k/SRR14730301"),
        genomeDir: resources_test.resolve("genomeDir/subset/Homo_sapiens/v0.0.3"),
        barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
        annotation: resources_test.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        project_id: "foo",
        experiment_id: "bar",
        fastq_publish_dir: params.fastq_publish_dir,
        results_publish_dir: params.results_publish_dir,
    ],
    [
        id: "run_3",
        input:resources_test.resolve("10k/SRR14730302"),
        genomeDir: resources_test.resolve("genomeDir/subset/Homo_sapiens/v0.0.3"),
        barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
        annotation: resources_test.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        project_id: "foo",
        experiment_id: "bar",
        fastq_publish_dir: params.fastq_publish_dir,
        results_publish_dir: params.results_publish_dir,
    ]
  ])
  | map { state -> [state.id, state]}
  | runner.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }

  tosortedlistch = output_ch
    | toSortedList()
    | map {events ->
        assert events.size() == 1, "Expected one events to be output, found ${events.size()}"
        events
    }
    | map {states -> 
        def output_state = states[0][1]
        assert output_state.eset_dir.listFiles().collect{it.name}.toSet() == ["VH02001612.rds", "VH02001614.rds"].toSet()
        assert output_state.star_output_dir.listFiles().collect{it.name}.toSet() == ["VH02001612", "VH02001614"].toSet()
        ["VH02001612", "VH02001614"].each{it ->
           assert output_state.star_output_dir.resolve(it).listFiles().collect{it.name}.toSet() == ["ACACCGAATT", "GGCTATTGAT"].toSet()
        }
        assert output_state.star_qc_metrics_dir.listFiles().collect{it.name}.toSet() == ["VH02001612.txt", "VH02001614.txt"].toSet()
        assert output_state.nrReadsNrGenesPerChrom_dir.listFiles().collect{it.name}.toSet() == ["VH02001612.txt", "VH02001614.txt"].toSet()
    }


    workflow.onComplete = {
        try {
            // Nexflow only allows exceptions generated using the 'error' function (which throws WorkflowScriptErrorException).
            // So in order for the assert statement to work (or allow other errors to let the tests to fail)
            // We need to wrap these in WorkflowScriptErrorException. See https://github.com/nextflow-io/nextflow/pull/4458/files
            // The error message will show up in .nextflow.log
            def fastq_subdir = file("${params.fastq_publish_dir}")
            assert fastq_subdir.isDirectory()
            def found_fastq_folders = fastq_subdir.listFiles().findAll{it.isDirectory()}.collect{it.name}.toSet()
            def expected_run_folders = ["run_1", "run_2", "run_3"].toSet()
            assert found_fastq_folders == expected_run_folders, "Expected correct run folders to be present. Found: ${found_fastq_folders}"
            unique_dirs = [
                "run1": files("${fastq_subdir.toUriString()}/run_1/*_htrnaseq_${pipeline_version}", type: 'any'),
                "run2": files("${fastq_subdir.toUriString()}/run_2/*_htrnaseq_${pipeline_version}", type: 'any'), 
                "run3": files("${fastq_subdir.toUriString()}/run_3/*_htrnaseq_${pipeline_version}", type: 'any'),  
            ]
            assert unique_dirs.every{it.value.size() == 1}
            unique_dirs = unique_dirs.collectEntries{k, v -> [k, v[0]]}

            assert unique_dirs.every{it.value.isDirectory()}
            assert unique_dirs.collect{_key, _value -> _value.name}.toSet().size() == 1
            def expected_samples = [
                "run1": "VH02001612",
                "run2": "VH02001612",
                "run3": "VH02001614"
            ]

            unique_dirs.each{_key, _value ->
                def expected_sample = expected_samples[_key]
                def expected_sample_dir = _value.resolve(expected_sample)
                assert expected_sample_dir.isDirectory(), "Expected ${expected_sample} to be present in ${_value}"
                def expected_fastq_files = [
                    "A1_R1_001.fastq", "A1_R2_001.fastq", 
                    "B1_R1_001.fastq", "B1_R2_001.fastq",
                    "unknown_R1_001.fastq", "unknown_R2_001.fastq"]
                def found_files = files("${expected_sample_dir}/*.fastq", type: 'any')
                assert found_files.every{it.isFile()}
                assert found_files.collect{it.name}.toSet() == expected_fastq_files.toSet()
            }

            def results_subdir = file("${params.results_publish_dir}")
            assert fastq_subdir.isDirectory()
            def expected_subdir = file("${results_subdir}/foo/bar/data_processed", type: 'any')
            assert expected_subdir.isDirectory()
            def expected_result_dir = files("${expected_subdir}/*_htrnaseq_${pipeline_version}", type: 'any')
            assert expected_result_dir.size() == 1
            expected_result_dir = expected_result_dir[0]
            assert expected_result_dir.isDirectory()
            def expected_esets = ["VH02001612.rds", "VH02001614.rds"]
            def found_esets = files("${expected_result_dir}/esets/*.rds", type: 'any')
            assert found_esets.size() == 2
            assert found_esets.collect{it.name}.toSet() == expected_esets.toSet()
            expected_table_filenames = ["VH02001612.txt", "VH02001614.txt"]
            def found_pdata = files("${expected_result_dir}/pData/*.txt", type: 'any')
            assert found_pdata.size() == 2
            assert found_pdata.collect{it.name}.toSet() == expected_table_filenames.toSet()
            def found_nr_genes_nr_reads = files("${expected_result_dir}/nrReadsNrGenesPerChrom/*.txt", type: 'any')
            assert found_nr_genes_nr_reads.size() == 2
            assert found_nr_genes_nr_reads.collect{it.name}.toSet() == expected_table_filenames.toSet() 
            def found_star_logs = files("${expected_result_dir}/starLogs/*.txt", type: 'any')
            assert found_star_logs.size() == 2
            assert found_star_logs.collect{it.name}.toSet() == expected_table_filenames.toSet()
            def star_output = file("${expected_result_dir}/star_output", type: 'any')
            assert star_output.isDirectory()
            
            assert files("${star_output}/*", type: 'any').collect{it.name}.toSet() == ["VH02001612", "VH02001614"].toSet()
            assert files("${star_output}/VH02001612/*", type: 'any').collect{it.name}.toSet() == ["ACACCGAATT", "GGCTATTGAT"].toSet()
            assert files("${star_output}/VH02001614/*", type: 'any').collect{it.name}.toSet() == ["ACACCGAATT", "GGCTATTGAT"].toSet()
            assert file("${expected_result_dir}/report.html").isFile()
            assert file("${expected_result_dir}/params.yaml").isFile()
            assert file("${expected_result_dir}/fData/fData.gencode.v41.annotation.gtf.gz.txt").isFile()

        } catch (Exception e) {
            throw new WorkflowScriptErrorException("Integration test failed!", e)
        }
    }
}


workflow test_wf_with_lanes {
  pipeline_version = get_version(viash_config)
  resources_test = file(params.resources_test)

  // results_publish_dir and results_publish_dir are inherited using params
  // but they must be defined in the state as well because viash will check
  // if all arguments are present in the hashmap
  output_ch = Channel.fromList([
    [
        id: "run_1",
        input: resources_test.resolve("10k_with_lanes/SRR14730301"),
        genomeDir: resources_test.resolve("genomeDir/subset/Homo_sapiens/v0.0.3"),
        barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
        annotation: resources_test.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        project_id: "foo",
        experiment_id: "bar",
        fastq_publish_dir: params.fastq_publish_dir,
        results_publish_dir: params.results_publish_dir,
    ],
    [
        id: "run_2",
        input:  resources_test.resolve("10k_with_lanes/SRR14730301"),
        genomeDir: resources_test.resolve("genomeDir/subset/Homo_sapiens/v0.0.3"),
        barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
        annotation: resources_test.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        project_id: "foo",
        experiment_id: "bar",
        fastq_publish_dir: params.fastq_publish_dir,
        results_publish_dir: params.results_publish_dir,
    ],
    [
        id: "run_3",
        input:resources_test.resolve("10k_with_lanes/SRR14730302"),
        genomeDir: resources_test.resolve("genomeDir/subset/Homo_sapiens/v0.0.3"),
        barcodesFasta: resources_test.resolve("2-wells-with-ids.fasta"),
        annotation: resources_test.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        project_id: "foo",
        experiment_id: "bar",
        fastq_publish_dir: params.fastq_publish_dir,
        results_publish_dir: params.results_publish_dir,
    ]
  ])
  | map { state -> [state.id, state]}
  | runner.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }

  tosortedlistch = output_ch
    | toSortedList()
    | map {events ->
        assert events.size() == 1, "Expected one events to be output, found ${events.size()}"
    }


    workflow.onComplete = {
        try {
            // Nexflow only allows exceptions generated using the 'error' function (which throws WorkflowScriptErrorException).
            // So in order for the assert statement to work (or allow other errors to let the tests to fail)
            // We need to wrap these in WorkflowScriptErrorException. See https://github.com/nextflow-io/nextflow/pull/4458/files
            // The error message will show up in .nextflow.log
            def fastq_subdir = file("${params.fastq_publish_dir}")
            assert fastq_subdir.isDirectory()
            def found_fastq_folders = fastq_subdir.listFiles().findAll{it.isDirectory()}.collect{it.name}.toSet()
            def expected_run_folders = ["run_1", "run_2", "run_3"].toSet()
            assert found_fastq_folders == expected_run_folders, "Expected correct run folders to be present. Found: ${found_fastq_folders}"
            unique_dirs = [
                "run1": files("${fastq_subdir.toUriString()}/run_1/*_htrnaseq_${pipeline_version}", type: 'any'),
                "run2": files("${fastq_subdir.toUriString()}/run_2/*_htrnaseq_${pipeline_version}", type: 'any'), 
                "run3": files("${fastq_subdir.toUriString()}/run_3/*_htrnaseq_${pipeline_version}", type: 'any'),  
            ]
            assert unique_dirs.every{it.value.size() == 1}
            unique_dirs = unique_dirs.collectEntries{k, v -> [k, v[0]]}

            assert unique_dirs.every{it.value.isDirectory()}
            assert unique_dirs.collect{_key, _value -> _value.name}.toSet().size() == 1
            def expected_samples = [
                "run1": "VH02001612",
                "run2": "VH02001612",
                "run3": "VH02001614"
            ]

            unique_dirs.each{_key, _value ->
                def expected_sample = expected_samples[_key]
                def expected_sample_dir = _value.resolve(expected_sample)
                assert expected_sample_dir.isDirectory(), "Expected ${expected_sample} to be present in ${_value}"
                def expected_fastq_files = [
                    "A1_R1_001.fastq", "A1_R2_001.fastq", 
                    "B1_R1_001.fastq", "B1_R2_001.fastq",
                    "unknown_R1_001.fastq", "unknown_R2_001.fastq"]
                def found_files = files("${expected_sample_dir}/*.fastq", type: 'any')
                assert found_files.every{it.isFile()}
                assert found_files.collect{it.name}.toSet() == expected_fastq_files.toSet()
            }

            def results_subdir = file("${params.results_publish_dir}")
            assert fastq_subdir.isDirectory()
            def expected_subdir = file("${results_subdir}/foo/bar/data_processed", type: 'any')
            assert expected_subdir.isDirectory()
            def expected_result_dir = files("${expected_subdir}/*_htrnaseq_${pipeline_version}", type: 'any')
            assert expected_result_dir.size() == 1
            expected_result_dir = expected_result_dir[0]
            assert expected_result_dir.isDirectory()
            def expected_esets = ["VH02001612.rds", "VH02001614.rds"]
            def found_esets = files("${expected_result_dir}/esets/*.rds", type: 'any')
            assert found_esets.size() == 2
            assert found_esets.collect{it.name}.toSet() == expected_esets.toSet()
            expected_table_filenames = ["VH02001612.txt", "VH02001614.txt"]
            def found_pdata = files("${expected_result_dir}/pData/*.txt", type: 'any')
            assert found_pdata.size() == 2
            assert found_pdata.collect{it.name}.toSet() == expected_table_filenames.toSet()
            def found_nr_genes_nr_reads = files("${expected_result_dir}/nrReadsNrGenesPerChrom/*.txt", type: 'any')
            assert found_nr_genes_nr_reads.size() == 2
            assert found_nr_genes_nr_reads.collect{it.name}.toSet() == expected_table_filenames.toSet() 
            def found_star_logs = files("${expected_result_dir}/starLogs/*.txt", type: 'any')
            assert found_star_logs.size() == 2
            assert found_star_logs.collect{it.name}.toSet() == expected_table_filenames.toSet()
            def star_output = file("${expected_result_dir}/star_output", type: 'any')
            assert star_output.isDirectory()
            
            assert files("${star_output}/*", type: 'any').collect{it.name}.toSet() == ["VH02001612", "VH02001614"].toSet()
            assert files("${star_output}/VH02001612/*", type: 'any').collect{it.name}.toSet() == ["ACACCGAATT", "GGCTATTGAT"].toSet()
            assert files("${star_output}/VH02001614/*", type: 'any').collect{it.name}.toSet() == ["ACACCGAATT", "GGCTATTGAT"].toSet()
            assert file("${expected_result_dir}/report.html").isFile()
            assert file("${expected_result_dir}/params.yaml").isFile()
            assert file("${expected_result_dir}/fData/fData.gencode.v41.annotation.gtf.gz.txt").isFile()

        } catch (Exception e) {
            throw new WorkflowScriptErrorException("Integration test failed!", e)
        }
    }
}