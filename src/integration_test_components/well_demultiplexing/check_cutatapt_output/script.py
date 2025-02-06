import dnaio
from operator import itemgetter
## VIASH START
par = {
}
## VIASH END

def assert_number_of_reads(reads):
    expected_number_of_reads = {
        "SRR14730301__A1": 165,
        "SRR14730301__B1": 194,
        "SRR14730302__A1": 141,
        "SRR14730302__B1": 213,
        "SRR14730302__unknown": 99646,
        "SRR14730301__unknown": 99641,
    }
    for input_id, expected_reads in expected_number_of_reads.items():
        num_reads = len(reads[input_id]) 
        assert num_reads == expected_reads, \
            f"Expected number of ouput reads for {input_id} to be {expected_reads}, was {num_reads}." 


def string_difference(string1, string2):
    result = 0
    for char1, char2 in zip(string1, string2, strict=True):
        if char1.lower() != char2.lower():
            result += 1
    return result


def assert_barcodes_not_removed(reads):
    barcodes = {
        "SRR14730301__A1": "ACACCGAATT",
        "SRR14730302__A1": "ACACCGAATT",
        "SRR14730301__B1": "GGCTATTGAT",
        "SRR14730302__B1": "GGCTATTGAT" 
    }
    for sample_id, barcode in barcodes.items():
        sample_reads = reads[sample_id]
        forward_reads = map(itemgetter(0), sample_reads)
        for i, forward_read in enumerate(forward_reads):
            read_sequence = forward_read.sequence
            read_barcode_start = read_sequence[: len(barcode)]
            # A 10% difference is allowed.
            assert string_difference(read_barcode_start, barcode) <= (0.1 * len(barcode)), \
                (f"Expected barcode {barcode} to be present for sample {sample_id} "
                 f"in read {i}. Found {read_barcode_start}")

def create_input_mapping(sample_ids, inputs_r1, inputs_r2):
    return {sample_id: [input_r1, input_r2] 
            for sample_id, input_r1, input_r2 
            in zip(sample_ids, inputs_r1, inputs_r2, strict=True)}

def read_input_files(input_mapping):
    expected_keys = {"SRR14730301__A1", "SRR14730301__B1",
                     "SRR14730302__A1", "SRR14730302__B1",
                     "SRR14730301__unknown", "SRR14730302__unknown"}
    difference = set(input_mapping.keys()) - expected_keys
    assert not difference, f"Found unexpected output id(s): {difference}"
    result = {}
    for input_id, input_files in input_mapping.items():
        input_r1, input_r2 = input_files
        # This reads the files into memory,
        # but they are reasonably small
        with dnaio.open(input_r1) as r1_reads, dnaio.open(input_r2) as r2_reads:
            for r1_read, r2_read in zip(r1_reads, r2_reads, strict=True):
                result.setdefault(input_id, []).append((r1_read, r2_read))
    return result


def main(par):
    inputs = create_input_mapping(par["ids"], par["fastq_r1"], par["fastq_r2"])
    reads = read_input_files(inputs)
    assert_number_of_reads(reads)
    assert_barcodes_not_removed(reads)

if __name__ == "__main__":
    main(par)