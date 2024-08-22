set -eo pipefail

## VIASH START
meta_executable="target/executable/parallel_map/parallel_map"
## VIASH END

# Some helper functions
assert_directory_exists() {
  [ -d "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}

assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}

assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}

assert_file_contains_regex() {
  grep -q -E "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}

echo "> Prepare test data in $meta_temp_dir"
TMPDIR=$(mktemp -d --tmpdir="$meta_temp_dir")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Sample 1, barcode ACAGTCACAG, UMI CTACGGATGA
cat > "$TMPDIR/sample1_R1.fastq" <<'EOF'
@SAMPLE_1_SEQ_ID1
ACAGTCACAGCTACGGATGAGCCTCATAAGCCTCACACATCCGCGCCTATGTTGTGACTCTCTGTGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SAMPLE_1_SEQ_ID2
ACAGTCACAGCTACGGATGAGCCTCATAAGCCTCACACATCCGCGCCTATGTTGTGACTCTCTGTGAG
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > "$TMPDIR/sample1_R2.fastq" <<'EOF'
@SAMPLE_1_SEQ_ID1
CTCACAGAGAGTCACAACATAGGCGCGGATGTGTGAGGCTTATGAGGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SAMPLE_1_SEQ_ID2
CTCACAGAGAGTCACAACATAGGCGCGGATGTGTGAGGCTTATGAGGC
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

# Sample 2, barcode CGGGTTTACC, UMI GCTAGCTAGC
cat > "$TMPDIR/sample2_R1.fastq" << 'EOF'
@SAMPLE_2_SEQ_ID1
CGGGTTTACCGCTAGCTAGCCACCACTATGGTTGGCCGGTTAGTAGTGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SAMPLE_2_SEQ_ID2
CGGGTTTACCGCTAGCTAGCCACCACTATGGTTGGCCGGTTAGTAGTGT
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > "$TMPDIR/sample2_R2.fastq" <<'EOF'
@SAMPLE_2_SEQ_ID1
ACACTACTAACCGGCCAACCATAGTGGTG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SAMPLE_2_SEQ_ID2
ACACTACTAACCGGCCAACCATAGTGGTG
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

# Note that there is a sjdbGTFchrPrefix argument for STAR:
# prefix for chromosome names in a GTF file (default: '-')
cat > "$TMPDIR/genome.fasta" <<'EOF'
>1
TGGCATGAGCCAACGAACGCTGCCTCATAAGCCTCACACATCCGCGCCTATGTTGTGACTCTCTGTGAGCGTTCGTGGG
GCTCGTCACCACTATGGTTGGCCGGTTAGTAGTGTGACTCCTGGTTTTCTGGAGCTTCTTTAAACCGTAGTCCAGTCAA
TGCGAATGGCACTTCACGACGGACTGTCCTTAGCTCAGGGGA
EOF

cat > "$TMPDIR/genes.gtf" <<'EOF'
1    example_source  gene       0    72   .   +   .   gene_id "gene1"; gene_name: "GENE1;
1    example_source  exon       20   71   .   +   .   gene_id "gene1"; gene_name: "GENE1"; exon_id: gene1_exon1;
1    example_source  gene       80   160   .   +   .   gene_id "gene2"; gene_name: "GENE2;
1    example_source  exon       80   159   .   +   .   gene_id "gene2"; gene_name: "GENE2"; exon_id: gene2_exon1;

EOF

echo "> Generate index"
STAR \
  ${meta_cpus:+--runThreadN $meta_cpus} \
  --runMode genomeGenerate \
  --genomeDir "$TMPDIR/index/" \
  --genomeFastaFiles "$TMPDIR/genome.fasta" \
  --sjdbGTFfile "$TMPDIR/genes.gtf" \
  --genomeSAindexNbases 2 > /dev/null 2>&1


echo "> Run test 1"
run_1_dir="$TMPDIR/run_1"
mkdir -p "$run_1_dir"
pushd "$run_1_dir" > /dev/null
"$meta_executable" \
    --input_r1 "$TMPDIR/sample1_R1.fastq;$TMPDIR/sample2_R1.fastq" \
    --input_r2 "$TMPDIR/sample1_R2.fastq;$TMPDIR/sample2_R2.fastq" \
    --genomeDir "$TMPDIR/index/" \
    --barcodes "ACAGTCACAG;CGGGTTTACC" \
    --wellBarcodesLength 10 \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$TMPDIR/output_*" > /dev/null 2>&1 
popd

echo ">> Check if output directories exists"
sample1_out="$TMPDIR/output_ACAGTCACAG"
sample2_out="$TMPDIR/output_CGGGTTTACC"
assert_directory_exists "$sample1_out"
assert_directory_exists "$sample2_out"

echo ">> Check if output files have been created"
for sample in "$sample1_out" "$sample2_out"; do
  assert_file_exists "$sample/Aligned.sortedByCoord.out.bam" 
  assert_file_exists "$sample/Unmapped.out.mate1"
  assert_file_exists "$sample/Unmapped.out.mate2"
  assert_file_exists "$sample/Log.out"
  assert_file_exists "$sample/Log.final.out"
  assert_file_exists "$sample/ReadsPerGene.out.tab"
done 


echo ">> Check if Solo output is present"
for sample in "$sample1_out" "$sample2_out"; do
  assert_directory_exists "$sample1_out/Solo.out"
  assert_directory_exists "$sample1_out/Solo.out/Gene"
  assert_file_exists "$sample1_out/Solo.out/Barcodes.stats"
  assert_file_exists "$sample1_out/Solo.out/Gene/raw/barcodes.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/raw/features.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/raw/matrix.mtx"
  assert_file_exists "$sample1_out/Solo.out/Gene/filtered/barcodes.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/filtered/features.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/filtered/matrix.mtx"
done

echo ">> Check contents of output"
echo ">>> Sample 1"
assert_file_contains "$sample1_out/Solo.out/Barcodes.stats" "yesWLmatchExact              2"
assert_file_contains "$sample1_out/Log.final.out" "Uniquely mapped reads number |	2"
assert_file_contains "$sample1_out/Log.final.out" "Number of input reads |	2"

cat << EOF | cmp -s "$sample1_out/Solo.out/Gene/filtered/barcodes.tsv" || { echo "Barcodes file is different"; exit 1; }
ACAGTCACAG
EOF

cat << EOF | cmp -s "$sample1_out/Solo.out/Gene/filtered/features.tsv" || { echo "Features file is different"; exit 1; }
gene1	gene1	Gene Expression
gene2	gene2	Gene Expression
EOF

cat << EOF | cmp -s "$sample1_out/Solo.out/Gene/filtered/matrix.mtx" || { echo "Matrix file is different"; exit 1; }
%%MatrixMarket matrix coordinate integer general
%
2 1 1
1 1 1
EOF

echo ">>> Sample 2"
assert_file_contains "$sample2_out/Solo.out/Barcodes.stats" "yesWLmatchExact              2"
assert_file_contains "$sample2_out/Log.final.out" "Uniquely mapped reads number |	2"
assert_file_contains "$sample2_out/Log.final.out" "Number of input reads |	2"

cat << EOF | cmp -s "$sample2_out/Solo.out/Gene/filtered/barcodes.tsv" || { echo "Barcodes file is different"; exit 1; }
CGGGTTTACC
EOF

cat << EOF | cmp -s "$sample2_out/Solo.out/Gene/filtered/features.tsv" || { echo "Features file is different"; exit 1; }
gene1	gene1	Gene Expression
gene2	gene2	Gene Expression
EOF

cat << EOF | cmp -s "$sample2_out/Solo.out/Gene/filtered/matrix.mtx" || { echo "Matrix file is different"; exit 1; }
%%MatrixMarket matrix coordinate integer general
%
2 1 1
2 1 1
EOF

echo "> Run test 2 (compressed input)"
gzip -c "$TMPDIR/sample1_R1.fastq" > "$TMPDIR/sample1_R1.fastq.gz"
gzip -c "$TMPDIR/sample2_R1.fastq" > "$TMPDIR/sample2_R1.fastq.gz"
gzip -c "$TMPDIR/sample1_R2.fastq" > "$TMPDIR/sample1_R2.fastq.gz"
gzip -c "$TMPDIR/sample2_R2.fastq" > "$TMPDIR/sample2_R2.fastq.gz"

run_2_dir="$TMPDIR/run_2"
mkdir -p "$run_2_dir" 
pushd "$run_2_dir" > /dev/null
"$meta_executable" \
    --input_r1 "$TMPDIR/sample1_R1.fastq.gz;$TMPDIR/sample2_R1.fastq.gz" \
    --input_r2 "$TMPDIR/sample1_R2.fastq.gz;$TMPDIR/sample2_R2.fastq.gz" \
    --genomeDir "$TMPDIR/index/" \
    --barcodes "ACAGTCACAG;CGGGTTTACC" \
    --wellBarcodesLength 10 \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$TMPDIR/output_gz_*" > /dev/null 2>&1
popd > /dev/null

echo ">> Check if output directories exists"
sample1_out="$TMPDIR/output_gz_ACAGTCACAG"
sample2_out="$TMPDIR/output_gz_CGGGTTTACC"
assert_directory_exists "$sample1_out"
assert_directory_exists "$sample2_out"

echo ">> Check if output files have been created"
for sample in "$sample1_out" "$sample2_out"; do
  assert_file_exists "$sample/Aligned.sortedByCoord.out.bam" 
  assert_file_exists "$sample/Unmapped.out.mate1"
  assert_file_exists "$sample/Unmapped.out.mate2"
  assert_file_exists "$sample/Log.out"
  assert_file_exists "$sample/Log.final.out"
  assert_file_exists "$sample/ReadsPerGene.out.tab"
done 


echo ">> Check if Solo output is present"
for sample in "$sample1_out" "$sample2_out"; do
  assert_directory_exists "$sample1_out/Solo.out"
  assert_directory_exists "$sample1_out/Solo.out/Gene"
  assert_file_exists "$sample1_out/Solo.out/Barcodes.stats"
  assert_file_exists "$sample1_out/Solo.out/Gene/raw/barcodes.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/raw/features.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/raw/matrix.mtx"
  assert_file_exists "$sample1_out/Solo.out/Gene/filtered/barcodes.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/filtered/features.tsv"
  assert_file_exists "$sample1_out/Solo.out/Gene/filtered/matrix.mtx"
done

echo ">> Check contents of output"
echo ">>> Sample 1"
assert_file_contains "$sample1_out/Solo.out/Barcodes.stats" "yesWLmatchExact              2"
assert_file_contains "$sample1_out/Log.final.out" "Uniquely mapped reads number |	2"
assert_file_contains "$sample1_out/Log.final.out" "Number of input reads |	2"

cat << EOF | cmp -s "$sample1_out/Solo.out/Gene/filtered/barcodes.tsv" || { echo "Barcodes file is different"; exit 1; }
ACAGTCACAG
EOF

cat << EOF | cmp -s "$sample1_out/Solo.out/Gene/filtered/features.tsv" || { echo "Features file is different"; exit 1; }
gene1	gene1	Gene Expression
gene2	gene2	Gene Expression
EOF

cat << EOF | cmp -s "$sample1_out/Solo.out/Gene/filtered/matrix.mtx" || { echo "Matrix file is different"; exit 1; }
%%MatrixMarket matrix coordinate integer general
%
2 1 1
1 1 1
EOF

echo ">>> Sample 2"
assert_file_contains "$sample2_out/Solo.out/Barcodes.stats" "yesWLmatchExact              2"
assert_file_contains "$sample2_out/Log.final.out" "Uniquely mapped reads number |	2"
assert_file_contains "$sample2_out/Log.final.out" "Number of input reads |	2"

cat << EOF | cmp -s "$sample2_out/Solo.out/Gene/filtered/barcodes.tsv" || { echo "Barcodes file is different"; exit 1; }
CGGGTTTACC
EOF

cat << EOF | cmp -s "$sample2_out/Solo.out/Gene/filtered/features.tsv" || { echo "Features file is different"; exit 1; }
gene1	gene1	Gene Expression
gene2	gene2	Gene Expression
EOF

cat << EOF | cmp -s "$sample2_out/Solo.out/Gene/filtered/matrix.mtx" || { echo "Matrix file is different"; exit 1; }
%%MatrixMarket matrix coordinate integer general
%
2 1 1
2 1 1
EOF


echo "> Check that wrong number of barcodes are detected."
run_3_dir="$TMPDIR/run_3"
mkdir -p "$run_3_dir" 
pushd "$run_3_dir" > /dev/null
set +eo pipefail
"$meta_executable" \
    --input_r1 "$TMPDIR/sample1_R1.fastq.gz;$TMPDIR/sample2_R1.fastq.gz" \
    --input_r2 "$TMPDIR/sample1_R2.fastq.gz;$TMPDIR/sample2_R2.fastq.gz" \
    --genomeDir "$TMPDIR/index/" \
    --barcodes "ACAGTCACAG" \
    --wellBarcodesLength 10 \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$TMPDIR/output_gz_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
set -eo pipefail
popd > /dev/null

echo "> Check that missing wildcard character is detected."
run_4_dir="$TMPDIR/run_4"
mkdir -p "$run_4_dir" 
pushd "$run_4_dir" > /dev/null
set +eo pipefail
"$meta_executable" \
    --input_r1 "$TMPDIR/sample1_R1.fastq.gz;$TMPDIR/sample2_R1.fastq.gz" \
    --input_r2 "$TMPDIR/sample1_R2.fastq.gz;$TMPDIR/sample2_R2.fastq.gz" \
    --genomeDir "$TMPDIR/index/" \
    --barcodes "ACAGTCACAG;CGGGTTTACC" \
    --wellBarcodesLength 10 \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$TMPDIR/output_run4" > /dev/null 2>&1 && echo "Expected non-zero exit code." && exit 1 
set -eo pipefail
popd > /dev/null

echo "> Check that a mismatch in the length of the input mates is detected."
empty_input_file="$TMPDIR/empty.fastq"
touch "$empty_input_file"
run_5_dir="$TMPDIR/run_5"
mkdir -p "$run_5_dir" 
pushd "$run_5_dir" > /dev/null
set +eo pipefail
"$meta_executable" \
    --input_r1 "$TMPDIR/sample1_R1.fastq;$empty_input_file" \
    --input_r2 "$TMPDIR/sample1_R2.fastq;$TMPDIR/sample2_R2.fastq" \
    --genomeDir "$TMPDIR/index/" \
    --barcodes "ACAGTCACAG;CGGGTTTACC" \
    --wellBarcodesLength 10 \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$TMPDIR/output_run5_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
set -eo pipefail
popd > /dev/null

echo "> Check that wrong number of input files is detected."
run_6_dir="$TMPDIR/run_6"
mkdir -p "$run_6_dir" 
pushd "$run_6_dir" > /dev/null
set +eo pipefail
"$meta_executable" \
    --input_r1 "$TMPDIR/sample1_R1.fastq" \
    --input_r2 "$TMPDIR/sample1_R2.fastq;$TMPDIR/sample2_R2.fastq" \
    --genomeDir "$TMPDIR/index/" \
    --barcodes "ACAGTCACAG;CGGGTTTACC" \
    --wellBarcodesLength 10 \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$TMPDIR/output_run_6_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
set -eo pipefail
popd > /dev/null


