set -eo pipefail

## VIASH START
meta_executable=$(realpath "target/executable/parallel_map/parallel_map")
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

cat > "$TMPDIR/barcodes.fasta" <<'EOF'
>sample1
ACAGTCACAG
>sample2
CGGGTTTACC
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
mkdir -p "$TMPDIR/index"
STAR \
  ${meta_cpus:+--runThreadN $meta_cpus} \
  --runMode genomeGenerate \
  --outFileNamePrefix "$TMPDIR/index" \
  --genomeDir "$TMPDIR/index" \
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
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_1_dir/output_*"
popd

echo ">> Check if output directories exists"
sample1_out="$run_1_dir/output_ACAGTCACAG"
sample2_out="$run_1_dir/output_CGGGTTTACC"
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
cat "$sample1_out/Solo.out/Barcodes.stats" 
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
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_2_dir/output_gz_*" > /dev/null 2>&1
popd > /dev/null

echo ">> Check if output directories exists"
sample1_out="$run_2_dir/output_gz_ACAGTCACAG"
sample2_out="$run_2_dir/output_gz_CGGGTTTACC"
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


cat > "$TMPDIR/wrong_number_of_barcodes.fasta" <<'EOF'
>A1
ACAGTCACAG
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
    --barcodesFasta "$TMPDIR/wrong_number_of_barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_3_dir/output_gz_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
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
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_4_dir/" > /dev/null 2>&1 && echo "Expected non-zero exit code." && exit 1 
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
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_5_dir/output_run5_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
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
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_6_dir/output_run_6_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
set -eo pipefail
popd > /dev/null


echo "> Check that wrong FASTQ order is detected."
run_6_dir="$TMPDIR/run_7"
mkdir -p "$run_6_dir" 
pushd "$run_6_dir" > /dev/null
set +eo pipefail
"$meta_executable" \
    --input_r1 "$TMPDIR/sample2_R1.fastq.gz;$TMPDIR/sample1_R1.fastq.gz" \
    --input_r2 "$TMPDIR/sample1_R2.fastq;$TMPDIR/sample2_R2.fastq" \
    --genomeDir "$TMPDIR/index/" \
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_6_dir/output_run_6_*" > /dev/null 2>&1 && echo "Expected non-zero exit code " && exit 1
set -eo pipefail
popd > /dev/null


echo "> Check that order of input FASTQ files must not match with the order of barcodes"
run_8_dir="$TMPDIR/run_8"
mkdir -p "$run_8_dir"
pushd "$run_8_dir" > /dev/null
"$meta_executable" \
    --input_r1 "$TMPDIR/sample2_R1.fastq;$TMPDIR/sample1_R1.fastq" \
    --input_r2 "$TMPDIR/sample2_R2.fastq;$TMPDIR/sample1_R2.fastq" \
    --genomeDir "$TMPDIR/index/" \
    --barcodesFasta "$TMPDIR/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_8_dir/output_*" > /dev/null 2>&1 
popd

echo ">> Check if output directories exists"
sample1_out="$run_8_dir/output_ACAGTCACAG"
sample2_out="$run_8_dir/output_CGGGTTTACC"
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

echo "> Test reproducibility of output format"
# These files were also known to cause segfault issues in 2.7.5c
# See https://github.com/alexdobin/STAR/issues/1071
cat > "$TMPDIR/barcode_fail.fasta" <<'EOF'
>GAGACAAGGC
GAGACAAGGC
EOF
run_9_dir="$TMPDIR/run_9"
mkdir -p "$run_9_dir"

pushd "$run_9_dir" > /dev/null
"$meta_executable" \
    --input_r1 "$meta_resources_dir/test_output_v2_7_6/GAGACAAGGC_R1_001.fastq.gz" \
    --input_r2 "$meta_resources_dir/test_output_v2_7_6/GAGACAAGGC_R2_001.fastq.gz" \
    --genomeDir "$meta_resources_dir/gencode.v41.star.sparse" \
    --barcodesFasta "$meta_resources_dir/test_output_v2_7_6/barcodes.fasta" \
    --umiLength 10 \
    --runThreadN 2 \
    --output "$run_9_dir/output_*"
popd

echo ">> Check if output directories exists"
sample1_out="$run_9_dir/output_GAGACAAGGC"
assert_directory_exists "$sample1_out"

expected_output_compressed="$meta_resources_dir/test_output_v2_7_6/expected_output.tar.gz"
expected_output="$TMPDIR/expected_output"
tar -xvf "$expected_output_compressed" -C "$TMPDIR"

echo ">> Checking output barcodes"
barcodes_file="$sample1_out/Solo.out/Gene/filtered/barcodes.tsv"
if ! cmp -s "$barcodes_file" "$expected_output/Solo.out/Gene/filtered/barcodes.tsv"; then
  echo "Barcodes file is different."
  diff "$barcodes_file" "$expected_output/Solo.out/Gene/filtered/barcodes.tsv"
  exit 1
fi

echo ">> Checking count matrix"
matrix_file="$sample1_out/Solo.out/Gene/filtered/matrix.mtx"
if ! cmp -s "$matrix_file" "$expected_output/Solo.out/Gene/filtered/matrix.mtx"; then
  echo "Matrix file is different."
  diff "$matrix_file" "$expected_output/Solo.out/Gene/filtered/matrix.mtx" 
  exit 1
fi

echo ">> Checking features"
features_file="$sample1_out/Solo.out/Gene/filtered/features.tsv"
if ! cmp -s "$features_file" "$expected_output/Solo.out/Gene/filtered/features.tsv"; then
  echo "Features file is different."
  diff "$features_file" "$expected_output/Solo.out/Gene/filtered/features.tsv"
  exit 1
fi

features_stats="$sample1_out/Solo.out/Gene/Features.stats"
# In version 2.7.9a of STAR, better multi-gene read support was added.
# Replace :
## nMatch           1115
## nMatchUnique           1025
# With:
## nMatch           1025

sed -i -e 's/^\(\s*nMatch\s*\)1115/\11025/; /^\s*nMatchUnique.*/d' "$features_stats"
# In version 2.7.10b; some more fields were changed
## 1,11c1,10
## <                                         noUnmapped            649
## <                                        noNoFeature            301
## <                                       MultiFeature             90
## <                        subMultiFeatureMultiGenomic              0
## <                                 noTooManyWLmatches              0
## <                               noMMtoWLwithoutExact              0
## <                                         yesWLmatch           1025
## <                                 yessubWLmatchExact           1000
## <                        yessubWLmatch_UniqueFeature           1025
## <                                    yesCellBarcodes              1
## <                                            yesUMIs            986
## ---
## >                                          nUnmapped            649
## >                                         nNoFeature            301
## >                                      nAmbigFeature             90
## >                              nAmbigFeatureMultimap              0
## >                                           nTooMany              0
## >                                      nNoExactMatch              0
## >                                        nExactMatch           1000
## >                                             nMatch           1025
## >                                      nCellBarcodes              1
## >                                              nUMIs            986

sed -i -e 's/noUnmapped/ nUnmapped/' "$features_stats"
sed -i -e 's/noNoFeature/ nNoFeature/' "$features_stats"
sed -i -e 's/ MultiFeature/nAmbigFeature/' "$features_stats"
sed -i -e 's/subMultiFeatureMultiGenomic/      nAmbigFeatureMultimap/' "$features_stats"
sed -i -e 's/noTooManyWLmatches/          nTooMany/' "$features_stats"
sed -i -e 's/noMMtoWLwithoutExact/       nNoExactMatch/' "$features_stats"
sed -i -e 's/yesWLmatch/    nMatch/' "$features_stats"
sed -i -e 's/yesCellBarcodes/  nCellBarcodes/' "$features_stats"
sed -i -e 's/yesUMIs/  nUMIs/' "$features_stats"
sed -i -e 's/yessubWLmatchExact/       nExactMatch/' "$features_stats"
# yessubWLmatch_UniqueFeature corresponds to nMatchUnique from  2.7.9a
# sed -i -e 's/yessubWLmatch_UniqueFeature/               nMatchUnique/' "$features_stats"
sed -i -e '/yessubWLmatch_UniqueFeature/d' "$features_stats"
# Swap the order of the lines
printf '7m8\nw\n' | ed "$features_stats"

if ! cmp -s "$features_stats" "$expected_output/Solo.out/Gene/Features.stats"; then
  echo "Features stats are different."
  cat $features_stats
  cat "$expected_output/Solo.out/Gene/Features.stats" 
  diff -u "$features_stats" "$expected_output/Solo.out/Gene/Features.stats"
  exit 1
fi

echo ">> Checking summary"
# In version 2.7.8a; a bug was fixed that "slightly underestimated value of Q30 'Bases in RNA read'"
# See https://github.com/alexdobin/STAR/releases/tag/2.7.8a
summary="$sample1_out/Solo.out/Gene/Summary.csv"
sed -i -e 's/Q30 Bases in RNA read,0.934069/Q30 Bases in RNA read,0.924911/' "$summary"

# In version 2.7.9a; some fields in the summary were renamed
sed -i -e 's/Reads Mapped to Gene: Unique+Multipe Gene/Reads Mapped to Transcriptome: Unique+Multipe Genes/' "$summary"
sed -i -e 's/Reads Mapped to Gene: Unique Gene/Reads Mapped to Transcriptome: Unique Genes/' "$summary"

sed -i -e 's/Unique Reads in Cells Mapped to Gene/Reads in Cells Mapped to Unique Genes/' "$summary"
sed -i -e 's/Fraction of Unique Reads in Cells/Fraction of Reads in Cells/' "$summary"

sed -i -e 's/Mean Gene per Cell/Mean Genes per Cell/' "$summary"
sed -i -e 's/Median Gene per Cell/Median Genes per Cell/' "$summary"
sed -i -e 's/Total Gene Detected/Total Genes Detected/' "$summary"

# In version 2.7.10b; this change was introduced
sed -i -e 's/Reads Mapped to Gene: Unique+Multiple Gene,NoMulti/Reads Mapped to Transcriptome: Unique+Multipe Genes,0.539952/' "$summary"

if ! cmp -s "$summary" "$expected_output/Solo.out/Gene/Summary.csv"; then
  echo "Features stats are different."
  diff "$summary" "$expected_output/Solo.out/Gene/Summary.csv"
  exit 1
fi

echo ">> Checking ReadsPerGene.out.tab"
reads_per_gene="$sample1_out/ReadsPerGene.out.tab"
if ! cmp -s "$reads_per_gene" "$expected_output/ReadsPerGene.out.tab"; then
  echo "ReadsPerGene.out.tab is different."
  diff "$reads_per_gene" "$expected_output/ReadsPerGene.out.tab" 
  exit 1
fi

echo ">> Checking SJ.out.tab"
sj_out="$sample1_out/SJ.out.tab"
if ! cmp -s "$sj_out" "$expected_output/SJ.out.tab"; then
  echo  "SJ.out.tab is different."
  diff "$sj_out" "$expected_output/SJ.out.tab" 
  exit 1
fi

echo ">> Checking Log.final.out file"
log_file="$sample1_out/Log.final.out"
if ! cmp -s <(sed -e '1,/^$/ d' "$log_file") <(sed -e '1,/^$/ d' "$expected_output/Log.final.out"); then
  echo  "Log.final.out is different."
  diff <(sed -e '1,/^$/ d' "$log_file") <(sed -e '1,/^$/ d' "$expected_output/Log.final.out")
  exit 1
fi

echo ">> Checking alignment"
sam_file="$sample1_out/Aligned.sortedByCoord.out.bam"
edited_sam_file="$sample1_out/edit.Aligned.sortedByCoord.out.bam"
# In version 2.7.10b; a change was introduced that filled missing tags with '-'
awk 'BEGIN{OFS="\t"} /^@/ {print; next} {
    printf "%s", $1;
    for (i = 2; i <= 11; i++) printf "\t%s", $i;
    for (i = 12; i <= NF; i++) {
        if ($i != "GX:Z:-" && $i != "GN:Z:-") printf "\t%s", $i;
    }
    printf "\n";
}' <(samtools view -h "$sam_file") > "$edited_sam_file"

if ! cmp -s <(samtools view "$edited_sam_file") <(samtools view "$expected_output/Aligned.sortedByCoord.out.bam"); then
  echo "BAM file seems to be different"
  diff <(samtools view "$edited_sam_file" | head -n 10) <(samtools view "$expected_output/Aligned.sortedByCoord.out.bam" | head -n 10)
  exit 1
fi

echo ">> Checking barcodes.stats"

# In version 2.7.10b; some of these fields were renamed
## <                                        noNoAdapter              0
## <                                            noNoUMI              0
## <                                             noNoCB              0
## <                                            noNinCB              0
## <                                           noNinUMI              0
## <                                   noUMIhomopolymer              0
## <                                        noNoWLmatch              0
## <                                        noTooManyMM              0
## <                                 noTooManyWLmatches              0
## <                                    yesWLmatchExact           2028
## <                                yesOneWLmatchWithMM             37
## <                               yesMultWLmatchWithMM              0
## ---
## >                                         nNoAdapter              0
## >                                             nNoUMI              0
## >                                              nNoCB              0
## >                                             nNinCB              0
## >                                            nNinUMI              0
## >                                    nUMIhomopolymer              0
## >                                           nTooMany              0
## >                                           nNoMatch              0
## >                                nMismatchesInMultCB              0
## >                                        nExactMatch           2028
## >                                     nMismatchOneWL             37
## >                                  nMismatchToMultWL              0

stats_file="$sample1_out/Solo.out/Barcodes.stats"

sed -i -e 's/noNoAdapter/ nNoAdapter/' "$stats_file"
sed -i -e 's/noNoUMI/ nNoUMI/' "$stats_file"
sed -i -e 's/noNoCB/ nNoCB/' "$stats_file"
sed -i -e 's/noNinCB/ nNinCB/' "$stats_file"
sed -i -e 's/noNinUMI/ nNinUMI/' "$stats_file"
sed -i -e 's/noUMIhomopolymer/ nUMIhomopolymer/' "$stats_file"
sed -i -e 's/noTooManyMM/   nTooMany/' "$stats_file"
sed -i -e 's/noNoWLmatch/   nNoMatch/' "$stats_file"
sed -i -e 's/ noTooManyWLmatches/nMismatchesInMultCB/' "$stats_file"
sed -i -e 's/yesWLmatchExact/    nExactMatch/' "$stats_file"
sed -i -e 's/yesOneWLmatchWithMM/     nMismatchOneWL/' "$stats_file"
sed -i -e 's/yesMultWLmatchWithMM/   nMismatchToMultWL/' "$stats_file"
printf '7m8\nw\n' | ed "$stats_file"

if ! cmp -s "$stats_file" "$expected_output/Solo.out/Barcodes.stats"; then
  echo "Barcodes.stats file seems to be different"
  cat "$stats_file"
  cat "$expected_output/Solo.out/Barcodes.stats" 
  diff "$stats_file" "$expected_output/Solo.out/Barcodes.stats"
  exit 1
fi

echo ">> Checking raw/matrix.mtx"
matrix_file="$sample1_out/Solo.out/Gene/raw/matrix.mtx"
if ! cmp -s "$matrix_file" "$expected_output/Solo.out/Gene/raw/matrix.mtx"; then
  echo "matrix.mtx file seems to be different"
  diff "$matrix_file" "$expected_output/Solo.out/Gene/raw/matrix.mtx"
  exit 1
fi

