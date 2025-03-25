import re
import yaml

## VIASH START
par = {
    "id": "sample_one",
    "params": ['workflow_id=sample_one',
              'umi_length=10',
              'fastq_output_r1[0]=fastq/*_R1_001.fastq',
              'fastq_output_r2[0]=fastq/*_R2_001.fastq',
              'star_output[0]=star.$id/*',
              'nrReadsNrGenesPerChrom=nrReadsNrGenesPerChrom.$id.txt',
              'star_qc_metrics=starLogs.$id.txt',
              'eset=eset.$id.rds',
              'f_data=fData.$id.tsv',
              'p_data=pData.$id.tsv',
              'html_report=report.$id.html',
              'input_r1[0]=/home/jakubmajercik/Data_Intuitive/htrnaseq/v1/100k/SRR14730301/VH02001612_S9_R1_001.fastq',
              'input_r1[1]=/home/jakubmajercik/Data_Intuitive/htrnaseq/v1/100k/SRR14730301/VH02001612_S9_R2_001.fastq',
              'barcodesFasta=/home/jakubmajercik/Data_Intuitive/htrnaseq/v1/360-wells-with-ids.fasta',
              'genomeDir=/home/jakubmajercik/Data_Intuitive/htrnaseq/v1/genomeDir/gencode.v41.star.sparse',
              'annotation=/home/jakubmajercik/Data_Intuitive/htrnaseq/v1/genomeDir/gencode.v41.annotation.gtf.gz'],
    "output": 'params_out.yaml'
}
## VIASH END

class Dumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(Dumper, self).increase_indent(flow, False)

params = {}

for param in par['params']:
    param = param.replace('$id', par['id'])
    key, value = param.split('=')
    
    array_match = re.match(r'(.+)\[(\d+)\]$', key)
    if array_match:
        base_key = array_match.group(1)
        index = int(array_match.group(2))
        
        if base_key not in params:
            params[base_key] = []
            
        while len(params[base_key]) <= index:
            params[base_key].append(None)
            
        params[base_key][index] = value
    else:
        params[key] = value
        
with open(par["output"], 'w') as f:
    yaml.dump(params, f, default_flow_style=False, Dumper=Dumper)