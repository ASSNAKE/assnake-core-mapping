from snakemake.shell import shell
import json
import yaml

def bwa_params(params_loc):
    params_dict = {}
    with open(params_loc, 'r') as stream:
        try:
            params_dict = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return params_str



shell('''(bwa mem -M -t {threads} {params.ind_prefix} {input.r1} {input.r2} > {output.sam}) >{log} 2>&1''')
