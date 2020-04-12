__author__ = "Dmitry Fedorov"
__copyright__ = "Copyright 2020, Dmitry Fedorov"
__email__ = "fedorov.de@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import os, yaml

def bw2_params(params_loc):
    params_str = ''
    params_dict = {}

    trail = params_dict.pop('TRAILING')
    if(len(trail.keys()) > 0):
        slw_str = 'TRAILING:{quality} '
        params_str += slw_str.format(**trail)
        
    
    return params_str

bowtie2_params = {}
with open(snakemake.input.params, 'r') as stream:
    try:
        bowtie2_params = yaml.load(stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print(exc)
bowtie2_input = ''

if bowtie2_params['as_single']:
    bowtie2_input = '-U {r1},{r2}'.format(r1 = snakemake.input.r1, r2 = snakemake.input.r2)
else:
    bowtie2_input = '-1 {r1} -2 {r2}'.format(r1 = snakemake.input.r1, r2 = snakemake.input.r2)


if bowtie2_params['use_single_reads']:
    if os.path.isfile(snakemake.params.r_s):
        if bowtie2_params['as_single']:
            bowtie2_input = bowtie2_input + ',' + snakemake.params.r_s
        else:
            bowtie2_input += ' -r {r_s}'.format(r_s = snakemake.input.r_s)
    else:
        print('NO SINGLE READS')

bowtie2_params_str = ''
if bowtie2_params['no-unal']:
    bowtie2_params_str += ' --no-unal'
if bowtie2_params['no-hd']:
    bowtie2_params_str += ' --no-hd'
if bowtie2_params['no-sq']:
    bowtie2_params_str += ' --no-sq'
if bowtie2_params['sam-no-qname-trunc']:
    bowtie2_params_str += ' --sam-no-qname-trunc'

mp2_cmd = \
'''bowtie2 -p {snakemake.threads} {bowtie2_params_str} \
-x {snakemake.params.ind_prefix} \
{bowtie2_input} \
-S {snakemake.output.sam} '''

shell("export PERL5LIB='';\n" + mp2_cmd + " > {snakemake.log} 2>&1")
