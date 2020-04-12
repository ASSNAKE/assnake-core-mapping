import assnake.api.loaders
import assnake
from tabulate import tabulate
import click, glob
from assnake.core.sample_set import generic_command_individual_samples, generate_result_list
from assnake.cli.cli_utils import sample_set_construction_options, add_options
from assnake.core.result import Result

parameters = [p.split('/')[-1].replace('.json', '') for p in glob.glob('/data11/bio/databases/ASSNAKE/params/tmtic/*.json')]
additional_options = [
    click.option('--params',
            help='Parameters id to use. Available parameter sets: ' + str(parameters), 
            required=False, 
            default = 'def'),
    click.option('--reference', 
            help='Reference to use', 
            required=True,
            type=click.STRING ),
    click.option('--version', 
            help='Version of Bowtie2 to use', 
            required=True,
            default = 'v2.4.1',
            type=click.STRING)
]


@click.command('map-bowtie2', short_help='Read mapping with Bowtie2')
@add_options(sample_set_construction_options)
@add_options(additional_options)
@click.pass_obj
def map_bowtie2(config, **kwargs):
    wc_str = '{fs_prefix}/{df}/mapped/bowtie2__{params}__{version}/{reference}/{df_sample}/{preproc}/{df_sample}.sam'
    sample_set, sample_set_name = generic_command_individual_samples(config,  **kwargs)
    config['requests'] += generate_result_list(sample_set, wc_str, **kwargs)
    config['requested_results'] += [{'result': 'map-bowtie2', 'sample_set': sample_set}]
