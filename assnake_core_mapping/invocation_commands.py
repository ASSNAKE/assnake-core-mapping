import assnake.api.loaders
import assnake
from tabulate import tabulate
import click
from assnake.cli.cli_utils import sample_set_construction_options, add_options, generic_command_individual_samples, generate_result_list
import os, datetime 
import pandas as pd 

@click.command('bbmap-coverage-stats', short_help='Map your samples on genome using BWA MEM')
@add_options(sample_set_construction_options)
@click.option('--reference', 
                help='Reference to use', 
                required=True,
                type=click.STRING )
@click.option('--params', help='Parameters to use', default='def', type=click.STRING )
@click.option('--version', help='Version of BWA', default='0.7.17', type=click.STRING )


@click.pass_obj
def map_bwa(config, reference, params, version, **kwargs):
    sample_set, sample_set_name = generic_command_individual_samples(config, **kwargs)

    res_list = []

    for s in sample_set.samples_pd.to_dict(orient='records'):
        preprocessing = s['preproc']
        res_list.append( '{fs_prefix}/{df}/mapped/bwa__{version}__{params}/{reference}/{sample}/{preproc}/{sample}.bam'.format(
            fs_prefix = s['fs_prefix'].rstrip('\/'),
            df = s['df'],
            preproc = preprocessing,
            sample = s['fs_name'],
            reference = reference,
            params = params,
            version = version
        ))

    config['requests'] += res_list