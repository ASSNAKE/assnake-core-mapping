import click, os
from assnake.core.result import Result

result = Result.from_location(name='bwa-mem',
                              description='Map your samples on reference using BWA MEM',
                              result_type='aligner',
                              location=os.path.dirname(os.path.abspath(__file__)), 
                              input_type='illumina_sample', 
                              additional_inputs=[
                                  click.option('--reference', help='Reference to use', required=True, type=click.STRING),
                                  click.option('--params', help='Parameters to use', default='def', type=click.STRING ),
                                  click.option('--version', help='Version of BWA', default='0.7.17', type=click.STRING )
                              ])
