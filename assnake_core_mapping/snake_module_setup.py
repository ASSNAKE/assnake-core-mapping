import os
from assnake.api.snake_module import SnakeModule
from assnake.utils import read_yaml


from assnake_core_mapping.bwa.invocation_commands import map_bwa
from assnake_core_mapping.bowtie2.invocation_commands import map_bowtie2

this_dir = os.path.dirname(os.path.abspath(__file__))
snake_module = SnakeModule(name = 'assnake-core-mapping', 
                           install_dir = this_dir,
                           snakefiles = [
                               './bwa/workflow_fasta_from_db.smk', 
                               './general_mapping_operations.smk',
                               './bowtie2/workflow.smk'
                            ],
                           invocation_commands = [map_bwa, map_bowtie2],
                           initialization_commands = [],
                           wc_configs = [
                            #    read_yaml(os.path.join(this_dir, './wc_config.yaml'))
                               ])