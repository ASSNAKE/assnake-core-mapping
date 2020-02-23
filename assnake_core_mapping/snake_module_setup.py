import os
from assnake.api.snake_module import SnakeModule
from assnake.utils import read_yaml


from assnake_core_mapping.bwa.invocation_commands import map_bwa

this_dir = os.path.dirname(os.path.abspath(__file__))
snake_module = SnakeModule(name = 'assnake-core-mapping', 
                           install_dir = this_dir,
                           snakefiles = ['./bwa/workflow_fasta_from_db.smk', './general_mapping_operations.smk'],
                           invocation_commands = [map_bwa],
                           initialization_commands = [],
                           wc_configs = [
                            #    read_yaml(os.path.join(this_dir, './wc_config.yaml'))
                               ])