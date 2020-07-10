from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install




setup(
    name='assnake-core-mapping',
    version='0.0.2',
    packages=find_packages(),
    entry_points = {
        'assnake.plugins': ['assnake-core-mapping = assnake_core_mapping.snake_module_setup:snake_module']
    }
)