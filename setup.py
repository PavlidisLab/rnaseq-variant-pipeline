from setuptools import setup

setup(name='rnaseq-variant-pipeline',
      packages=['rnaseq_variant_pipeline'],
      install_requires=['luigi==2.8.13', 'bioluigi'],
      dependency_links=[
          'git+https://github.com/PavlidisLab/bioluigi.git#egg=bioluigi-0.0.0'])
