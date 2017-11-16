from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license_ = f.read()

setup(
    name='cellmlmanip',
    version='0.0.1',
    description='CellML loading and model equation manipulation',
    long_description=readme,
    author='',
    author_email='',
    url='https://github.com/ModellingWebLab/cellmlmanip',
    license=license_,
    find_packages=find_packages(exclude=('tests', 'docs'))
)
