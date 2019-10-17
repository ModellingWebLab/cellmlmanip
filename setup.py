from setuptools import find_packages, setup


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
    packages=find_packages(exclude=('tests', 'docs')),
    package_data={'cellmlmanip': ['cellml_units.txt']},
    include_package_data=True,
    python_requires='>=3.5',
    install_requires=[
        'lxml>=4',
        'networkx>=2',
        'pint>=0.8.1',
        'rdflib>=4',
        'sympy>=1.4',
    ],
    extras_require={
        'test': [
            'codecov',
            'flake8',
            'isort',
            'pytest>=3.2',
            'pytest-cov',
        ],
    },
)
