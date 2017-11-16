# cellmlmanip
CellML loading and model equation manipulation

## Developer installation

Set up a virtual environment, and install requirements:
```sh
conda create -n cellmlmanip python=3.6
source activate cellmlmanip
pip install -r requirements/dev.txt
pip install -e .
```

## Testing

To run tests, just run
```sh
pytest
```
