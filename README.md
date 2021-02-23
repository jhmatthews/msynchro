## msynchro

Routines for calculating synchrotron emission from a population of electrons (under certain assumptions). Also includes a tridiagonal matrix (TDM) algorithm for evolving particle distributions. The TDM algorithm is compiled as a C extension. 

### Docs 

Documentation is hosted on [ReadTheDocs](https://msynchro.readthedocs.io/).

[![Documentation Status](https://readthedocs.org/projects/msynchro/badge/?version=latest)](https://msynchro.readthedocs.io/en/latest/?badge=latest)

### Prerequisites

The code requires numpy, scipy and a working C compiler.

### Installation

Typing

```
python setup.py install
```

should install msynchro as a module. You can then import it as usual.

### Using the code

If you use this code, please cite the corresponding paper, [Matthews \& Taylor 2021]()
