# mrspoc

[![Build Status](https://travis-ci.org/bmorris3/mrspoc.svg?branch=master)](https://travis-ci.org/bmorris3/mrspoc)
[![Documentation Status](https://readthedocs.org/projects/mrspoc/badge/?version=latest)](http://mrspoc.readthedocs.io/en/latest/?badge=latest)

M-dwarf Rotational Stellar Photocenter Offset Calculator (mrspoc).


### Installation

Clone this repository:
```bash
git clone https://github.com/bmorris3/mrspoc.git
```
Create an environment variable that contains the path to the `data/` directory, 
called `MRSPOC_DATA_DIR`. For example, add the following line to your `~/.cshrc`:
```
setenv MRSPOC_DATA_DIR '/path/to/mrspoc/data/'
```
Change directories into the repository, and run: 
```bash
python setup.py install
```
And you're good to go! You can run the unit tests with: 
```bash
python setup.py test
```
