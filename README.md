# mrspoc

[![Build Status](https://travis-ci.org/bmorris3/mrspoc.svg?branch=master)](https://travis-ci.org/bmorris3/mrspoc)
[![Documentation Status](https://readthedocs.org/projects/mrspoc/badge/?version=latest)](http://mrspoc.readthedocs.io/en/latest/?badge=latest) [![DOI](https://zenodo.org/badge/106725406.svg)](https://zenodo.org/badge/latestdoi/106725406)

M-dwarf Rotational Stellar Photocenter Offset Calculator (mrspoc).


### Installation

Clone this repository:
```bash
git clone https://github.com/bmorris3/mrspoc.git
```
Run this command to download a local copy of the required TGAS data, Hipparcos, 
and  interferometric radius tables:
```
cd mrspoc
python download_data.py
```
You can delete the data archive `data.tar.gz`.

Create an environment variable that contains the path to the `data/` directory, 
called `MRSPOC_DATA_DIR`. For example, add the following line to your `~/.cshrc`:
```
setenv MRSPOC_DATA_DIR '/path/to/mrspoc/data/'
```
or `~/.bashrc`: 
```
MRSPOC_DATA_DIR='/path/to/mrspoc/data/'
```
Change directories into the repository, and run: 
```bash
python setup.py install
```
And you're good to go! You can run the unit tests with: 
```bash
python setup.py test
```

Citation
--------
If you make use of this code, please cite [Morris et al 2018](http://adsabs.harvard.edu/abs/2018arXiv180209943M):
```
@ARTICLE{Morris2018,
   author = {{Morris}, B.~M. and {Agol}, E. and {Davenport}, J.~R.~A. and 
	{Hawley}, S.~L.},
    title = "{Spotting stellar activity cycles in Gaia astrometry}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1802.09943},
 primaryClass = "astro-ph.SR",
 keywords = {astrometry, Sun: activity, stars: activity, stars: individual: GJ 1243, stars: individual: KIC 7174505, stars: individual: AX Mic},
     year = 2018,
    month = jun,
   volume = 476,
    pages = {5408-5416},
      doi = {10.1093/mnras/sty568},
   adsurl = {http://adsabs.harvard.edu/abs/2018MNRAS.476.5408M},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
