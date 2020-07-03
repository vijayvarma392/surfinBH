[![github](https://img.shields.io/badge/GitHub-surfinBH-blue.svg)](https://github.com/vijayvarma392/surfinBH)
[![PyPI version](https://badge.fury.io/py/surfinBH.svg)](https://pypi.org/project/surfinBH/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/surfinbh.svg)](https://anaconda.org/conda-forge/surfinbh)
[![DOI](https://zenodo.org/badge/145179417.svg)](https://zenodo.org/badge/latestdoi/145179417)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/vijayvarma392/surfinBH/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/vijayvarma392/surfinBH.svg?branch=master)](https://travis-ci.org/vijayvarma392/surfinBH)

# Welcome to surfinBH!

<img src="https://raw.githubusercontent.com/vijayvarma392/surfinBH/master/images/interaction.jpeg" alt="BHScattering" width="600px"/>

<br/>

_**surfinBH**_ provides _**sur**rogate **fin**al **B**lack_ _**H**ole_
properties for mergers of binary black holes (BBH).

These fits are described in the following papers: <br/>
[1] Vijay Varma, D. Gerosa, L. C. Stein, F. Hébert and H. Zhang,
[arxiv:1809.09125](https://arxiv.org/abs/1809.09125).

[2] Vijay Varma, S. E. Field, M. A. Scheel, J. Blackman, D. Gerosa, L. C.
Stein, L. E. Kidder, H. P. Pfeiffer,
[arxiv:1905.09300](https://arxiv.org/abs/1905.09300).

If you find this package useful in your work, please cite reference [1] and,
if available, the relevant paper describing the particular model. Please also
cite this package, see the DOI badge at the top of this page for BibTeX keys.

This package is compatible with both python2 and python3.
This package lives on [GitHub](https://github.com/vijayvarma392/surfinBH) and
is tested every day with [Travis CI](https://travis-ci.org/). You can see the
current build status of the master branch at the top of this page.

## Installation

### PyPI
_**surfinBH**_ is available through [PyPI](https://pypi.org/project/surfinBH/):

```shell
pip install surfinBH
```

### Conda
_**surfinBH**_ is available on [conda-forge](https://anaconda.org/conda-forge/surfinbh):

```shell
conda install -c conda-forge surfinbh
```


### From source

```shell
git clone https://github.com/vijayvarma392/surfinBH
cd surfinBH
git submodule init
git submodule update
python setup.py install
```

If you do not have root permissions, replace the last step with
`python setup.py install --user`


## Dependencies
All of these can be installed through pip or conda.
* [numpy](https://docs.scipy.org/doc/numpy/user/install.html)
* [scipy](https://www.scipy.org/install.html)
* [h5py](http://docs.h5py.org/en/latest/build.html)
* [scikit-learn](http://scikit-learn.org/stable/install.html) (at least 0.19.1)
* [lalsuite](https://pypi.org/project/lalsuite) (at least 6.70)
* [gwsurrogate](https://pypi.org/project/gwsurrogate)
* [NRSur7dq2](https://pypi.org/project/NRSur7dq2) (only for surfinBH7dq2)

## Usage

```python
import surfinBH
```

### See list of available fits
```python
print(surfinBH.fits_collection.keys())
>>> ['NRSur3dq8Remnant', 'NRSur7dq4Remnant', 'surfinBH7dq2']
```

Pick your favorite fit and get some basic information about it.
```python
fit_name = 'NRSur7dq4Remnant'

surfinBH.fits_collection[fit_name].desc
>>> 'Fits for remnant mass, spin and kick veclocity for generically precessing BBH systems up to mass ratio 4.'

surfinBH.fits_collection[fit_name].refs
>>> 'arxiv:1905.09300'
```

### Load the fit
This only needs to be done **once** at the start of your script.
If the fit data is not already downloaded, this will also download the data.
```python
fit = surfinBH.LoadFits(fit_name)
>>> Loaded NRSur7dq4Remnant fit.
```
### Evaluation
The evaluation of each fit is different, so be sure to read the documentation.
This also describes the frames in which different quantities are defined.
```python
help(fit)
```

We also provide ipython examples for usage of different fits:

##### Current fits

* [NRSur7dq4Remnant](https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_7dq4.ipynb) (Ref. [2])

* [NRSur3dq8Remnant](https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_3dq8.ipynb) (called surfinBH3dq8 in Ref. [1])

##### Older fits

* [surfinBH7dq2](https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_7dq2.ipynb) (Ref. [1]. Superseded by NRSur7dq4Remnant)

## Animations

We also provide a tool to visualize the binary black hole scattering process,
see
[binary black hole explorer](https://vijayvarma392.github.io/binaryBHexp/).
Here's an example:

<img src="https://raw.githubusercontent.com/vijayvarma392/binaryBHexp/master/animations/video.gif" width="500"/>


## Making contributions
See this
[README](https://github.com/vijayvarma392/surfinBH/blob/master/README_developers.md)
for instructions on how to make contributions to this package.

You can find the list of contributors
[here](https://github.com/vijayvarma392/surfinBH/graphs/contributors).


## Credits
The code is developed and maintained by [Vijay
Varma](https://vijayvarma.com). Please report bugs by raising an
issue on our [GitHub](https://github.com/vijayvarma392/surfinBH) repository.
