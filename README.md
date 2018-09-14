[![github](https://img.shields.io/badge/GitHub-surfinBH-blue.svg)](https://github.com/vijayvarma392/surfinBH)
[![DOI](https://zenodo.org/badge/145179417.svg)](https://zenodo.org/badge/latestdoi/145179417)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/vijayvarma392/surfinBH/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/vijayvarma392/surfinBH.svg?branch=master)](https://travis-ci.org/vijayvarma392/surfinBH)

# Welcome to surfinBH!

<img src="https://raw.githubusercontent.com/vijayvarma392/surfinBH/master/images/point_break.jpeg" alt="Point Break" width="400px"/>


<br/>
<br/>

_**surfinBH**_ provides _**sur**rogate **fin**al **B**lack_ _**H**ole_
properties for mergers of binary black holes (BBH).

These fits are described in the following papers: <br/>
[1] Vijay Varma, Davide Gerosa, François Hébert, Leo C. Stein and Hao Zhang,
2018, in preparation.

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
* [lalsuite](https://pypi.org/project/lalsuite)
* [NRSur7dq2](https://pypi.org/project/NRSur7dq2) (at least 1.0.5)

## Usage

```python
import surfinBH
```

### See list of available fits
```python
print(surfinBH.fits_collection.keys())
>>> ['surfinBH3dq8', 'surfinBH7dq2']
```

Pick your favorite fit and get some basic information about it.
```python
fit_name = 'surfinBH7dq2'

surfinBH.fits_collection[fit_name].desc
>>> 'Fits for remnant mass, spin and kick veclocity for generically precessing BBH systems.'

surfinBH.fits_collection[fit_name].refs
>>> 'arxiv.2018.xxxx'
```

### Load the fit
This only needs to be done **once** at the start of your script.
```python
fit = surfinBH.LoadFits(fit_name)
>>> Loaded surfinBH7dq2 fit.
```
### Evaluation
The evaluation of each fit is different, so be sure to read the documentation.
This also describes the frames in which different quantities are defined.
```python
help(fit)
```

We also provide ipython examples for usage of different fits:

* [surfinBH3dq8](https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_3dq8.ipynb)

* [surfinBH7dq2](https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_7dq2.ipynb)

## Making contributions
See this
[README](https://github.com/vijayvarma392/surfinBH/blob/master/README_developers.md)
for instructions on how to make contributions to this package.

You can find the list of contributors
[here](https://github.com/vijayvarma392/surfinBH/graphs/contributors).


## Credits
The code is developed and maintained by [Vijay Varma](http://www.tapir.caltech.edu/~vvarma/). Please, report bugs to
[&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;](mailto:&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;).
