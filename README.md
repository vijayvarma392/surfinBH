[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/vijayvarma392/surfinBH/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/vijayvarma392/surfinBH.svg?branch=master)](https://travis-ci.org/vijayvarma392/surfinBH)

# Welcome to surfinBH!

<img src="https://raw.githubusercontent.com/vijayvarma392/surfinBH/master/images/point_break.jpeg" alt="Point Break" width="400px"/>


<br/>
<br/>

_**surfinBH**_ provides _**sur**rogate **fin**al **B**lack_ _**H**ole_
properties for mergers of binary black holes (BBH). Just like Point Break, but
with black holes!

These fits are described in the following papers: <br/>
[1] Vijay Varma, Davide Gerosa, Francois Hebert and Leo C. Stein, 2018, in
preparation.

If you find this package useful in your work, please cite reference [1] and,
if available, the relevant paper describing the particular model.

This package lives on [GitHub](https://github.com/vijayvarma392/surfinBH) and
is tested every day with [Travis CI](https://travis-ci.org/). You can see the
current build status of the master branch at the top of this page.

## Installation

### PyPi
_**surfinBH**_ is available through [PyPI](https://pypi.org/project/surfinBH/).

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

## Usage

```python
import surfinBH
```

See list of available fits
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

Get data for the fit. This only needs to done **once, ever**.
```python
surfinBH.DownloadData(fit_name)
>>> ################################################################ 100.0%
```

Load the fit. This only needs to be done **once** at the start of your script.
```python
fit = surfinBH.LoadFits(fit_name)
>>> Loaded surfinBH7dq2 fit.
```

The evaluation of each fit is different, so be sure to read the documentation.
This also defines the frames in which different quantities are defined.
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


## Credits
The code is developed and maintained by [Vijay Varma](http://www.tapir.caltech.edu/~vvarma/). Please, report bugs to
[&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;](mailto:&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;).
