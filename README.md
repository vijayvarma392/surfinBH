# Welcome to surfinBH!

<img src="https://raw.githubusercontent.com/vijayvarma392/surfinBH/master/images/point_break.jpeg" alt="Point Break" width="400px"/>


<br/>
<br/>

_**surfinBH**_ provides _**sur**rogate **fin**al **B**lack_ _**H**ole_
properties for mergers of binary black holes (BBH). Just like Point Break, but
with black holes! This package lives on
[GitHub](https://github.com/vijayvarma392/surfinBH).

These fits are described in the following papers:

[1] Vijay Varma, Davide Gerosa, Francois Hebert and Leo C. Stein, 2018, in
preparation.

If you find this package useful in your work, please cite reference [1] and,
if available, the relevant paper describing the particular model.

## Installation

### PyPi
_**surfinBH**_ is available through [PyPi](https://pypi.org/project/surfinBH/).

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
* numpy
* scipy
* scikit-learn (at least 0.19.1)
* h5py

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
>>> fit_7dq2.h5  100%[======================>]  42.85M  495KB/s  in 60s
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

* [surfinBH3dq8](http://htmlpreview.github.io/?https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_3dq8.html)

* [surfinBH7dq2](http://htmlpreview.github.io/?https://github.com/vijayvarma392/surfinBH/blob/master/examples/example_7dq2.html)


## Credits
The code is developed and maintained by [Vijay Varma](http://www.tapir.caltech.edu/~vvarma/). Please, report bugs to
[&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;](mailto:&#118;&#118;&#097;&#114;&#109;&#097;&#064;&#099;&#097;&#108;&#116;&#101;&#099;&#104;&#046;&#101;&#100;&#117;).
