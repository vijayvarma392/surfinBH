# Instructions for developers

## Contributing to surfinBH

The preferred method of making contributions is to
[fork](https://help.github.com/articles/fork-a-repo/) + [pull
request](https://help.github.com/articles/about-pull-requests/) from the main
[repo](https://github.com/vijayvarma392/surfinBH).

After cloning your fork, do:

```shell
cd surfinBH
git submodule init
git submodule update
```

Before doing a pull request, you should check that your changes don't break
anything by running `py.test` from the root directory of your check-out. Every
pull request will be automatically tested by github.

## Adding a new fit

The fits in this package have the naming format: NRSur\*Remnant <br/> Let's say your
fancy new fit has `fit_name = '23dModGR'`, the name to load and evaluate the
package would be `'NRSur23dModGRRemnant'`.

You need to do the following to add this fit to this package:

* Add `fit_23dModGR.py` in `surfinBH/_fit_evaluators/`; see
  `surfinBH/_fit_evaluators/fit_3dq8.py` for an example.
* Add `from fit_23dModGR import Fit23dModGR` to
  `surfinBH/_fit_evaluators/__init__.py`.
* Add `'NRSur23dModGRRemnant'` key to `fits_collection` in `surfinBH/_loadFits.py`.
  See example for `surfinBH3dq8` in the same file (`surfinBH*` was the old format).
* Add `example_23dModGR.ipynb` in `examples`. See `examples/example_3dq8.ipynb`
  for an (cough, cough) example.
* Generate regression data for testing using `test/generate_regression_data.py`.
  Run it from the root directory of your check-out.
* Test (see above), commit and push all of the above changes.

Note: Do not push the fit data itself, but push the regression data generated
in `test/regression_data/`.

## PyPI release

Note: This is currently under Vijay's account, so only he can do this.
First make sure :

1. There are no data files in surfinBH/data. The
   required data files will get downloaded automatically.

2. `surfinBH/_eval_pysur` is not empty. If empty, do `git submodule update --init .`.

Then:

```shell
python setup.py sdist bdist_wheel
twine upload dist/*
```

NOTE for the uploader: For the last step, you should have the following in your `~/.pypirc`:

```
[distutils]
index-servers =
    pypi

[pypi]
username = __token__
```

When asked for the password, get a fresh API token https://pypi.org/manage/account/.
