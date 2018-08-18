# Welcome to surfinBH!

<img src="https://raw.githubusercontent.com/vijayvarma392/surfinBH/master/images/point_break.jpeg" alt="Point Break" width="400px"/>


<br/>
<br/>

*surfinBH* provides fits for *Surrogate Final Black Hole* properties from mergers of binary black holes. Just like Point Break, but with black holes! This
package lives on [GitHub](https://github.com/vijayvarma392/surfinBH).

These fits are described in the following papers:

[1] Vijay Varma, Davide Gerosa, Francois Hebert and Leo C. Stein, in
preparation.

If you find this package useful in your work, please cite reference [1] and,
if available, the relevant paper describing the particular model.

## Installation

### PyPi
*surfinBH* is available through [PyPi](https://pypi.org/project/surfinBH/).

```shell
pip install surfinBH
```


### From source

```shell
git clone https://github.com/vijayvarma392/surfinBH
cd surfinBH
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
print(surfinBH.FIT_CLASSES.keys())
>>> ['7dq2', '3dq8']
```

