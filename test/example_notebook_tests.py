import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import platform
import surfinBH

python_version = platform.python_version()
if python_version[:2] == '2.':
    python_version = 'python2'
elif python_version[:2] == '3.':
    python_version = 'python3'

def test_example_notebooks():
    """ Tests that all the example ipython notebooks of format
    surfinBH/examples/example_*.ipynb are working. Since we expect these to be
    used by our users, it would be emabarassing if our own examples failed.
    """
    # List of all available fits
    fit_names = surfinBH.fits_collection.keys()

    for name in fit_names:
        # allow for both naming formats surfinBH7dq2 and NRSur7dq4Remnant
        if 'surfinBH' in name:
            name_tag = name.split('surfinBH')[-1]
        else:
            name_tag = name.split('NRSur')[-1].split('Remnant')[0]

        notebook_filename = 'examples/example_%s.ipynb'%name_tag
        print('testing %s'%notebook_filename)

        with open(notebook_filename) as f:
            nb = nbformat.read(f, as_version=4)

        ep = ExecutePreprocessor(timeout=None, kernel_name=python_version)

        ep.preprocess(nb, {'metadata': {'path': '.'}})
