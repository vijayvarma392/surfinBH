"""surfinBH
========

Surrogate Final BH properties.

Example usage
-------------

To evaluate fits for remnant mass, spin and kick velocity of nonprecessing
binary black hole systems as presented in Varma:2018_in_prep.

import surfinBH

# First load the model
fit_3dq8 = surfinBH.LoadFits('3dq8')

# Mass ratio and component spins along orbital angular momentum direction
q = 6.7
chi1z = 0.74
chi2z = -0.6
x = [q, chi1z, chi2z]

## Evaluate the fits and GPR error estimate.

# Final mass
mC, mC_err_est = fit_3dq8('mC', x)

# Final spin along orbital angular momentum direction
chiCz, chiCz_err_est = fit_3dq8('chiC', x)

# Final kick velocity (velCz is zero for nonprecessing systems)
velC, velC_err_est = fit_3dq8('velC', x)
velCx, velCy = velC
velCx_err_est, velCy_err_est = velC_err_est

Additional examples can be found in the accompanying ipython notebooks.
To get a list of all available fits do:
surfinBH.FIT_CLASSES.keys()
"""

__copyright__ = "Copyright (C) 2018 Vijay Varma"
__email__ = "vvarma@caltech.edu"
__status__ = "testing"
__author__ = "Vijay Varma"
__license__ = """
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from _loadFits import LoadFits
from _loadFits import fits_collection
from _loadFits import DownloadData
