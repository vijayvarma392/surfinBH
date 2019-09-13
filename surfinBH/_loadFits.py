import os
import errno
try:
    from urllib.request import urlretrieve # py 3
except ImportError:
    from urllib import  urlretrieve # py 2
from time import gmtime, strftime
from glob import glob

from . import _fit_evaluators
from ._dataPath import DataPath

#=============================================================================
class FitAttributes(object):
    """ Saves attributes of a particular fit.
    """
    #-------------------------------------------------------------------------
    def __init__(self, **kwargs):
        self.fit_class =  kwargs['fit_class']
        self.desc = kwargs['desc']
        self.data_url =  kwargs['data_url']
        self.refs =  kwargs['refs']

#-------------------------------------------------------------------------
def LoadFits(name):
    """ Loads data for a fit.
    If data is not available, downloads it before loading.
    """
    if name not in fits_collection.keys():
        raise Exception('Invalid fit name : %s'%name)
    else:
        testPath = DataPath() + '/' + fits_collection[name].data_url.split('/')[-1]
        if (not os.path.isfile(testPath)):
            DownloadData(name)

        # allow for both naming formats surfinBH7dq2 and NRSur7dq4Remnant
        if 'surfinBH' in name:
            name_tag = name.split('surfinBH')[-1]
        else:
            name_tag = name.split('NRSur')[-1].split('Remnant')[0]

        fit = fits_collection[name].fit_class(name_tag)
        print('Loaded %s fit.'%name)
        return fit

#-------------------------------------------------------------------------
def DownloadData(name='all', data_dir=DataPath()):
    """ Downloads fit data to DataPath() diretory.
        If name='all', gets all fit data.
    """
    if name == 'all':
        for tmp_name in fits_collection.keys():
            DownloadData(name=tmp_name, data_dir=data_dir)
        return

    if name not in fits_collection.keys():
        raise Exception('Invalid fit name : %s'%name)

    print('Downloading %s data'%name)
    data_url = fits_collection[name].data_url
    filename = data_url.split('/')[-1]
    try:
        os.makedirs(data_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(data_dir):
            pass
        else:
            raise
    urlretrieve(data_url, data_dir + '/' + filename)


##############################################################################

#### Add fits here, each name should start with surfinBH.
fits_collection = {}

fits_collection['NRSur3dq8Remnant'] = FitAttributes( \
    fit_class = _fit_evaluators.Fit3dq8,
    desc = 'Fits for remnant mass, spin and kick veclocity for nonprecessing'
        ' BBH systems. This model was called surfinBH3dq8 in the paper.',
    data_url = 'https://zenodo.org/record/1435751/files/remnant_fits/fit_3dq8.h5',
    refs = 'arxiv:1809.09125',
    )

fits_collection['surfinBH7dq2'] = FitAttributes( \
    fit_class = _fit_evaluators.Fit7dq2,
    desc = 'Fits for remnant mass, spin and kick veclocity for generically'
        ' precessing BBH systems.',
    data_url = 'https://zenodo.org/record/1435751/files/remnant_fits/fit_7dq2.h5',
    refs = 'arxiv:1809.09125',
    )


fits_collection['NRSur7dq4Remnant'] = FitAttributes( \
    fit_class = _fit_evaluators.Fit7dq4,
    desc = 'Fits for remnant mass, spin and kick veclocity for generically'
        ' precessing BBH systems up to mass ratio 4.',
    data_url = 'http://www.tapir.caltech.edu/~vvarma/fit_7dq4.h5',
    refs = 'arxiv:1905.09300',
    )
