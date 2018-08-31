import os
from time import gmtime, strftime
from glob import glob

import _fit_evaluators
from _dataPath import DataPath

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
    if name not in fits_collection.keys():
        raise Exception('Invalid fit name : %s'%name)
    else:
        fit = fits_collection[name].fit_class(name.split('surfinBH')[-1])
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

    print 'Downloading %s data'%name
    data_url = fits_collection[name].data_url
    os.system('mkdir -p {0}; cd {0}; curl -# -L -O {1}'.format(data_dir,
            data_url))


##############################################################################

#### Add fits here, each name should start with surfinBH.
fits_collection = {}

fits_collection['surfinBH3dq8'] = FitAttributes( \
    fit_class = _fit_evaluators.Fit3dq8,
    desc = 'Fits for remnant mass, spin and kick veclocity for nonprecessing'
        ' BBH systems.',
    data_url = 'https://www.dropbox.com/s/06mrxalxqjhzy9d/fit_3dq8.h5',
    refs = 'arxiv.2018.xxxx',
    )

fits_collection['surfinBH7dq2'] = FitAttributes( \
    fit_class = _fit_evaluators.Fit7dq2,
    desc = 'Fits for remnant mass, spin and kick veclocity for genrically'
        ' precessing BBH systems.',
    data_url = 'https://www.dropbox.com/s/8b4o7n5aswrnami/fit_7dq2.h5',
    refs = 'arxiv.2018.xxxx',
    )
