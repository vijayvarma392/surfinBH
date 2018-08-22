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
def DownloadData(name, data_dir=DataPath() ):
    """ Downloads fit data to DataPath() diretory.
    """
    if name not in fits_collection.keys():
        raise Exception('Invalid fit name : %s'%name)

    data_url = fits_collection[name].data_url
    fname = os.path.basename(data_url)

    # If file already exists, move it to backup dir with time stamp
    if os.path.isfile('%s/%s'%(data_dir, fname)):
        timestamp=strftime("%Y%b%d_%Hh:%Mm:%Ss", gmtime())
        backup_fname = '%s_%s'%(timestamp, fname)
        backup_dir = '%s/backup'%(data_dir)
        os.system('mkdir -p %s'%backup_dir)
        print('\n%s file exits, moving to %s/%s.'%(fname, backup_dir, \
            backup_fname))
        os.system('mv %s/%s %s/%s'%(data_dir, fname, backup_dir, backup_fname))
        number_of_backup_files = glob('%s/*_%s'%(backup_dir, fname))
        if len(number_of_backup_files) > 5:
            print('There are a lot of backup files in %s, consider removing'
                ' some.'%backup_dir)

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
