import os
from time import gmtime, strftime
from glob import glob


#-------------------------------------------------------------------------
def _data_path():
    """ Return the default path for fit data h5 files"""
    return os.path.abspath('%s/../data'%(os.path.dirname( \
        os.path.realpath(__file__))))

#-------------------------------------------------------------------------
def DownloadData(name):
    """ Downloads fit data to surfinBH/data diretory.
    """
    if name not in fits_collection.keys():
        raise Exception('Invalid fit name : %s'%name)

    data_url = fits_collection[name].url
    fname = os.path.basename(data_url)
    data_dir = _data_path()

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

    os.system('wget -q --show-progress --directory-prefix=%s %s'%(data_dir, \
        data_url))


#=============================================================================
class FitAttributes(object):
    """ Saves attributes of a particular fit.
    """
    def __init__(self, **kwargs):
        self.desc = kwargs['desc']
        self.url =  kwargs['url']
        self.refs =  kwargs['refs']


#### Add fits here
fits_collection = {}

fits_collection['3dq8'] = FitAttributes( \
    desc = 'Fits for remnant mass, spin and kick veclocity for nonprecessing'
        ' BBH systems.',
    url = 'https://www.dropbox.com/s/046g2eqvabjjtyl/fit_3dq8.h5',
    refs = 'Varma:2018_inprep',
    )

fits_collection['7dq2'] = FitAttributes( \
    desc = 'Fits for remnant mass, spin and kick veclocity for genrically'
        ' precessing BBH systems.',
    url = 'https://www.dropbox.com/s/np1rh9ijmdnu9ko/fit_7dq2.h5',
    refs = 'Varma:2018_inprep',
    )
