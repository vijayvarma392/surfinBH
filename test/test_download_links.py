import os, shutil
import surfinBH
import h5py

def test_download_links():
    """ Tests that the download links for fit data are working.
    """

    # dir to download data to
    out_dir = 'test/download_data'

    # remove out_dir if it already exists and make a new one
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.system('mkdir -p %s'%out_dir)

    # List of all available fits
    fit_names = surfinBH.fits_collection.keys()
    for name in fit_names:
        surfinBH.DownloadData(name=name, data_dir=out_dir)

        short_name = name.split('surfinBH')[-1]
        # check that it has the right name
        assert(os.path.isfile('%s/fit_%s.h5'%(out_dir, short_name)))
        # check that the fit_name matches with the name in the attributes
        # of h5 file.
        h5file = h5py.File('%s/fit_%s.h5'%(out_dir, short_name))
        assert(short_name == h5file.attrs['name'].decode('utf-8'))
        h5file.close()
