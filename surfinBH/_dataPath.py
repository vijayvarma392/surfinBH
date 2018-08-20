import os

#-------------------------------------------------------------------------
def DataPath():
    """ Return the default path for fit data h5 files"""
    return os.path.abspath('%s/../surfinBH_data'%(os.path.dirname( \
        os.path.realpath(__file__))))
