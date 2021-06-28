from utils import Constellation, starlink_constellation
import numpy as np
import glob


def read_ss_objs(filename):
    result = np.genfromtxt(filename, names=True, skip_header=14)
    return result


def run_check(path='baseline_nexp2_v1.7_10yrs'):

    files = glob.glob(path+'/*.txt')
    arrays = [read_ss_objs(filename) for filename in files]
    ss_data = np.hstack(arrays)
    # XXX---starting with small constellation for speed to start
    sats = starlink_constellation(fivek=True, alt_limit=20.)#(supersize=True)
    const = Constellation(sats)
    

    # I think I want to output mjd, ra, dec for the satellites that get hit.


if __name__ == '__main__':
    run_check()


    


