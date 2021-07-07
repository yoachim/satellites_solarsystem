import numpy as np
import glob
import pandas as pd



if __name__ == '__main__':

    filenames = glob.glob('survived_*')
    filenames.sort()

    array_list = []
    for fn in filenames:
        print(fn)
        array_list.append(pd.read_csv(fn))
    ss_data = pd.concat(array_list).to_csv('survived_merged.txt', index=False, sep=' ')
