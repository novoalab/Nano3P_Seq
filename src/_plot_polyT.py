#!/usr/bin/env python
import sys
import h5py
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use ('seaborn-whitegrid')
import ont_fast5_api
from ont_fast5_api.fast5_interface import get_fast5_file
import numpy as np
import pandas as pd
import os 
import time 


desc = '''
Author: Huanle Liu
Date: April 08 2021 
Place: Mercat de Ninot, BCN

Usage:
input 1: csv file from tailrFinder
input 2: single- | multi- fast5 file
output directory: folder name to output plots
output: pdf files (one for each read) with plots of squiggle with polyT marked
'''

if len (sys.argv) != 4:
    print (desc)
    sys.exit()

def raw_from_f5 (fast5handle, readid, f5type):
    if f5type == 'multi':
        read = fast5handle.get_read(readid)
        raw = read.get_raw_data()
        return raw 
    else:
        read = fast5handle.get_read (readid)
        squiggle = read.get_raw_data()
        return squiggle

def f5_type (f5handle):
    return 'multi' if len (f5handle.get_read_ids()) > 1 else 'single'

def main():
    df = pd.read_csv (sys.argv[1])
    print ('there are', df.shape[0],'reads')
    df = df[(df['tail_length'] > 0) & (df['tail_is_valid'] ) ]
    df.drop(columns=['read_type', 'tail_is_valid','file_path'], inplace=True)
    print (df.shape[0], 'has valid tails')
    indexes = df.index.to_numpy()
    f5 = get_fast5_file (sys.argv[2], mode='r')
    f5type = f5_type (f5)

    outdir = sys.argv[3]
    if not os.path.exists(outdir):
        os.mkdir (outdir)
    else:
        timestr = time.strftime("%Y%m%d-%H%M%S")
        outdir = 'tmp_tail_plots_{}'.format(timestr)
        os.mkdir (outdir)

    cnt = 0 

    for i in indexes:
        row = df.loc[i]
        readid = df.loc[i]['read_id']
        sigs = raw_from_f5 (f5, readid, f5type)
        outpdf = "{}/{}.pdf".format(outdir, readid)
        
        try:

            fig, ax = plt.subplots(figsize=(12,4))
            start, end = row['tail_start'], row['tail_end']
            ax.plot(sigs[:int(end)+4000], '-k')
            seg = sigs[int(start):int(end)]
            ax.plot(range(int(start), int(end)), seg, '-r')
            title = "Poly(T) | Tail Length [nt]:{} | Tail start: {} | Tail end: {} | Sample per nt: {}".format(
            row['tail_length'], row['tail_start'], row['tail_end'], row['samples_per_nt'])
            ax.set( xlabel='Sample index', ylabel='pA')
            fig.suptitle(title, fontsize=13)
            plt.axis('tight')
            fig.savefig(outpdf)
            cnt += 1 
        except:
            print ("plot", readid, "failed!")
    print ('in the end', cnt, "reads were plotted!")
if __name__ == '__main__':
    main()









