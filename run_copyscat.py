import subprocess
import pathlib
import os
import sys, getopt
import glob
import numpy as np

import chromograph
from chromograph.pipeline import config
config = config.load_config()

import logging

logger = logging.getLogger()
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%H:%M:%S')

binsize = 1000000
min_filter = 5000
chrom_sizefile = os.path.join(chromograph.__path__[0], f'references/male.GRCh38.chrom.sizes')
try:
   opts, args = getopt.getopt(sys.argv[1:],"hs:b:f:g:",["help","sample=","binsize=","frags=","genome="])
   #opts, args = getopt.getopt(sys.argv[1:],"h",["help"])
   print("process_fragment_file",len(opts))
except getopt.GetoptError as err:
   print('process_fragment_file.py -s <sample> -b <binsize> -f <minfrags> -g <genomefile.chrom.sizes>')
   print(str(err),"process_fragment_file")
   sys.exit(2)
for opt, arg in opts:
   print("process_fragment_file",opt)
   if opt in ("-h","help"):
    print('process_fragment_file.py  -s <sample> -b <binsize> -f <minfrags> -g <genomefile>')
    sys.exit()
   elif opt in ("-s", "--sample"):
     sample = arg
   elif opt in ("-b", "--binsize"):
     binsize = int(arg)
   elif opt in ("-f", "--frags"):
     min_filter = int(arg)
   elif opt in ("-g", "--genome"):
     chrom_sizefile = arg

logging.info(f'Looking for {sample} Cellranger outputs at {config.paths.cell_ranger}')
dirs = glob.glob(f'{config.paths.cell_ranger}/{sample}*')

## Find the right 10X folder
if len(dirs) > 0:
    dirs = np.array([d for d in dirs if os.path.isdir(d)])
    dirs = dirs[[len(d.split('_')) == 4 for d in dirs]]
    n_flowcells = [len(d.split('_')[-2]) for d in dirs]
    dirs = np.array(dirs)[np.where(n_flowcells==np.max(n_flowcells))[0]]
    ID = [int(d.split('_')[-1]) for d in dirs]
    
    indir = dirs[np.where(ID==np.max(ID))[0]][0]
    inputfile = os.path.join(indir, 'outs', 'fragments.tsv.gz')
    meta_files = [os.path.join(indir, 'outs', x) for x in ['singlecell.csv', 'per_barcode_metrics.csv']]
    meta_file = [x for x in meta_files if os.path.isfile(x)][0]

    outdir = f"{config.paths.samples}/{sample}"
    outfile = os.path.join(outdir, 'fragments.tsv')
    logging.info(f'Using cellranger output: {indir}')
    logging.info(f'Saving to {outdir}')

else:
    logging.info(f'Could not find sample')

path_copyscat = config.paths.copyscat
if not os.path.isdir(outdir):
    os.mkdir(outdir)

## Process fragments
subprocess.run(['python', os.path.join(path_copyscat, 'process_fragment_file.py'), 'i', inputfile, 'o', outfile, 'b', binsize, 'f', min_filter, 'g', chrom_sizefile, 'm', meta_file])

## Find neoplastic cells
subprocess.run([config.paths.R, os.path.join(path_copyscat, 'find_clones.R'), outdir])
