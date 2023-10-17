import pickle
import tempfile
from Bio import SeqIO
import os
import sys
import subprocess
import multiprocessing
from joblib import Parallel, delayed
import shutil

def process_pickle(pickle_file):

    with open(pickle_file, 'rb') as f:
      data = pickle.load(f)
      all_scores = [ data['ptm'], data['iptm'], data['pitm']['score'], data['interface']['score']]

    return all_scores

def collect_scores(folder):

  scores = dict()
  pkls = [ p for p in os.listdir(f'{folder}/{folder}/') if p.endswith('.pkl') and p.startswith('model') ]

  for pkl in pkls:
    p  = pkl.split("/")[-1].replace('.pkl', '.pdb')
    scores[ p ] = process_pickle(pkl) 

  with open(f'{folder}/scores.out', 'w') as o:
      o.write('PDB\tpTM\tipTM\tpiTM\ti-score\n')
      for k in sorted( scores.keys() ):
        s = "\t".join([ str(x) for x in scores[k] ])
        o.write(f'{k}\t{s}\n')
  print(f'{folder} is processed')


def check_outputs(out_dir):

  pdbs = sorted( [ p for p in os.listdir(out_dir) if p.endswith('.pdb') and p.startswith('model') ] )
  pkls = sorted( [ p for p in os.listdir(out_dir) if p.endswith('.pkl') and p.startswith('model') ] )

  pdb_prefix = sorted( [ p.split('.pdb')[0] for p in pdbs])
  pkl_prefix = sorted( [ p.split('.pkl')[0] for p in pkls])

  correct = [ f'model_{i}_multimer_v3' for i in range(1, 6) ]

  if pdb_prefix == correct and pkl_prefix == correct:
    return True

  else: 
    print(f'{out_dir}: Prediction needs to be completed')
    return False

def check_file(filename):
    with open(filename, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

        # Check the number of lines
        if len(lines) != 5:
            return False

        # Check the number of columns in each line
        for line in lines:
            # Split the line by whitespace (assuming columns are separated by whitespace)
            columns = line.split()

            # Check the number of columns in the line
            if len(columns) != 5:
                return False

    return True


def main():
  num_process = int(multiprocessing.cpu_count()/4*3)

  os.chdir(f'structure')
  folders_initial = [ f for f in os.listdir() if f.startswith('AT') ]
  folders = []

  for f in folders_initial:
    out_dir = os.path.join(f, f)
    score_out = os.path.join(f, 'scores.out')
    if check_outputs(out_dir):
        if os.path.exists(score_out) and os.path.exists(os.path.join(f, 'alphafold.done')):
            if not check_file(score_out):
                folders.append(f)
        else:
                folders.append(f)

  print(f'{len(folders)} will be processed')

  Parallel(n_jobs=num_process, timeout=300)(delayed(collect_scores)(f) for f in folders)


if __name__ == '__main__':
    main()
