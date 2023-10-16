import multiprocessing
import subprocess
import tempfile
import os, sys
from pathlib import Path
from joblib import Parallel, delayed
import traceback
import shutil
import logging

# Paths to databases
uniref90_database_path = Path('/global/scratch/users/skyungyong/Software/alphafold/Database/uniref90/uniref90.fasta')
mgnify_database_path = Path('/global/scratch/users/skyungyong/Software/alphafold/Database/mgnify/mgy_clusters.fa')
uniprot_database_path = Path('/global/scratch/users/skyungyong/Software/alphafold/Database/uniprot/uniprot.fasta')
pdb_database_path = Path('/global/scratch/users/skyungyong/Software/AF2-DB_2023-03-31/pdb_seqres.txt')

# Set Jackhmmer parameters
jackhmmer_cmd_flags = ['-o', '/dev/null',
                      '--noali',
                      '--F1', str(0.0005),
                      '--F2', str(0.00005),
                      '--F3', str(0.0000005),
                      '--incE', str(0.0001),
                      '-E', str(0.0001),
                      '--cpu', '1',
                      '-N', str(1)]

# Set HMMER parameters
hmmsearch_cmd_flags = ['--F1', '0.1',
                       '--F2', '0.1',
                       '--F3', '0.1',
                       '--incE', '100',
                       '-E', '100',
                       '--domE', '100',
                       '--incdomE', '100',
                       '--noali',
                       '--cpu', '1']

def run_jackhmmer(input_file, database_path, output_path):
    """Run jackhmmer and return the output in sto format"""
    with tempfile.TemporaryDirectory(dir='/tmp/') as tmpdir:
      sto_path = Path(tmpdir) / 'output.sto'
      cmd = ['jackhmmer'] + jackhmmer_cmd_flags + ['-A', str(sto_path)] + [str(input_file), str(database_path)]

      logging.info("Running command: %s", " ".join(cmd))
      process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

      if process.returncode == 0:
        shutil.copy(sto_path, output_path)
      else:
        logging.error(f"Failed: {' '.join(cmd)}")
        logging.error(process.stderr.decode())

def run_hmmsearch(database_path, output_path, msa_output_dir):
    """Run hmmsearch and return the output in sto format"""
    # Build hmm profile with uniref90_hits
    with tempfile.TemporaryDirectory(dir='/tmp/') as tmpdir:
      query_file = Path(tmpdir) / 'query.hmm'
      subprocess.run(['hmmbuild', str(query_file), f'{msa_output_dir}/uniref90_hits.sto'])

      # Run hmmsearch
      sto_path = Path(tmpdir) / 'output.sto'
      cmd = ['hmmsearch'] + hmmsearch_cmd_flags + ['-A', str(sto_path), str(query_file), str(database_path)]

      logging.info("Running command: %s", " ".join(cmd))
      process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

      if process.returncode == 0:
        shutil.copy(sto_path, output_path)
      else:
        logging.error(f"Failed: {' '.join(cmd)}")
        logging.error(process.stderr.decode())

def compute_msa(input_fasta_path, msa_output_dir, mode):

    """Run jackhmmer and hhsearch and return outputs"""
    # Run jackhmmer against uniref
    uniref90_out_path = Path(msa_output_dir) / 'uniref90_hits.sto'
    if not uniref90_out_path.exists():
        run_jackhmmer(input_fasta_path, uniref90_database_path, uniref90_out_path)

    # Run jackhmmer against mgnify
    mgnify_out_path = Path(msa_output_dir) / 'mgnify_hits.sto'
    if not mgnify_out_path.exists():
        run_jackhmmer(input_fasta_path, mgnify_database_path, mgnify_out_path)

    if mode == 'Multimer':
      # Run jackhmmer against uniprot
      uniprot_out_path = Path(msa_output_dir) / 'uniprot_hits.sto'
      if not uniprot_out_path.exists():
        run_jackhmmer(input_fasta_path, uniprot_database_path, uniprot_out_path)

      # Run hmmer against seqres
      pdb_out_path = Path(msa_output_dir) / 'pdb_hits.sto'
      if not pdb_out_path.exists() and uniref90_out_path.exists():
        run_hmmsearch(pdb_database_path, pdb_out_path, msa_output_dir)


def msa_parallel(sequence, mode):

    seq_path = Path(os.path.abspath(sequence))
    fa_path  = str(seq_path) + '/' + f'{sequence}.fasta'
    running_file = seq_path / "alphafold.msa.running"
    done_file = seq_path / "alphafold.msa.done"

    if not (running_file.exists() or done_file.exists()):
        with open(running_file, "w"): pass

        msa_output_dir = seq_path / sequence / 'msas'
        msa_output_dir.mkdir(parents=True, exist_ok=True)

        compute_msa(fa_path, msa_output_dir, mode)

        with open(done_file, "w"): pass
        os.unlink(running_file)

def find_sequences(prefix):
    folders = list(Path(".").glob(f"{prefix}*"))
    return [str(f) for f in folders if os.path.isdir(f)]

def main():
    try:
        prefix = sys.argv[1]

        # Default == Multimer
        # Enable extra database searches
        if len(sys.argv) == 3: mode = sys.argv[2]
        else: mode = 'Multimer'

        if "," in prefix:
            sequences = []
            for p in prefix[:-1].split(","):
                sequences += find_sequences(p)
        else:
            sequences = find_sequences(prefix)

        if not sequences:
            print(f"No folders found for prefix '{prefix}'.")
            return

        sequences = set(sequences)
        num_cores = multiprocessing.cpu_count() # Reduce for longer sequences or those that have a lot of matches
        Parallel(n_jobs=num_cores)(delayed(msa_parallel)(item, mode) for item in sequences)

    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    main()
