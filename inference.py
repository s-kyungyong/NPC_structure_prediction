from pathlib import Path
import os
import subprocess
import shutil
import traceback
import re
import logging
import datetime
from Bio import SeqIO

# Input/output directories
target_dir = Path(os.getcwd())
model_out_dir = Path('/global/scratch/users/skyungyong/NPC/structure')

# Accessory scripts
feature = '/global/scratch/users/skyungyong/CO_Yangnan/run_fea_gen.sh'
af2complex = '/global/scratch/users/skyungyong/CO_Yangnan/af2complex.sh'

def log_message(message):
    """Prints a message with the current time to the console and log file."""
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_message = f"[{now}] {message}"
    print(log_message)

def search_pair():
    # Get all folders
    records = list(SeqIO.parse("NPC.all.fasta", "fasta"))
    records.sort(key=lambda x: len(x.seq))

    targets = [f for f in os.listdir(target_dir) if os.path.isdir(target_dir / f)]

    # Iterate each pair, pick pairs that satisfy seqence length limit
    pairs = dict()
    min_limit, max_limit = 1, 2000

    for i, record in enumerate(records):
        if len(record.seq) <= max_limit:
            for j in range(i + 1, len(records)):
                target_record = records[j]
                total_length = len(record.seq) + len(target_record.seq)

                if min_limit < total_length <= max_limit:
                    pairs[f'{record.id}_and_{target_record.id}'] = total_length
        else:
            break

    # A list of keys based on sorted length
    return sorted(pairs.items(), key=lambda x: x[1])

def inference(pair_info):

    # Get the pair information
    pair, pair_len = pair_info
    query, hit = pair.split('_and_')

    # Set input file directories
    query_dir = target_dir / query
    hit_dir = target_dir / hit

    # Set output directories
    out_path = model_out_dir / pair
    fa_path = out_path / f'{pair}.fasta'
    running_file = out_path / "alphafold.running"
    done_file = out_path / "alphafold.done"

    # Generate folders if missing
    msa_output_dir = out_path / pair / 'msas'
    if not msa_output_dir.exists():
        msa_output_dir.mkdir(parents=True, exist_ok=True)
        for l in ['A', 'B']:
            (msa_output_dir / l).mkdir(parents=True, exist_ok=True)

    # Check if process is already running in the directory
    if not (running_file.exists() or done_file.exists()):
        with open(running_file, "w") as f: pass
        log_message(f'Starting inference for pair {pair}')

        # Generate input
        input_list = out_path / 'input.list'
        with open(input_list, 'w') as f: f.write(f'A/B {pair_len} {pair}')

        # Create a fasta file
        with open(fa_path, 'w') as f:
            with open(query_dir / f'{query}.fasta', 'r') as r: f.write(r.read())
            with open(hit_dir / f'{hit}.fasta', 'r') as r: f.write(r.read())

        # Copy files
        shutil.copytree(str(query_dir / query / 'msas'), str(msa_output_dir / 'A'), dirs_exist_ok=True)
        shutil.copytree(str(hit_dir / hit / 'msas'), str(msa_output_dir / 'B'), dirs_exist_ok=True)
        log_message(f"Input files generated for {pair}")

        # Generate features
        feature_file = out_path / pair / 'features.pkl'
        if not feature_file.exists():
          subprocess.run([feature, str(fa_path), str(out_path)])
          logging.info(f'Generated features for pair {pair}')

        # Run af2complex
        subprocess.run([af2complex, str(input_list), str(out_path)])
        log_message(f'Generated predictions for pair {pair}')

        # Clean up the MSAs
        if msa_output_dir.exists():
          shutil.rmtree(msa_output_dir)

        # Finish the inference
        with open(done_file, "w") as f: pass
        if running_file.exists():
          running_file.unlink()
        log_message(f'Finished inference for pair {pair}')

def main():
    """Iterate over a list of complexes and model each one."""
    # Get list of complexes to model in the given GPU
    complex2model = search_pair()

    # Iterate over the list and model each complex
    for pair in complex2model:
        inference(pair)
    log_message(f'Finished all runs')

if __name__ == "__main__":
    main()
