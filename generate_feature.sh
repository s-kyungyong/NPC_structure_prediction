#!/bin/bash

#parameters
DATA_DIR=$HOME/scratch/afold/data

#databases
UNIREF=/global/scratch/users/skyungyong/Software/alphafold/Database/uniref90/uniref90.fungalDB.added.fasta
MGNIFY=/global/scratch/users/skyungyong/Software/alphafold/Database/mgnify/mgy_clusters.fa
UNIPROT=/global/scratch/users/skyungyong/Software/alphafold/Database/uniprot/uniprot.fasta
PDB=/global/scratch/users/skyungyong/Software/AF2-DB_2023-03-31/pdb_seqres.txt
BFD=/global/scratch/users/skyungyong/Software/alphafold/Database/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
UNICLUST=/global/scratch/users/skyungyong/Software/alphafold/Database/uniclust30/uniclust30_2018_08/uniclust30_2018_08
MMCIF=/global/scratch/users/skyungyong/Software/AF2-DB_2023-03-31/pdb_mmcif/
SMALLBFD=/global/scratch/users/skyungyong/Software/alphafold/Database/bfd_small/bfd-first_non_consensus_sequences.fasta

#software
export HHLIB=/global/scratch/users/skyungyong/Software/hhsuite-3.3.0/
export HMMER=/global/scratch/users/skyungyong/Software/anaconda3/envs/alphafold/bin/
export KALIGN=/global/scratch/users/skyungyong/Software/anaconda3/envs/alphafold/bin/

#af2complex
af_dir=/global/scratch/users/skyungyong/Software/af2complex/src

if [ $# -eq 0 ]
  then
    echo "Usage: $0 <seq_file>"
    exit 1
fi

#input parameters
fasta_path=$1
out_dir=$2
db_preset='full_dbs'
feature_mode='multimer'
max_template_date=2030-01-01


echo "Info: sequence file is $fasta_path"
echo "Info: out_dir is $out_dir"
echo "Info: db_preset is $db_preset"
echo "Info: feature mode is $feature_mode"
echo "Info: max_template_date is $max_template_date"


##########################################################################################


if [ "$feature_mode" = "multimer" ] || [ "$feature_mode" = "monomer+fullpdb" ]; then
  python $af_dir/run_af2c_fea.py --fasta_paths=$fasta_path --db_preset=$db_preset \
    --data_dir=$DATA_DIR --output_dir=$out_dir      \
    --uniprot_database_path=$UNIPROT \
    --uniclust30_database_path=$UNICLUST \
    --uniref90_database_path=$UNIREF \
    --mgnify_database_path=$MGNIFY   \
    --pdb_seqres_database_path=$PDB \
    --bfd_database_path=$BFD \
    --template_mmcif_dir=$MMCIF/mmcif_files/ \
    --max_template_date=$max_template_date   \
    --obsolete_pdbs_path=$MMCIF/obsolete.dat \
    --hhblits_binary_path=$HHLIB/bin/hhblits   \
    --hhsearch_binary_path=$HHLIB/bin/hhsearch \
    --jackhmmer_binary_path=$HMMER/jackhmmer \
    --hmmsearch_binary_path=$HMMER/hmmsearch \
    --hmmbuild_binary_path=$HMMER/hmmbuild \
    --kalign_binary_path=$KALIGN/kalign \
    --feature_mode=$feature_mode \
    --use_precomputed_msas=True
fi
