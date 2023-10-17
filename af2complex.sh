#!/bin/bash
#multimer parameters
DATA_DIR=/global/scratch/users/skyungyong/Software/alphafold/alphafold

### input targets
target_lst_file=$1
fea_dir=$2   # input feature pickle files of individual monomers under $inp_dir/$monomer
out_dir=$2   # model output files will be under $out_dir/$target

preset=deepmind
model=model_1_multimer_v3,model_2_multimer_v3,model_3_multimer_v3,model_4_multimer_v3,model_5_multimer_v3
model_preset=multimer

echo "Info: input feature directory is $fea_dir"
echo "Info: result output directory is $out_dir"
echo "Info: model preset is $model_preset"

# AF2Complex source code directory
af_dir=/global/scratch/users/skyungyong/Software/af2complex/src


python -u $af_dir/run_af2c_mod.py --target_lst_path=$target_lst_file \
  --data_dir=$DATA_DIR --output_dir=$out_dir --feature_dir=$fea_dir \
  --model_names=$model \
  --preset=$preset \
  --model_preset=$model_preset \
