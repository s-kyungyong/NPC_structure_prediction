# Predicting NPC heterodimer complexes

This section outlines our method to use [AF2Complex](https://github.com/FreshAirTonight/af2complex) v1.4.0 and [AlphaFold-Multimer](https://github.com/google-deepmind/alphafold) v2.3.0 to predict the heterodimer complex of the candidate NPC proteins. 

### Input sequences

We initially have 109 candidate sequences that can be found in [NPC.fasta](https://github.com/s-kyungyong/NPC_structure_prediction/blob/main/sequences.fasta) in this repository. Generate a folder and a fasta file for each accession.

```
grep ">" NPC.fasta | cut -d ">" -f 2 | while read accession; do
  mkdir -p "$accession" && awk -v acc="$accession" 'BEGIN {RS=">"} $0 ~ "^"acc"\\n" {print RS$0}' NPC.fasta > "$accession/$accession.fasta"
done

```

In the current working directory, we will have all the folders listed as below. 
```
ls -d AT*
AT1G07410  AT1G15290  AT1G55540  AT1G78300  AT2G30050  AT2G45640  AT3G14120  AT3G56900  AT4G16143  AT5G08450  AT5G48810
AT1G07970  AT1G15750  AT1G59660  AT1G79280  AT2G32080  AT3G01340  AT3G15690  AT3G57350  AT4G30840  AT5G10860  AT5G49880
AT1G09270  AT1G23170  AT1G63000  AT1G79920  AT2G32240  AT3G06720  AT3G15970  AT3G57410  AT4G31430  AT5G14060  AT5G51200
AT1G09620  AT1G24310  AT1G64350  AT1G80670  AT2G36130  AT3G06910  AT3G16310  AT3G57940  AT4G32285  AT5G16070  AT5G51280
AT1G09640  AT1G26630  AT1G67230  AT1G80680  AT2G37040  AT3G10650  AT3G18790  AT3G58110  AT4G32910  AT5G20200  AT5G52800
AT1G10390  AT1G27430  AT1G68790  AT2G05120  AT2G38770  AT3G11830  AT3G20370  AT4G02150  AT4G34660  AT5G38480  AT5G64930
AT1G13120  AT1G33410  AT1G73240  AT2G19520  AT2G39630  AT3G11910  AT3G25980  AT4G11660  AT4G37130  AT5G40480  AT5G65770
AT1G13160  AT1G35160  AT1G73660  AT2G19860  AT2G39810  AT3G12080  AT3G53110  AT4G11790  AT4G38760  AT5G42080  AT5G67240
AT1G13220  AT1G48610  AT1G75340  AT2G25970  AT2G41620  AT3G13200  AT3G53260  AT4G15880  AT5G05680  AT5G42950  ATCG00180
AT1G14850  AT1G52380  AT1G77180  AT2G27140  AT2G45000  AT3G13870  AT3G55610  AT4G15900  AT5G05970  AT5G46070
```

There should be a sequence file in each folder. 
```
cat AT1G07410/AT1G07410.fasta

>AT1G07410
MANRIDHEYDYLFKIVLIGDSGVGKSNILSRFTRNEFCLESKSTIGVEFATRTLQVEGKTVKAQIWDTAGQERYRAITSAYYRGAVGALLVYDITKRQTFENVLRWLRELRDHADSNIVIMMAGNKSDLNHLRSVADEDGRSLAEKEGLSFLETSALEATNIEKAFQTILSEIYHIISKKALAAQEAAGNLPGQGTAINISDSSATNRKGCCST
```


### Multiple sequence aligments

For each sequence, we collect homologous sequences and construct multiple sequence alignments (MSAs), following the standard AlphaFold pipeline. To speed up this step, we use customized scripts derived from AlphaFold. These scripts will look for any folders that start with 'AT' in the current working directory.

```
compute_sto.py AT
compute_a3m.py AT
```

After these scripts are run, the input MSAs are generated in each accession folder. 

```
ls AT1G07410/AT1G07410/msas/
bfd_uniref_hits.a3m  hhblits.out  hmm_output.sto  mgnify_hits.sto  uniprot_hits.sto  uniref90_hits.sto
```

### Structure inference

All heterodimer complexes that are equal to or smaller than 2,500 amino acids were predicted with AF2Complex with AlphaFold-Multimer models (v3). This script also relies on generate_feature.sh and af2complex.sh. We then collect the prediction confidence scores for the five PDB files for each prediction.  

```
mkdir structure
python inference.py
python collect_scores.py
```

[NPC_prediction_scores.txt](https://github.com/s-kyungyong/NPC_structure_prediction/blob/main/NPC_prediction_scores.txt) contains the best scores from the five models for each heterodimer complex prediction. 
### Data visualization

The putative protein-protein interaction network and the prediction confidence scores were initially visualized with R and modified in Illustrator. 

```
generate_graph.r heatmap.r
```
