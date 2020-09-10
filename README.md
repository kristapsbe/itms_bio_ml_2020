# Tools

## Kraken 2

### Setup

``` bash
PATH=$PATH:$(pwd)/

kraken2-build --download-taxonomy --db itms --use-ftp
kraken2-build --add-to-library /home/kristaps/Projs/itms/tools/MetaVW/data/train-dataset/train_small-db.fasta --db itms --no-masking 
kraken2-build --build --db itms --threads 12
```

### Usage

``` bash
for f in ../../samples/merged/*.fq.gz
do
    fname=${f/..\/..\/samples\/merged\//}
    fname=${fname/.fq.gz/}
    date
    kraken2 --paired --db itms --threads 12 ../../samples/${fname}_1.fq.gz ../../samples/${fname}_2.fq.gz > $f.kraken2_paired.classif
    date
done
```

## MetaVW

### Setup

``` bash
wget http://projects.cbio.mines-paristech.fr/largescalemetagenomics/large-scale-metagenomics-1.0.tar.gz

git clone https://github.com/lucren/MetaVW/

cd MetaVW/tools

bash INSTALL.sh

chmod +x drawfrag
chmod +x fasta2vw

cd ext/gdl-1.2/GDL/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/
cd ../../..

wget http://gnu.mirror.vexxhost.com/gsl/gsl-1.6.tar.gz
tar -zxvf gsl-1.6.tar.gz
cd gsl-1.6
./configure
make
sudo make install
cd .libs
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/
cd ../..

wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.86.tar.gz
tar -zxvf BBMap_38.86.tar.gz
cd bbmap/
PATH=$PATH:$(pwd)/
cd ..

wget https://github.com/lucren/vowpal_wabbit/archive/7.7.tar.gz
tar -zxvf 7.7.tar.gz
cd vowpal_wabbit-7.7
make 
cd vowpalwabbit
PATH=$PATH:$(pwd)/

cd ../../../..

cd src/2-build-models/src/
bash 01.main.sh
```

### Usage

``` bash
for f in ../../samples/merged/*.fq.gz
do
    date
    fasta2vw -i $f -k 12 | vw -t -i src/2-build-models/output/train_small-db/vw-model_batch-10.model -p $f.preds.vw
    # convert vw class to taxid
    Rscript src/3-make-predictions/src/vw-class-to-taxid.R $f.preds.vw src/2-build-models/output/train_small-db/vw-dico.txt $f.preds.vw
    date
done
```

## fastDNA

### Setup

``` bash
awk '/^>/{printf(toupper("%s|%d\n"),$0,++counter);next}{print $0}' train_small-db.fasta > train_small-db.fasta.fixed
awk '/^>/{printf(toupper("%s|%d\n"),$0,++counter);next}{print $0}' train_small-db.fasta > train_small-db.fasta.fixed

awk '{print toupper($0)}' train_small-db.fasta > train_small-db.fasta.fixed

./fastdna supervised -input train_small-db.fasta -labels train_small-db.species-level.taxid -output model -minn 12
```

### Usage

``` bash
for f in ../../samples/merged/*.fq.gz
do
    date
    ./fastdna predict model.bin $f > $f.fastdna.classif
    date
done
```

## GeNet

### Setup

``` bash
for f in label*.fasta; do echo "${f/label_/};$(grep -o -a -m 1 -h -r -oP '\|\K[^|]+' $f | head -1)"; done


cd data
tar -xvzf nodes.dmp.tar.gz
rm nodes.dmp.tar.gz

# pārkopē deepmicrobes kopu un jauno csv - datiem jāiet data/saved_models mapē
# NB jātiek skaidrībā ar dep install'u

# remember - ir vajadzīgi drusku modi, lai tas viss actually strādātu

# some sort of weird cudnn issue - disabling for now
# NB - remeber to not let it install cudnn 7.1.2 in the conda env
conda uninstall cudnn

export CUDA_VISIBLE_DEVICES=0
python genet_train.py --dir_path=../data #--read_length=150
```

## DeepMicrobes

### Setup

``` bash
PATH=$PATH:$(pwd)/
cd bin
rm parallels # doesn't work for some reason - install via apt instead
PATH=$PATH:$(pwd)/
cd ../scripts
PATH=$PATH:$(pwd)/
cd ../pipelines
PATH=$PATH:$(pwd)/
cd ..

conda activate DeepMicrobes

mkdir train_data
# copy the make scripts and the data
cd train_data

python3.8 generate_deepmicrobes.py
cd test_db_split
mkdir outfiles

# the script's a bit unintuitive to use
python3 fna_label.py -m label_file.txt -o outfiles

wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
tar -zxvf artbinmountrainier2016.06.05linux64.tgz

folds=1 # really low to try and fit in the gpu

mkdir -p single/outfiles
#mkdir -p paired/outfiles

#for f in outfiles/*.fasta
#do
#    # simulating both single and paired-end data
#    art_bin_MountRainier/art_illumina -ss HS25 -i $f -l 150 -f $folds -o single/$f
#    art_bin_MountRainier/art_illumina -ss HS25 -i $f -p -l 150 -f $folds -m 200 -s 10 -o paired/$f
#done

for f in outfiles/*.fasta
do
    # simulating both single and paired-end data
    art_bin_MountRainier/art_illumina -ss HS25 -i $f -l 150 -f $folds -o single/$f
done

cd single/outfiles
#cd paired/outfiles

for f in *.fq
do
    python3 ../../../../scripts/random_trim.py -i $f -o trimmed_$f -f fastq -l 150 -min 0 -max 75
done 

cat trimmed_*.fq > train.fa
#wget https://github.com/MicrobeLab/DeepMicrobes-data/raw/master/vocabulary/tokens_merged_12mers.txt.gz
#gunzip tokens_merged_12mers.txt.gz

wget https://github.com/MicrobeLab/DeepMicrobes-data/raw/master/vocabulary/tokens_merged_11mers.txt.gz
gunzip tokens_merged_11mers.txt.gz

# NB remeber that this will not work with tf 2
#bash tfrec_train_kmer.sh -i train.fa -v $(pwd)/tokens_merged_12mers.txt -o train.tfrec -s 20480000 -k 12

conda activate DeepMicrobes
bash tfrec_train_kmer.sh -i train.fa -v $(pwd)/tokens_merged_11mers.txt -o train.tfrec -s 20480000 -k 11 && date


mkdir model
#python3 ../../../../DeepMicrobes.py --input_tfrec=train.tfrec --model_name=attention --model_dir=$(pwd)/model

export CUDA_VISIBLE_DEVICES=0
python3 ../../../../DeepMicrobes.py --input_tfrec=train.tfrec --model_name=attention --model_dir=$(pwd)/model --kmer 11 --vocab_size 2097154


tfrec_predict_kmer.sh -f control_1_1.fq.gz -r control_1_2.fq.gz -t fastq -v $(pwd)/train_data/test_db_split/single/outfiles/tokens_merged_11mers.txt -o control_1 -s 4000000 -k 11


bash tfrec_train_kmer.sh -i train.fa -v $(pwd)/tokens_merged_11mers.txt -o train.tfrec -s 20480000 -k 11 && date 
mkdir model && export CUDA_VISIBLE_DEVICES=0 && python3 ../../../../DeepMicrobes.py --input_tfrec=train.tfrec --model_name=attention --model_dir=$(pwd)/model --kmer 11 --vocab_size 2097154 && date

########
#>>><<<#
########

cd ../../paired/outfiles

for f in *.fq
do
    python3 ../../../../scripts/random_trim.py -i $f -o trimmed_$f -f fastq -l 150 -min 0 -max 75
done 
```

### Usage

``` bash
wget https://github.com/MicrobeLab/DeepMicrobes-data/raw/master/vocabulary/tokens_merged_11mers.txt.gz
gunzip tokens_merged_11mers.txt.gz

for f in ../../samples/merged/*.fq.gz
do
    # just use the paired stuff
    fname=${f/..\/..\/samples\/merged\//}
    fname=${fname/.fq.gz/}
    cp ../../samples/${fname}* .
    gunzip ${fname}_1.fq.gz
    gunzip ${fname}_2.fq.gz
    echo ""
    date
    echo ""
    tfrec_predict_kmer.sh -f ${fname}_1.fq -r ${fname}_2.fq -t fastq -v $(pwd)/tokens_merged_11mers.txt -o ${fname} -s 4000000 -k 11 
    echo ""
    date
    echo ""
    predict_DeepMicrobes.sh -i ${fname}.tfrec -b 2048 -l species -p 8 -m $(pwd)/train_data/test_db_split/single/outfiles/model/ -o $fname 
    echo ""
    date
    echo ""
    rm ${fname}_1.fq
    rm ${fname}_2.fq
    rm ${fname}.tfrec
done
```

``` bash
for i in 0 1 2 3 4 5 6 7 8 9
do
    sample_name="V300024282_L02_97"
    fname=../../samples/${i}_${sample_name}
    tfrec_predict_kmer.sh -f ${fname}_1.fq -r ${fname}_2.fq -t fastq -v $(pwd)/tokens_merged_11mers.txt -o ${i}_${sample_name} -s 4000000 -k 11 
    predict_DeepMicrobes.sh -i ${i}_${sample_name}.tfrec -b 2048 -l species -p 8 -m $(pwd)/train_data/test_db_split/single/outfiles/model/ -o ${i}_${sample_name} 
    rm *.tfrec
done
```

# Notebooks and Utilities

* coverage_and_precision.ipynb
* get_taxo_levels.ipynb
* split_up_samples.ipynb
* generate_deepmicrobes.py

# Article