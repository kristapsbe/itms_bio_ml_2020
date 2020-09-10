######### MetaVW #########

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


######## Deepmicrobes #########

#Download the genomes in the complete bacterial repertoire of the human gut microbiota from this FTP site.
#Assign a category label for each species.
#Read simulation with ART simulator for each genome.
#Trim reads to variable lengths with the custom script random_trim.py.
#Shuffle all the reads and convert them to TFRecord. 

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


# GeNet

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

# fastdna

awk '/^>/{printf(toupper("%s|%d\n"),$0,++counter);next}{print $0}' train_small-db.fasta > train_small-db.fasta.fixed
awk '/^>/{printf(toupper("%s|%d\n"),$0,++counter);next}{print $0}' train_small-db.fasta > train_small-db.fasta.fixed

awk '{print toupper($0)}' train_small-db.fasta > train_small-db.fasta.fixed

./fastdna supervised -input train_small-db.fasta -labels train_small-db.species-level.taxid -output model -minn 12

# kraken2 

PATH=$PATH:$(pwd)/

date && kraken2-build --download-taxonomy --db itms --use-ftp && date 
date && kraken2-build --add-to-library /home/kristaps/Projs/itms/tools/MetaVW/data/train-dataset/train_small-db.fasta --db itms --no-masking && kraken2-build --build --db itms --threads 12 && date