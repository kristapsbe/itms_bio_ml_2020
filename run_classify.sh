#metavw

cd MetaVW/src/2-build-models/src

# clear out if not first run
rm -rf ../output/train_small-db/
rm -rf ../output/TMP

bash 01.main.sh

# the thing doesn't support paired end reads - merge both files
bbmerge.sh in1=control_1_1.fq.gz in2=control_1_2.fq.gz out=control_1.fq.gz


# TODO: make a custom prediction script


/home/kristaps/Projs/itms/samples/control_1.fq.gz








# fastdna 

for f in ../../samples/merged/*.fq.gz
do
    date
    ./fastdna predict model.bin $f > $f.fastdna.classif
    date
done

trešdiena, 2020. gada 19. augusts, 20:49:37 EEST
trešdiena, 2020. gada 19. augusts, 20:50:06 EEST

trešdiena, 2020. gada 19. augusts, 20:50:06 EEST
trešdiena, 2020. gada 19. augusts, 20:50:50 EEST

trešdiena, 2020. gada 19. augusts, 20:50:50 EEST
trešdiena, 2020. gada 19. augusts, 20:51:05 EEST

trešdiena, 2020. gada 19. augusts, 20:51:05 EEST
trešdiena, 2020. gada 19. augusts, 20:51:30 EEST

trešdiena, 2020. gada 19. augusts, 20:51:30 EEST
trešdiena, 2020. gada 19. augusts, 20:52:04 EEST

# kraken2
for f in ../../samples/merged/*.fq.gz
do
    date
    kraken2 --db itms $f --threads 12 > $f.kraken2.classif
    date
done

for f in ../../samples/merged/*.fq.gz
do
    fname=${f/..\/..\/samples\/merged\//}
    fname=${fname/.fq.gz/}
    date
    kraken2 --paired --db itms --threads 12 ../../samples/${fname}_1.fq.gz ../../samples/${fname}_2.fq.gz > $f.kraken2_paired.classif
    date
done

trešdiena, 2020. gada 19. augusts, 17:53:58 EEST
trešdiena, 2020. gada 19. augusts, 17:54:29 EEST

trešdiena, 2020. gada 19. augusts, 17:54:29 EEST
trešdiena, 2020. gada 19. augusts, 17:55:15 EEST

trešdiena, 2020. gada 19. augusts, 17:55:16 EEST
trešdiena, 2020. gada 19. augusts, 17:55:32 EEST

trešdiena, 2020. gada 19. augusts, 17:55:32 EEST
trešdiena, 2020. gada 19. augusts, 17:55:58 EEST

trešdiena, 2020. gada 19. augusts, 17:55:59 EEST
trešdiena, 2020. gada 19. augusts, 17:56:35 EEST


# metavw
# PATH=$PATH:$(pwd)/ 
# a bunch of stuff to make this work
for f in ../../samples/merged/*.fq.gz
do
    date
    fasta2vw -i $f -k 12 | vw -t -i src/2-build-models/output/train_small-db/vw-model_batch-10.model -p $f.preds.vw
    # convert vw class to taxid
    Rscript src/3-make-predictions/src/vw-class-to-taxid.R $f.preds.vw src/2-build-models/output/train_small-db/vw-dico.txt $f.preds.vw
    date
done

trešdiena, 2020. gada 19. augusts, 18:47:09 EEST
trešdiena, 2020. gada 19. augusts, 19:10:42 EEST

trešdiena, 2020. gada 19. augusts, 19:10:42 EEST
trešdiena, 2020. gada 19. augusts, 19:46:12 EEST

trešdiena, 2020. gada 19. augusts, 19:46:12 EEST
trešdiena, 2020. gada 19. augusts, 19:58:35 EEST

trešdiena, 2020. gada 19. augusts, 19:58:35 EEST
trešdiena, 2020. gada 19. augusts, 20:18:52 EEST

trešdiena, 2020. gada 19. augusts, 20:18:52 EEST
trešdiena, 2020. gada 19. augusts, 20:46:44 EEST

# DeepMicrobes
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



for i in 0 1 2 3 4 5 6 7 8 9
do
    sample_name="V300024282_L02_97"
    fname=../../samples/${i}_${sample_name}
    tfrec_predict_kmer.sh -f ${fname}_1.fq -r ${fname}_2.fq -t fastq -v $(pwd)/tokens_merged_11mers.txt -o ${i}_${sample_name} -s 4000000 -k 11 
    predict_DeepMicrobes.sh -i ${i}_${sample_name}.tfrec -b 2048 -l species -p 8 -m $(pwd)/train_data/test_db_split/single/outfiles/model/ -o ${i}_${sample_name} 
    rm *.tfrec
done

ls ../../samples/*.fq

f=../../samples/0_V300024282_L02_104_1.fq
fname=${f/..\/..\/samples\/merged\//}
fname=${fname/.fq.gz/}
cp ../../samples/${fname}* .
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

C1
pirmdiena, 2020. gada 31. augusts, 00:31:30 EEST
pirmdiena, 2020. gada 31. augusts, 03:20:21 EEST
pirmdiena, 2020. gada 31. augusts, 08:56:06 EEST
pirmdiena, 2020. gada 31. augusts, 14:55:35 EEST
C2
otrdiena, 2020. gada  1. septembris, 09:51:47 EEST
otrdiena, 2020. gada  1. septembris, 13:31:41 EEST
otrdiena, 2020. gada  1. septembris, 21:10:44 EEST
1
trešdiena, 2020. gada  2. septembris, 15:16:10 EEST
trešdiena, 2020. gada  2. septembris, 17:18:39 EEST
trešdiena, 2020. gada  2. septembris, 21:09:22 EEST
