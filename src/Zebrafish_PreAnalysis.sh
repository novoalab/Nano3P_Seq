##################################################
### PRE-ANALYSIS OF THE pA Selected Zebrafish Run ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############

###################################################
## PART 1 : BASECALLING AND DEMUX #################
###################################################

# BASECALLING USING GUPPY 3.6.1
guppy_3_6_1=/users/enovoa/lpryszcz/src/ont-guppy_3.6.1/bin/guppy_basecaller
input=/no_backup_isis/enovoa/data/ont/template_switching/cDNA852361_2_Zebrafish_pA_TGIRT_Flongle
output=/users/enovoa/boguzhan/NanoCapture/cDNA852361_2_Zebrafish_pA_TGIRT_Flongle
$guppy_3_6_1 --device cuda:0 -c dna_r9.4.1_450bps_hac.cfg  --fast5_out -ri $input -s $output

#Count how many reads are mapped
wc -l *fastq | awk '{x=$1/4; print x; print $2}'
# Run porechop for the unclassified


####################################################
## PART 2 : RUNNING TAILFINDR     ##################
####################################################

#Copy the Image for Singularity
rsync /no_backup_isis/enovoa/users/andelgado/master_of_pores_dev/singularity/biocorecrg-mopmod-1.0.img .

#Create fast5_dir.sh file
echo 'current_dir=`pwd`
f5file=$1
f5dir=${f5file/fast5/dir}
mkdir -p $f5dir 
mv $f5file $f5dir
csv=${f5file/fast5/csv} ' > fast5_dir.sh

#Make it executable
chmod +x fast5_dir.sh
#Run the script to split fast5 into folders
for i in *fast5; do ./fast5_dir.sh $i ; done

#Create the sh files to run the tailfindR
for INPUT in *.dir; do echo 'R --slave -e "library(tailfindr); find_tails(fast5_dir ='\'$INPUT\'',
save_dir = '\'$INPUT\'',
csv_filename = '\'$INPUT.csv\'' ,
num_cores = 10)" ' > $INPUT.sh ; done

#Make them exevutable
for i in *dir.sh; do chmod +x $i ; done


#Make singularity sh files
for i in *.dir.sh; do echo "module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load Singularity/3.2.1/ 
singularity exec -e biocorecrg-mopmod-1.0.img" './'$i > $i.sing.sh  ; done

#Submit the jobs
for i in *sing.sh; do qsub -cwd  $i ; done


#Merge the tail files
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.dir/*.csv > tails.csv



####################################################
## PART 3 : MAPPING     ############################
####################################################
ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_rRNA_sequins.fa


minimap2 -ax splice -k14 -uf --MD $ref cDNA8523612.fastq > cDNA8523612.bam


samtools view -hSb -F 3844 cDNA8523612.bam >  cDNA8523612.sam


samtools sort cDNA8523612.sam cDNA8523612.sorted && samtools index cDNA8523612.sorted.bam




for i in *.sorted.bam;do samtools view -F 4 $i | cut -f1 | sort | uniq | wc -l;echo $i; done

