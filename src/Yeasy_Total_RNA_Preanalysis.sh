##################################################
### PRE-ANALYSIS OF THE Yeast TEMPLATE SWT RUN ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############

###################################################
## PART 1 : BASECALLING AND DEMUX #################
###################################################

# BASECALLING USING GUPPY 3.6.1
guppy_3_6_1=/users/enovoa/lpryszcz/src/ont-guppy_3.6.1/bin/guppy_basecaller
input=/no_backup_isis/enovoa/data/ont/template_switching/yeast_RNA_template_switching/
output=/users/enovoa/boguzhan/NanoCapture/cDNA231532_Yeast_Total
$guppy_3_6_1 --device cuda:0 -c dna_r9.4.1_450bps_hac.cfg  --fast5_out -ri $input -s $output


#Count how many reads are mapped
wc -l *fastq | awk '{x=$1/4; print x; print $2}'



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




ref=/users/enovoa/boguzhan/references/s.cerevisiae/saccer_all_rRNA.fa

minimap2 -ax splice -k14 -uf --MD $ref cDNA231532.fastq > cDNA231532_SacCer_Genome_rRNA.bam


samtools view -hSb -F 3844 cDNA231532_SacCer_Genome_rRNA.bam>  cDNA231532_SacCer_Genome_rRNA.sam


samtools sort cDNA231532_SacCer_Genome_rRNA.sam cDNA231532_SacCer_Genome_rRNA.sorted && samtools index cDNA231532_SacCer_Genome_rRNA.sorted.bam

for i in *.sorted.bam;do samtools view -F 4 $i | cut -f1 | sort | uniq | wc -l;echo $i; done



