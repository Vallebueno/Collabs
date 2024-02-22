#!/bin/bash

# Samtools Variant calling pipeline SLURM job script
# input mapped reads in bam format
# Author: Miguel Vallebueno  CC BY-NC-SA 4.0
# Date: 2024-02-22


##########example how to run:
# &sbatch MAVE_Samtools_Varcall.sh input_file.bam name_sample reference.fasta path2OUTDIR Quality_phred_cutoff



# === SLURM Config ===

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G
#SBATCH --qos=long        #-time +priority


jram=30
pcrs=5


# === Load Modules ===

module load build-env/2020
module load samtools/1.10-foss-2018b
module load bcftools/1.9-foss-2018b

# === INPUT Config ===

# where are the input files?

filei=$1
name=$2
REF=$3
WDIR=$4
QLT=$5

if [[ ! -d "${WDIR}" ]] ; then
 mkdir ${WDIR}
 fi
cd $WDIR


file=${WDIR}/${name}

echo "<<Caronte><Samtools_Hydra>>      <$INDIR>   <$name>  <$REF> <$OUT> <$QLT> <$WDIR> "


# === INSTRUCTIONS ===

####filter unmmaped reads and mark duplicates

samtools view -Sh -F 4 -q $QLT -b $filei > ${file}.onlymapped.bam

echo sort

samtools sort -n -@ $pcrs ${file}.onlymapped.bam -o ${file}.sortedn.bam

rm ${file}.onlymapped.bam

#https://github.com/samtools/samtools/issues/765  #fix mate pairs

echo fixmate

samtools fixmate -m ${file}.sortedn.bam ${file}.fixm.bam

rm ${file}.sortedn.bam

echo sort2

samtools sort -@ $pcrs ${file}.fixm.bam -o ${file}.sort2.bam

rm ${file}.fixm.bam

echo index

samtools index ${file}.sort2.bam

echo markdup

samtools markdup -@ $pcrs ${file}.sort2.bam ${file}.markdup.bam

rm ${file}.sort2.bam

echo index
samtools index ${file}.markdup.bam

echo mpileup
####Raw_variants call GVCF#####
bcftools mpileup -f ${REF} ${file}.markdup.bam -o ${file}.mpileup --threads $pcrs -g 1 -a "AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR"

echo call
bcftools call ${file}.mpileup -o ${file}.samtools.g.vcf --threads $pcrs -m -g 1


gzip ${file}.mpileup
gzip ${file}.samtools.g.vcf

#bcftools call ${file}.mpileup -o ${file}.samtools.v.vcf  --threads $pcrs -m -v
###############################

rm ${file}.markdup.bam
rm ${file}.markdup.bam.bai

###Generate bed
#zcat ${OUTDIR}/$file\_onlymapped$RNAM.sorted.markdup.g.vcf | grep -v "#" |gawk '{print$1"	"$2}' > ${OUTDIR}/$file\_onlymapped$RNAM.sorted.markdup.g.vcf.bed

################################

echo "<<Caronte><Samtools_Hydra>><DONE>"
