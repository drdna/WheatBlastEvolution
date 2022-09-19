#!/bin/bash

#SBATCH --time 12:00:00
#SBATCH --job-name=bowtie2-GATK
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=<partition>  
#SBATCH --mem=100GB
#SBATCH --mail-type ALL
#SBATCH	-A <account-info>
#SBATCH --mail-type ALL
#SBATCH --mail-user >user@email.com>

echo "SLURM_NODELIST: "$SLURM_NODELIST
echo "PWD :" $PWD
echo "index: $1"
echo "path to reads: $2/$3"


#### RUN BOWTIE2

# run bowtie2-build

  index=$1
  readsdir=$2
	
	indexdir=${index}_index

#reads 

	reads=$3

#bowtie2-build $ref $indexdir/${ref/.fasta/}

# create output directory

	outdir=${index}_${reads/*\//}_ALIGN
	mkdir $outdir

# run bowtie2 alignment

if [ $4 == 'SE' ]
  then
	bowtie2 --threads 16 --very-sensitive-local --phred33 --no-unal -x $indexdir/$index \
	 -U $readsdir/$reads*\.f*q* 2>$outdir/alignment_summary.txt |  samtools view -bS - > $outdir/accepted_hits.bam
else
	bowtie2 --threads 16 --very-sensitive-local --phred33 --no-unal -x $indexdir/$index \
 	-1 $readsdir/${reads}_*1.f*q* -2 $readsdir/${reads}_*2.f*q* 2>$outdir/alignment_summary.txt | samtools view -bS - > $outdir/accepted_hits.bam
fi

# deleted: --no-unal

#bowtie2 --threads 16 --very-sensitive-local --phred33 -x $indexdir/$index \
# -1 $readsdir/${reads}_*1.f*q* -2 $readsdir/${reads}_*2.f*q* | samtools view -bS - > $outdir/accepted_hits.bam

# sort bamfile and remove original

	samtools sort $outdir/accepted_hits.bam -o $outdir/accepted_hits_sorted.bam
	rm $outdir/accepted_hits.bam

# add read group info to bam header and remove source file

	o=${reads/*\//}
	conda activate gatk
	java -jar picard.jar AddOrReplaceReadGroups I=$outdir/accepted_hits_sorted.bam \
	O=$outdir/accepted_hits_sortedRG.bam RGID=$o RGSM=$o RGLB=$o RGPI=50 RGPL=illumina RGPU=unit1
	
	rm $outdir/accepted_hits_sorted.bam

# index the bamfile

  samtools index $outdir/accepted_hits_sortedRG.bam
	


##### RUN GATK

	# specify index information

	fasta=${indexdir}/$index.fasta


# give reads prefix (base ID)

	readsprefix=$3


# specify input/output directory

	indir=$outdir


# index reference fasta file

	samtools faidx $fasta


# Create dict and fai for reference

# remove existing dictionary

	rm ${fasta/fasta/dict}

#create dictionary

	java -jar picard.jar CreateSequenceDictionary R=$fasta O=${fasta/fasta/dict}

# Call haplotype

	gatk --java-options "-Xmx35g" HaplotypeCaller \
        	--native-pair-hmm-threads 16 \
	      	-R $fasta \
        	-ploidy 1 \
	        -I $indir/accepted_hits_sortedRG.bam \
        	--emit-ref-confidence GVCF \
	        -O $indir/${readsprefix}.vcf \

# Determine genotypes

	gatk --java-options "-Xmx35g" GenotypeGVCFs \
	 -R $fasta \
	 -V $indir/$readsprefix.vcf \
	 --output $indir/${readsprefix}_genotyped.vcf


# Filter SNPs only

	gatk SelectVariants \
	 -R $fasta \
	 -select-type SNP \
	 -V $indir/${readsprefix}_genotyped.vcf \
	 --output $indir/${readsprefix}_genotyped-snps.vcf


	vcftools --vcf $indir/${readsprefix}_genotyped.snp-only.filters.vcf \
	        --remove-filtered-all --recode \
	        --out $indir/${readsprefix}_genotyped.snp-only.filters.subset.vcf
