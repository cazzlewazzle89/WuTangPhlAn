# Enter the Wu-TangPhlAn
The goal is to turn this into a (w)rapper script for **Ph**y**l**ogenetic **An**alysis of metagenomic data.  
At the moment it is a work-in-progress tutorial mainly applicable for use on the hcux400 HPC at Teagasc but the general code should work on other machines.  

### Table of Contents
[Setup](https://github.com/cazzlewazzle89/WuTangPhlAn#setup)  
[Quality Control](https://github.com/cazzlewazzle89/WuTangPhlAn#quality-control)  
[Microbiome Profiling](https://github.com/cazzlewazzle89/WuTangPhlAn#microbiome-profiling)  
[Metagenome Assembly](https://github.com/cazzlewazzle89/WuTangPhlAn#metagenome-assembly-using-metaspades)  
[Strain-level Analysis](https://github.com/cazzlewazzle89/WuTangPhlAn#strain-level-analysis)  
[Recovery of Metagenome-Assembled Genomes (MAGs)](https://github.com/cazzlewazzle89/WuTangPhlAn#recovery-of-metagenome-assembled-genomes-mags)

## Setup  
Your working directory should contain
* a `RawFastQ/` directory with demultiplexed, lane-merged, paired-end metagenomic reads
  * named in the format `samplename_R[1-2].fastq.gz`
* a text file named `samplenames.txt` which lists each unique `samplename`, one per line.

For example, in an experiment examining two samples, a mock community, and a negative control - your RawFastQ/ directory would contain the following files

```
sampleA_R1.fastq.gz  
sampleA_R2.fastq.gz  
sampleB_R1.fastq.gz  
sampleB_R2.fastq.gz  
mock_R1.fastq.gz  
mock_R2.fastq.gz  
neg_R1.fastq.gz  
neg_R2.fastq.gz
```

and your samplenames.txt file would look like this

```
sampleA  
sampleB  
mock  
neg
```

## Quality control
The two main steps in quality control of metagenomic sequencing data are:  
* adapter removal and quality trimming
* removal of contaminant reads (usually host DNA is we are studying human, animal, or food microbiomes)

One option is to run these steps separately (the first code block below) or use a wrapper script that combines these steps into one command (the second code block)

If running these separately, I will usually perform adapter removal and quality trimming using [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)  followed by contamination removal using [Bowtie2](https://github.com/BenLangmead/bowtie2).  
TrimGalore is a wrapper tool around [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
Bowtie2 is a very popular aligner which we can use to identify contaminant DNA by aligning paried-end metagenomic reads against a database (eg. the human or bovine genome).  

The steps in this script are as follows:
* create output directories to store quality-trimming reports, read quality summaries, and quality-trimmed reads
* run TrimGalore, tidy and rename output files, and delete intermediate files
* create an output directory for quality-trimmed, non-contaminant reads
* Run Bowtie2 and [Samtools](http://www.htslib.org/) to remove contamination by
  * aligning quality-trimmed reads to the human genome (you may need to use a more appropriate genome/database depending on your sample source)
  * converting the alignment output to binary output to save space
  * removing reads which aligned to the contaminant database
  * converting the alignment back into fastq format
 The files are stored in the MicrobialFastQ/ directory and are the metagenomic reads which you will use for downstream analysis.  
 Note: the steps in the example below all specify that 24 CPUs are used for every task. You will need to specify this in your SLURM batch script when submitting the job using `#SBATCH --cpus-per-task=8` or adjust the number of however many CPUs you would like to use.  
 They also use default parameters. You may need to change these based on your data. 


```bash
mkdir TrimmingReports/
mkdir FastQC/
mkdir TrimmedFastQ/

for i in $(cat samplenames.txt)
do
    module load cutadapt/2.6
    module load fastqc/0.11.8
    ~/TrimGalore-0.6.0/trim_galore --paired RawFastQ/"$i"_R1.fastq.gz RawFastQ/"$i"_R2.fastq.gz --fastqc -j 8 -o "$i"_trimout
    module load cutadapt/2.6
    module unload fastqc/0.11.8

    mv "$i"_trimout/"$i"_R1_val_1.fq.gz TrimmedFastQ/"$i"_trimmed_R1.fastq.gz
    mv "$i"_trimout/"$i"_R2_val_2.fq.gz TrimmedFastQ/"$i"_trimmed_R2.fastq.gz
    mv "$i"_trimout/*_trimming_report.txt TrimmingReports/
    mv "$i"_trimout/*.html FastQC/

    rm -r "$i"_trimout
done

mkdir MicrobialFastQ/

for i in $(cat samplenames.txt)
do
    module load bowtie2/2.3.4
    bowtie2 -x /data/databases/hostremoval/Homo_sapiens/Bowtie2/Homo_sapiens -1 TrimmedFastQ/"$i"_trimmed_R1.fastq.gz -2 TrimmedFastQ/"$i"_trimmed_R2.fastq.gz -S "$i".sam -p 8
    module unload bowtie2/2.3.4

    module load samtools/1.10
    samtools view -bS "$i".sam -@ 8 > "$i".bam
    rm "$i".sam
    samtools view -b -f 12 -F 256 "$i".bam -@ 8 > "$i"_nonhost.bam
    rm "$i".bam
    samtools sort -n "$i"_nonhost.bam -T "$i"_nonhost_temp -o "$i"_nonhost_sorted.bam -@ 8
    rm "$i"_nonhost.bam
    module unload samtools/1.10

    module load bedtools/2.27.1
    bedtools bamtofastq -i "$i"_nonhost_sorted.bam -fq MicrobialFastQ/"$i"_microbial_R1.fastq -fq2 MicrobialFastQ/"$i"_microbial_R2.fastq
    rm "$i"_nonhost_sorted.bam
    module unload bedtools/2.27.1

    gzip MicrobialFastQ/"$i"_microbial_R1.fastq MicrobialFastQ/"$i"_microbial_R2.fastq
done
```

If using a wrapper script for both steps I will usually use [KneadData](https://github.com/biobakery/kneaddata).  
This wrapper script uses [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) for adapter removal and quality trimming, folowed by contaminant read removal using either Bowtie2 (default) or [BMTagger](https://www.westgrid.ca/support/software/bmtagger).

```bash
mkdir Knead_Outputs/
mkdir MicrobialFastQ/
mkdir Kneaddata_Logs/

module load kneaddata/0.6.1
source activate kneaddata_0.6.1

for i in $(cat samplenames.txt)
do
    kneaddata -i "$i"_R1.fastq.gz -i "$i"_R2.fastq.gz -db /data/databases/hostremoval/Homo_sapiens/Bowtie2/ --output Knead_Outputs/"$i" --threads 8 --output-prefix "$i" --remove-intermediate-output
    
    mv Knead_Outputs/"$i"/"$i"_paired_1.fastq MicrobialFastQ/"$i"_microbial_R1.fastq
    mv Knead_Outputs/"$i"/"$i"_paired_2.fastq MicrobialFastQ/"$i"_microbial_R2.fastq 
    mv Knead_Outputs/"$i"/"$i".log Kneaddata_Logs/"$i"_log.txt
    rm -r Knead_Outputs/"$i"/
done

conda deactivate
module unload kneaddata_0.6.1

gzip MicrobialFastQ/*.fastq
```

## Microbiome Profiling
This script perfoms taxonomic and functional profiling using [HUMANn3](https://github.com/biobakery/humann) and will
* create directories for each output
* create interleaved fastq files using [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
* run the HUMANn3 pipeline
* remove the interleaved fastq files
* convert the MetaPhlAn3 alignment output to BAM format
* tidy and rename output directories
* join, regroup, and rename genefamily-level functional outputs based on different protein/enzyme family categorisation systems
  * In the example below I use the default uniref90 scheme but also regroup based on EC number and Gene Ontology terms. There are more systems supported and the naming/grouping files are located in `data/databases/Humann3/utility_mapping/`
* normalise functional outputs to "copies per million reads" to account for inter-sample vairation in sequencing depth
* split functional outputs into two files
  * `stratified` outputs describe the predicted functional capacity of the community separated by species
  * `unstratified`outputs describe the function capacity of the community as a whole

The final output files are listed here and can be imported into your analysis/visualision software of choice to get started on the fun stuff
* `Metaphlan3outputs/metaphlan.tsv`
* `Humann3outputs/genefamilies_uniref90_name_cpm_stratified.tsv`
* `Humann3outputs/genefamilies_uniref90_name_cpm_unstratified.tsv`
* `Humann3outputs/genefamilies_level4ec_name_cpm_stratified.tsv`
* `Humann3outputs/genefamilies_level4ec_name_cpm_unstratified.tsv`
* `Humann3outputs/genefamilies_go_name_cpm_stratified.tsv`
* `Humann3outputs/genefamilies_go_name_cpm_unstratified.tsv`
* `Humann3outputs/pathabundance_cpm_stratified.tsv`
* `Humann3outputs/pathabundance_cpm_unstratified.tsv`
```
mkdir Metaphlan3outputs/
mkdir Metaphlan3outputs/Sam/
mkdir Metaphlan3outputs/MPAoutput/
mkdir Humann3outputs/
mkdir Humann3outputs/Genefamilies/
mkdir Humann3outputs/Pathabundance/
mkdir Humann3outputs/Pathcoverage/
mkdir Humann3outputs/Logs/

for i in $(cat samplenames.txt)
do
    module load bbmap/38.22
    reformat.sh in=MicrobialFastQ/"$i"_microbial_R1.fastq.gz in2=MicrobialFastQ/"$i"_microbial_R2.fastq.gz out=MicrobialFastQ/"$i"_interleaved.fastq.gz
    module unload bbmap/38.22

    module load humann3/3.0
    module load diamond/0.9.24
    humann -i MicrobialFastQ/"$i"_interleaved.fastq.gz -o "$i"_humann3out --nucleotide-database /data/databases/Humann3/chocophlan --protein-database /data/databases/Humann3/uniref --search-mode uniref90 --metaphlan-options "-s "$i"_metaphlan3.sam --bowtie2db /data/databases/MetaPhlAn3 -x mpa_v30_CHOCOPhlAn_201901 --bowtie2out "$i"_metaphlan3bowtie2out.txt" --threads 8 --memory-use maximum --output-basename "$i"
    rm MicrobialFastQ/"$i"_interleaved.fastq.gz
    module unload humann3/3.0
    module unload diamond/0.9.24

    module load samtools/1.10
    mv "$i"_metaphlan3.sam Metaphlan3outputs/Sam/"$i".sam
    samtools view -bS Metaphlan3outputs/Sam/"$i".sam -@ 8 > Metaphlan3outputs/Sam/"$i".bam
    rm Metaphlan3outputs/Sam/"$i".sam
    module unload samtools/1.10

    mv "$i"_humann3out/"$i"_genefamilies.tsv Humann3outputs/Genefamilies/"$i".tsv
    mv "$i"_humann3out/"$i"_pathabundance.tsv Humann3outputs/Pathabundance/"$i".tsv
    mv "$i"_humann3out/"$i"_pathcoverage.tsv Humann3outputs/Pathcoverage/"$i".tsv
    mv "$i"_humann3out/"$i"_humann_temp/"$i".log Humann3outputs/Logs/"$i".txt
    mv "$i"_humann3out/"$i"_humann_temp/"$i"_metaphlan_bugs_list.tsv Metaphlan3outputs/MPAoutput/"$i".tsv
    rm -r "$i"_humann3out/
done

module load humann3/3.0

humann_join_tables -i Humann3outputs/Pathabundance/ -o Humann3outputs/pathabundance.tsv -s --file_name '.tsv'
humann_join_tables -i Humann3outputs/Pathcoverage/ -o Humann3outputs/pathcoverage.tsv -s --file_name '.tsv'
humann_join_tables -i Humann3outputs/Genefamilies/ -o Humann3outputs/genefamilies.tsv -s --file_name '.tsv'
merge_metaphlan_tables.py Metaphlan3outputs/MPAoutput/*.tsv -o Metaphlan3outputs/metaphlan.tsv

humann_rename_table --input Humann3outputs/genefamilies.tsv -c /data/databases/Humann3/utility_mapping/map_uniref90_name.txt.gz -o Humann3outputs/genefamilies_uniref90_name.tsv
humann_renorm_table --input Humann3outputs/genefamilies_uniref90_name.tsv --units cpm --output Humann3outputs/genefamilies_uniref90_name_cpm.tsv
humann_split_stratified_table --input Humann3outputs/genefamilies_uniref90_name_cpm.tsv --output Humann3outputs/

humann_regroup_table -i Humann3outputs/genefamilies.tsv -c /data/databases/Humann3/utility_mapping/map_level4ec_uniref90.txt.gz --output Humann3outputs/genefamilies_level4ec.tsv -u N
humann_rename_table --input Humann3outputs/genefamilies_level4ec.tsv -c /data/databases/Humann3/utility_mapping/map_ec_name.txt.gz -o Humann3outputs/genefamilies_level4ec_name.tsv
humann_renorm_table --input Humann3outputs/genefamilies_level4ec_name.tsv --units cpm --output Humann3outputs/genefamilies_level4ec_name_cpm.tsv
humann_split_stratified_table --input Humann3outputs/genefamilies_level4ec_name_cpm.tsv --output Humann3outputs/

humann_regroup_table -i Humann3outputs/genefamilies.tsv -c /data/databases/Humann3/utility_mapping/map_go_uniref90.txt.gz --output Humann3outputs/genefamilies_go.tsv -u N
humann_rename_table --input Humann3outputs/genefamilies_go.tsv -c /data/databases/Humann3/utility_mapping/map_ko_name.txt.gz --output Humann3outputs/genefamilies_go_name.tsv
humann_renorm_table --input Humann3outputs/genefamilies_go_name.tsv --units cpm --output Humann3outputs/genefamilies_go_name_cpm.tsv
humann_split_stratified_table --input Humann3outputs/genefamilies_go_name_cpm.tsv --output Humann3outputs/

humann_renorm_table --input Humann3outputs/pathabundance.tsv --units cpm --output Humann3outputs/pathabundance_cpm.tsv
humann_split_stratified_table --input Humann3outputs/pathabundance_cpm.tsv --output Humann3outputs/

humann_renorm_table --input Humann3outputs/pathcoverage.tsv --units cpm --output Humann3outputs/pathcoverage_cpm.tsv
humann_split_stratified_table --input Humann3outputs/pathcoverage_cpm.tsv --output Humann3outputs/

module unload humann3/3.0
```
## Strain-level Analysis
You can use the SAM output files from [MetaPhlAn3](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0) or [HUMANn3](https://github.com/biobakery/humann) to perform strain-level profiling using [StrainPhlAn3](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-3.0)
```bash
mkdir StrainPhlAn/
mkdir StrainPhlAn/ConsensusMarkers/
mkdir StrainPhlAn/Output/
mkdir StrainPhlAn/CladeMarkers/

module load metaphlan2/3.0

for i in $(cat samplenames.txt)
do
    sample2markers.py -i Metaphlan3outputs/Sam/"$i".bam -o StrainPhlAn/ConsensusMarkers/ --nproc 24 -f bam
done

strainphlan -s StrainPhlAn/ConsensusMarkers/*.pkl -d /data/databases/MetaPhlAn3/mpa_v30_CHOCOPhlAn_201901.pkl --print_clades_only -o . --nproc 24 > strainphlan_clades.txt

for i in $(grep 's__' strainphlan_clades.txt | sed 's/.*s__/s__/' | sed 's/: .*//')
do
    extract_markers.py -c "$i" -d /data/databases/MetaPhlAn3/mpa_v30_CHOCOPhlAn_201901.pkl -o StrainPhlAn/CladeMarkers/
    strainphlan -s StrainPhlAn/ConsensusMarkers/*.pkl -m StrainPhlAn/CladeMarkers/"$i".fna -o StrainPhlAn/Output/ -c "$i" -d /data/databases/MetaPhlAn3/mpa_v30_CHOCOPhlAn_201901.pkl --nproc 24
done

module unload metaphlan2/3.0

rm -r StrainPhlAn/ConsensusMarkers/ StrainPhlAn/CladeMarkers/
```
## Metagenome Assembly using Metaspades
```bash
mkdir Metaspades_Assemblies/

for i in $(cat samplenames.txt)
do
    module load spades/3.13
    spades.py --meta -1 MicrobialFastQ/"$i"_microbial_R1.fastq.gz -2 MicrobialFastQ/"$i"_microbial_R2.fastq.gz --threads 24 -o Metaspades_Assemblies/"$i"_metaspades
    module unload spades/3.13

    mv Metaspades_Assemblies/"$i"_metaspades/contigs.fasta Metaspades_Assemblies/"$i"_contigs.fasta
    mv Metaspades_Assemblies/"$i"_metaspades/spades.log Metaspades_Assemblies/"$i"_metaspadeslog.txt
    rm -r Metaspades_Assemblies/"$i"_metaspades/
done
```
## Recovery of Metagenome-Assembled Genomes (MAGs)
 Metabat2
```bash
mkdir MetaBat2_Bins/

for i in $(cat samplenames.txt)
do
    module load bowtie2/2.3.4
    bowtie2-build Metaspades_Assemblies/"$i"_contigs.fasta Metaspades_Assemblies/"$i" -p 24
    bowtie2 -x Metaspades_Assemblies/"$i" -1 MicrobialFastQ/"$i"_microbial_R1.fastq.gz -2 MicrobialFastQ/"$i"_microbial_R2.fastq.gz -p 24 -S MetaBat2_Bins/"$i"_mapped.sam
    rm Metaspades_Assemblies/"$i"*.bt2
    module unload bowtie2/2.3.4

    module load samtools/1.10
    samtools view -bS MetaBat2_Bins/"$i"_mapped.sam -@ 24 > MetaBat2_Bins/"$i"_mapped.bam
    rm MetaBat2_Bins/"$i"_mapped.sam
    samtools sort MetaBat2_Bins/"$i"_mapped.bam -o MetaBat2_Bins/"$i"_mapped_sorted.bam -@ 24
    rm MetaBat2_Bins/"$i"_mapped.bam
    module unload samtools/1.10

    module load metabat2
    jgi_summarize_bam_contig_depths --outputDepth MetaBat2_Bins/"$i"_depth.txt MetaBat2_Bins/"$i"_mapped_sorted.bam
    metabat2 -i Metaspades_Assemblies/"$i"_contigs.fasta -a MetaBat2_Bins/"$i"_depth.txt -o MetaBat2_Bins/"$i"_bin -t 24
    module unload metabat2
    rm MetaBat2_Bins/"$i"_mapped_sorted.bam
    rm MetaBat2_Bins/"$i"_depth.txt
done

rename '_bin.'  '_bin_' MetaBat2_Bins/*.fa

module load checkm
checkm lineage_wf -t 24 -x fa MetaBat2_Bins/ CheckM_Output/ -f MetaBat2_Bins/checkm.txt
module unload checkm
rm -r CheckM_Output/
```
Alternatively, you can use [metaWRAP](https://github.com/bxlab/metaWRAP) to recover MAGs.  This is a wrapper that uses multiple different binning algorithms (in this case I use [MetaBat2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/), [CONCOCT](https://pubmed.ncbi.nlm.nih.gov/25218180/), and [MaxBin2](https://pubmed.ncbi.nlm.nih.gov/26515820/) but you can use any combination and number of binners), and combines and refines their outputs.  
In my (admittedly limited) testing it appears to recover more MAGs than MetaBat2 alone.
It can also perform quality trimming of raw metagenomic reads and removal of contamination (eg. human or bovine DNA) but here I only use it for recovering MAGs.  

The script below will
* create an output directory for each sample
* remove metagenomic contigs shorter than 1kbp using [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
* copy metagenomic data to the current working directory, decompress, and rename input files to match the format expected by metaWrap
* bin metagenomic contigs using MetaBat2, MaxBin2, and Concoct
* refine these bins and retain only those with >= 90% completeness and < 5% contamination as evaluated by CheckM
* reassemble these bins to improve assembly length
* tidy up the working directory by removing temporary and intermediary files

```bash
mkdir MetaWrap_Outputs/

for i in $(cat samplenames.txt)
do

    mkdir MetaWrap_Outputs/"$i"/

    module load bbmap/38.22
    reformat.sh in=Assemblies/"$i"_contigs.fasta out="$i"_filteredcontigs.fasta minlength=1000
    module unload bbmap/38.22

    cp FastQ/"$i"_R1.fastq.gz .
    gunzip "$i"_R1.fastq.gz
    mv "$i"_R1.fastq "$i"_1.fastq

    cp FastQ/"$i"_R2.fastq.gz .
    gunzip "$i"_R2.fastq.gz
    mv "$i"_R2.fastq "$i"_2.fastq

    module load metawrap/1.3.2
    module load concoct/1.1.0
    source activate concoct_env
    metawrap binning -o "$i"_binning -t 30 -a "$i"_filteredcontigs.fasta --metabat2 --maxbin2 --concoct "$i"_1.fastq "$i"_2.fastq
    conda deactivate
    module unload concoct/1.1.0
    module load checkm/1.0.18
    metawrap bin_refinement -o "$i"_binrefining -t 30 -A "$i"_binning/metabat2_bins/ -B "$i"_binning/maxbin2_bins/ -C "$i"_binning/concoct_bins/ -c 90 -x 5
    metawrap reassemble_bins -o "$i"_binreassembly -1 "$i"_1.fastq -2 "$i"_2.fastq -t 30 -m 800 -c 90 -x 5 -b "$i"_binrefining/metawrap_90_5_bins
    module unload checkm/1.0.18
    module unload metawrap/1.3.2
    
    mv "$i"_binreassembly/reassembled_bins/*.fa MetaWrap_Outputs/"$i"/
    mv "$i"_binreassembly/reassembled_bins.stats MetaWrap_Outputs/"$i"/
    rm -r "$i"_1.fastq "$i"_2.fastq "$i"_filteredcontigs.fasta "$i"_binning "$i"_binrefining "$i"_binreassembly
 
done


```


Using fastANI to get pairwise ANI values for all MAGs in a set
```bash
module load fastani/1.1

ls FNA/*.fa > fastani_querylist.txt
fastANI --ql fastani_querylist.txt --rl fastani_querylist.txt -o allMAGs_fastANIoutput.txt -t 10
rm fastani_querylist.txt

module unload fastani/1.1
```
Screening contigs (metagenomic or genomic) for genes of interest eg. antimicrobial resistance determinants, plasmid-associated genes, virulence factors
```bash
mkdir AbricateOutputs/

module load abricate/0.8

for i in $(cat samplenames.txt)
do
    abricate --db card "$i"_contigs.fasta --mincov 50 --threads 10 > AbricateOutputs/"$i"_card.txt
    abricate --db vfdb "$i"_contigs.fasta --mincov 50 --threads 10 > AbricateOutputs/"$i"_vfdb.txt
    abricate --db plasmidfinder "$i"_contigs.fasta --mincov 50 --threads 10 > AbricateOutputs/"$i"_plasmid.txt
done

module unload abricate/0.8
```
