# Enter the Wu-TangPhlAn
The goal is to turn this into a (w)rapper script for **Ph**y**l**ogenetic **An**alysis of metagenomic data.  
At the moment it is a work-in-progress tutorial mainly applicable for use on the hcux400 HPC at Teagasc but the general code should work on other machines.

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

Adapter removal and quality trimming using TrimGalore
```
mkdir TrimmingReports/
mkdir FastQC/
mkdir TrimmedFastQ/
```
```bash
for i in $(cat samplenames.txt)
do
    module load cutadapt/2.6
    module load fastqc/0.11.8
    ~/TrimGalore-0.6.0/trim_galore --paired RawFastQ/"$i"_R1.fastq.gz RawFastQ/"$i"_R2.fastq.gz --fastqc -j 24 -q 30 -o "$i"_trimout
    module load cutadapt/2.6
    module unload fastqc/0.11.8

    mv "$i"_trimout/"$i"_R1_val_1.fq.gz TrimmedFastQ/"$i"_trimmed_R1.fastq.gz
    mv "$i"_trimout/"$i"_R2_val_2.fq.gz TrimmedFastQ/"$i"_trimmed_R2.fastq.gz
    mv "$i"_trimout/*_trimming_report.txt TrimmingReports/
    mv "$i"_trimout/*.html FastQC/

    rm -r "$i"_trimout
done
```
Host contamination removal using Bowtie2
```
mkdir MicrobialFastQ/
```
```
for i in $(cat samplenames.txt)
do
    module load bowtie2/2.3.4
    bowtie2 -x /data/databases/hostremoval/Homo_sapiens/Bowtie2/Homo_sapiens -1 TrimmedFastQ/"$i"_trimmed_R1.fastq.gz -2 TrimmedFastQ/"$i"_trimmed_R2.fastq.gz -S "$i".sam -p 24
    module unload bowtie2/2.3.4

    module load samtools/1.10
    samtools view -bS "$i".sam -@ 24 > "$i".bam
    rm "$i".sam
    samtools view -b -f 12 -F 256 "$i".bam -@ 24 > "$i"_nonhost.bam
    rm "$i".bam
    samtools sort -n "$i"_nonhost.bam -T "$i"_nonhost_temp -o "$i"_nonhost_sorted.bam -@ 24
    rm "$i"_nonhost.bam
    module unload samtools/1.10

    module load bedtools/2.27.1
    bedtools bamtofastq -i "$i"_nonhost_sorted.bam -fq MicrobialFastQ/"$i"_microbial_R1.fastq -fq2 MicrobialFastQ/"$i"_microbial_R2.fastq
    rm "$i"_nonhost_sorted.bam
    module unload bedtools/2.27.1

    gzip MicrobialFastQ/"$i"_microbial_R1.fastq MicrobialFastQ/"$i"_microbial_R2.fastq
done
```
Species-level and functional profiling using HUMANn3
```
mkdir Metaphlan3outputs/
mkdir Metaphlan3outputs/Sam/
mkdir Metaphlan3outputs/MPAoutput/
mkdir Humann3outputs/
mkdir Humann3outputs/Genefamilies/
mkdir Humann3outputs/Pathabundance/
mkdir Humann3outputs/Pathcoverage/
mkdir Humann3outputs/Logs/
```
```
for i in $(cat samplenames.txt)
do
    module load bbmap/38.22
    reformat.sh in=MicrobialFastQ/"$i"_microbial_R1.fastq.gz in2=MicrobialFastQ/"$i"_microbial_R2.fastq.gz out=MicrobialFastQ/"$i"_interleaved.fastq.gz
    module unload bbmap/38.22

    module load humann3/3.0
    module load diamond/0.9.24
    humann -i MicrobialFastQ/"$i"_interleaved.fastq.gz -o "$i"_humann3out --nucleotide-database /data/databases/Humann3/chocophlan --protein-database /data/databases/Humann3/uniref --search-mode uniref90 --metaphlan-options "-s "$i"_metaphlan3.sam --bowtie2db /data/databases/MetaPhlAn3 -x mpa_v30_CHOCOPhlAn_201901 --bowtie2out "$i"_metaphlan3bowtie2out.txt" --threads 24 --memory-use maximum --output-basename "$i"
    rm MicrobialFastQ/"$i"_interleaved.fastq.gz
    module unload humann3/3.0
    module unload diamond/0.9.24

    module load samtools/1.10
    mv "$i"_metaphlan3.sam Metaphlan3outputs/Sam/"$i".sam
    samtools view -bS Metaphlan3outputs/Sam/"$i".sam -@ 24 > Metaphlan3outputs/Sam/"$i".bam
    rm Metaphlan3outputs/Sam/"$i".sam
    module unload samtools/1.10

    mv "$i"_humann3out/"$i"_genefamilies.tsv Humann3outputs/Genefamilies/"$i".tsv
    mv "$i"_humann3out/"$i"_pathabundance.tsv Humann3outputs/Pathabundance/"$i".tsv
    mv "$i"_humann3out/"$i"_pathcoverage.tsv Humann3outputs/Pathcoverage/"$i".tsv
    mv "$i"_humann3out/"$i"_humann_temp/"$i".log Humann3outputs/Logs/"$i".txt
    mv "$i"_humann3out/"$i"_humann_temp/"$i"_metaphlan_bugs_list.tsv Metaphlan3outputs/MPAoutput/"$i".tsv
    rm -r "$i"_humann3out/
done
```
```
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

humann_regroup_table -i Humann3outputs/genefamilies.tsv -c /data/databases/Humann3/utility_mapping/map_ko_uniref90.txt.gz --output Humann3outputs/genefamilies_ko.tsv -u N
humann_rename_table --input Humann3outputs/genefamilies_ko.tsv -c /data/databases/Humann3/utility_mapping/map_ko_name.txt.gz --output Humann3outputs/genefamilies_ko_name.tsv
humann_renorm_table --input Humann3outputs/genefamilies_ko_name.tsv --units cpm --output Humann3outputs/genefamilies_ko_name_cpm.tsv
humann_split_stratified_table --input Humann3outputs/genefamilies_ko_name_cpm.tsv --output Humann3outputs/

humann_renorm_table --input Humann3outputs/pathabundance.tsv --units cpm --output Humann3outputs/pathabundance_cpm.tsv
humann_split_stratified_table --input Humann3outputs/pathabundance_cpm.tsv --output Humann3outputs/

humann_renorm_table --input Humann3outputs/pathcoverage.tsv --units cpm --output Humann3outputs/pathcoverage_cpm.tsv
humann_split_stratified_table --input Humann3outputs/pathcoverage_cpm.tsv --output Humann3outputs/

module unload humann3/3.0
```
Strain-level phylogenetic analysis using StrainPhlAn3
```
mkdir StrainPhlAn/
mkdir StrainPhlAn/ConsensusMarkers/
mkdir StrainPhlAn/Output/
mkdir StrainPhlAn/CladeMarkers/
```
```
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
Metagenome Assembly using Metaspades
```
mkdir Metaspades_Assemblies/
```
```
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
MAG recovery using Metabat2
```
mkdir MetaBat2_Bins/
```
```
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
Using fastANI to get pairwise ANI values for all MAGs in a set
```
module load fastani/1.1

ls FNA/*.fa > fastani_querylist.txt
fastANI --ql fastani_querylist.txt --rl fastani_querylist.txt -o allMAGs_fastANIoutput.txt -t 10
rm fastani_querylist.txt

module unload fastani/1.1
```
Screening contigs (metagenomic or genomic) for 
```
mkdir AbricateOutputs/
```
```
module load abricate/0.8

for i in $(cat samplenames.txt)
do
    abricate --db card "$i"_contigs.fasta --mincov 50 --threads 10 > AbricateOutputs/"$i"_card.txt
    abricate --db vfdb "$i"_contigs.fasta --mincov 50 --threads 10 > AbricateOutputs/"$i"_vfdb.txt
    abricate --db plasmidfinder "$i"_contigs.fasta --mincov 50 --threads 10 > AbricateOutputs/"$i"_plasmid.txt
done

module unload abricate/0.8
```
