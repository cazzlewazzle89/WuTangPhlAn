# WuTangPhlAn
A work-in-progress (w)rapper script for Phylogenetic Analysis of metagenomic data

Mainly applicable for use on the hcux400 HPC at Teagasc but the general code should work on other machines.

Your working directory should contain
* a `RawFastQ/` directory with demultiplexed, lane-merged, paired-end metagenomic reads
  * named in the format `samplename_R[1-2].fastq.gz`
* a text file named `samplenames.txt` which lists each unique `samplename`, one per line.

For example, in an experiment examining two samples, a mock community, and a negative control - your RawFastQ/ directory would contain the following files

sampleA_R1.fastq.gz  
sampleA_R2.fastq.gz  
sampleB_R1.fastq.gz  
sampleB_R2.fastq.gz  
mock_R1.fastq.gz  
mock_R2.fastq.gz  
neg_R1.fastq.gz  
neg_R2.fastq.gz

and your samplenames.txt file would look like this

sampleA  
sampleB  
mock  
neg
