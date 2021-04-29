# WuTangPhlAn
A work-in-progress (w)rapper script for Phylogenetic Analysis of metagenomic data

Mainly applicable for use on the hcux400 HPC at Teagasc but the general code should work on other machines.

Your working directory should contain a `RawFastQ/` directory with demultiplexed, lane-merged, paired-end metagenomic reads named in the format `samplename_R[1-2].fastq.gz` and a text file named `samplenames.txt` which lists each unique `samplename`, one per line.

An example for these can be found in the `example_inputs` directory above.
