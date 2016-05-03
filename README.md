#CEGR Tools
Collection of java command-line utilities used in the Pugh Lab (Center for Eukaryotic Gene Regulation @ Penn State) and as part of 
the Chip-exo Galaxy flavor. Written by Will Lai. 

##Build Instructions
(after cloning this repo):
```
cd cegrtools;
ant makealljars;
```

The `makealljars` build target builds jar files for each of the tools described below, and outputs them into the `build/dist` directory. 

##BAMtoscIDX
Converts BAM data to ScIdx, the Strand-specific coordinate count format, which is used by 
tools within the Chip-exo Galaxy flavor. ScIdx files are 1-based. 

The format consists of 5 columns:
the chromosome, the position of the genomic coordinate, the number of tags on the forward strand, 
the number of tags on the reverse strand and the number of total tags on the position.

With pair-end reads, only the 5' end of READ1 will be used to create the ScIdx data file.
Tools that use this format include GeneTrack and MultiGPS.

##FourColorPlot
Produces a graphical representation of FASTA data with each nucleotide represented by a selected color.

##PEHistogram 
Creates a histogram of inferred fragment size lengths from mapped paired-end sequencing data.

##TagPileup
Creates a composite histogram of mapped sequencing tag 5' ends around a set of genomic coordinate points. 

##PairedendCrossPlot
Creates a cross-plot from paired end sequencing data. 

