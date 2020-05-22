interSummit_plotter for transcription factor binding sites based on ChIP-seq data
=======================


Data source
-----------
Gene Transcription Regulation Database (GTRD) is a collection of uniformly processed ChIP-seq data to identify transcription factor binding sites for human and mouse.


Shiny app
-----------
This shiny app takes selected genes from the 709 transcription factors in Gene Transcription Regulation Database (GTRD) and their possible 6.7 billion pairwise intersummit distance, create boxplots, density plots and statistical tests results on the distribution of inter summit distances. 

This allows us to visualize and test if candidate transcription factors might be cofactors that work together based on the proximity of binding sites on the genome. Studies have found many transcription factors play coordinated role during embryo development and cell fate determination. For example, transcription factor Sp5 has been observed to have very close binding sites to transcription factor Oct4. 


Usage
-----------

- Example datasets is in *_intersummit.txt files. 
- Complete datasets (290,000 files) at biowulf.nih.gov:/data/CCBR_Pipeliner/db/PipeDB/db/GTRD/
- Input
  - Pair 1:
    - gene1: transcription factor A
    - gene2: transcription factor B
  - Pair 2:
    - gene3: transcription factor A
    - gene4: transcription factor B
- Output
  - boxplot
  - densityplot
  - statistical test (Kolmogorov–Smirnov test, is a nonparametric test of the equality of one-dimensional probability distributions)

![workflow chart](https://github.com/da-yin/ccbr872_ChIPseq/blob/master/UI.PNG)