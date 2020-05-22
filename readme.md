interSummit_plotter for transcription factor binding sites based on ChIP-seq data
=======================


Data source
-----------
Gene Transcription Regulation Database (GTRD) is the most complete collection of uniformly processed ChIP-seq data to identify transcription factor binding sites for human and mouse.


Shiny app
-----------
This shiny app takes selected genes from the 709 transcription factors in Gene Transcription Regulation Database (GTRD) and their possible 6.7 billion pairwise intersummit distance, create boxplots, density plots and statistical tests results on the distribution of inter summit distances. 

This allows us to visualize and test if candidate transcription factors might be cofactors that work together based on the proximity of binding sites on the genome. Studies have found many transcription factors play coordinated role during embryo development and cell fate determination. For example, transcription factor Sp5 has been observed to have very close binding sites to transcription factor Oct4. 


Usage
-----------

- Example datasets is in *_intersummit.txt files. 
- Complete datasets (290,000 files) at biowulf.nih.gov:/data/CCBR_Pipeliner/db/PipeDB/db/GTRD/
- Input
  - Pair A:
    - gene1: transcription factor1
    - gene2: transcription factor2
  - Pair B:
    - gene3: transcription factor3
    - gene4: transcription factor4
- Output
  - boxplot
  - densityplot
  - statistical test

![workflow chart](https://github.com/da-yin/ccbr872_ChIPseq/UI.PNG)