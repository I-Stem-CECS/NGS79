# Synergism of dual AAV gene therapy and rapamycin rescues GSDIII phenotype in muscle and liver

Louisa Jauze, Mallaury Vie, Quentin Miagoux, Lucille Rossiaud, Patrice Vidal, Valle Montalvo-Romeral, Hanadi Saliba, Margot Jarrige, Helene Polveche, Justine Nozi, Pierre-Romain Le Brun, Luca Bocchialini, Amandine Francois, Jeremie Cosette, Jérémy Rouillon, Fanny Collaud, Fanny Bordier, Emilie Bertil-Froidevaux, Christophe Georger, Laetitia Van Wittenberghe, Adeline Miranda, Nathalie Daniele, David Gross, Lucile Hoch, Xavier Nissan & Giuseppe Ronzitti

DOI: [10.1172/jci.insight.172614](https://insight.jci.org/articles/view/172614)

## Description

- Team : [Pharmacologie des dystrophies musculaires](https://www.istem.eu/randd/pharmacologie-des-dystrophies-musculaires/) 
- Sequencing type : Quantseq
- Samples : 36
- Organism : Mouse GRCm38.99

## Objectives  

- Raw & normalized counts, 
- List of filtered DE genes (p-adj < 0.05 & |log2FC| >= 0.4 & BM >= 20). 

## Tools 

- FastQC ( v0.11.9 ) 
- cutadapt ( 1.18 ) 
- PRINSEQ-lite ( 0.20.4 ) 
- STAR ( 2.7.6a ) 
- samtools ( 0.1.18 )  
- HTSeq-count ( 0.12.4 ) 
- R Packages ( 4.1.2 ) 
( RColorBrewer, gplots, DESeq2, pheatmap , tidyverse, gridExtra, reshape, ggforce, cowplot, reshape2, viridis)

## DEG results 

### All tissus

|  | Down  | Up | Total |
| :------------------------- | :----:  | :----: | :----: |
| CTL WT vs CTL KO | 27 | 5 | 32 |
| CTL AAV KO vs CTL KO | 1 | 8 | 9 |
| KO AAV RAPA vs CTL KO | 8 | 0 | 8 |
| KO AAV RAPA vs CTL AAV KO | 0 | 0 | 0 |
| KO AAV RAPA vs RAPA KO | 3 | 3 | 6 |
| CTL WT vs RAPA KO | 85 | 18 | 103 |

### Only Triceps 

|  | Down  | Up | Total |
| :------------------------- | :----:  | :----: | :----: |
| CTL WT vs CTL KO | 181 | 24 | 205 |
| CTL AAV KO vs CTL KO | 0 | 0 | 0 |
| KO AAV RAPA vs CTL KO | 30 | 3 | 33 |
| KO AAV RAPA vs CTL AAV KO | 0 | 1 | 1 |
| KO AAV RAPA vs RAPA KO | 23 | 2 | 25 |
| CTL WT vs RAPA KO | 175 | 54 | 229 |

### Only Liver 

|  | Down  | Up | Total |
| :------------------------- | :----:  | :----: | :----: |
| CTL WT vs CTL KO | 21 | 11 | 32 |
| CTL AAV KO vs CTL KO | 24 | 78 | 102 |
| KO AAV RAPA vs CTL KO | 22 | 18 | 40 |
| KO AAV RAPA vs CTL AAV KO | 9 | 5 | 14 |
| KO AAV RAPA vs RAPA KO | 1 | 3 | 4 |
| CTL WT vs RAPA KO | 60 | 35 | 95 |


## Details samples 

- Size : 76 bp  
- Encoding : Sanger / Illumina 1.9 
- Date : 06/07/2022 

| Names | #Reads  | #Reads post-filters | STAR : % mapped reads  | # mapped reads | # Filtered mapped reads | STAR : % Filtered mapped reads |
| :----: | :--------:  | :--------: | :----:  | :--------: | :--------: | :----: |
| 1_S1 | 12 327 255 | 12 275 923 | 98,20 | 12 105 373 | 9 098 096 | 73,80 |
| 3_S2 | 12 248 834 | 12 165 186 | 96,68 | 11 842 510 | 9 145 699 | 74,67 |
| 5_S3 | 13 135 852 | 13 069 485 | 98,13 | 12 889 858 | 9 870 398 | 75,14 |
| 11_S4 | 10 787 406 | 10 721 206 | 97,75 | 10 544 518 | 8 207 843 | 76,09 |
| 12_S5 | 12 402 330 | 12 306 458 | 97,49 | 12 090 692 | 9 672 122 | 77,99 |
| 13_S6 | 14 490 263 | 14 415 162 | 98,07 | 14 210 186 | 11 231 611 | 77,51 |
| 14_S7 | 13 305 003 | 13 214 874 | 97,86 | 13 019 945 | 10 310 793 | 77,50 |
| 15_S8 | 13 571 362 | 13 487 144 | 97,92 | 13 289 141 | 10 437 726 | 76,91 |
| 16_S9 | 11 112 174 | 11 029 043 | 97,26 | 10 808 091 | 8 837 076 | 79,53 |
| 17_S10 | 12 867 347 | 12 781 855 | 97,64 | 12 564 316 | 9 940 610 | 77,25 |
| 18_S11 | 12 442 825 | 12 372 521 | 97,88 | 12 179 096 | 9 348 520 | 75,13 |
| 19_S12 | 12 637 899 | 12 570 230 | 98,00 | 12 385 242 | 9 476 576 | 74,99 |
| 20_S13 | 13 581 473 | 13 527 345 | 98,10 | 13 323 282 | 10 265 135 | 75,58 |
| 22_S14 | 13 439 393 | 13 381 456 | 98,25 | 13 203 804 | 10 135 109 | 75,41 |
| 23_S15 | 14 145 293 | 14 085 184 | 98,18 | 13 887 736 | 10 620 687 | 75,08 |
| 24_S16 | 14 067 552 | 13 987 881 | 97,64 | 13 735 175 | 10 386 196 | 73,83 |
| 25_S17 | 15 424 319 | 15 378 394 | 98,41 | 15 179 098 | 11 360 048 | 73,65 |
| 26_S18 | 13 782 515 | 13 713 954 | 98,10 | 13 521 150 | 10 338 003 | 75,01 |
| 47_S19 | 13 290 336 | 13 214 181 | 92,36 | 12 274 490 | 9 454 155 | 71,14 |
| 49_S20 | 12 805 545 | 12 706 900 | 94,00 | 12 036 977 | 9 532 876 | 74,44 |
| 51_S21 | 13 637 086 | 13 557 829 | 94,15 | 12 838 809 | 10 137 978 | 74,34 |
| 57_S22 | 13 799 433 | 13 708 904 | 92,82 | 12 808 373 | 10 056 585 | 72,88 |
| 58_S23 | 15 154 044 | 15 073 973 | 94,31 | 14 291 730 | 11 563 832 | 76,31 |
| 59_S24 | 13 988 141 | 13 915 512 | 93,85 | 13 128 221 | 10 566 973 | 75,54 |
| 60_S25 | 14 745 430 | 14 641 000 | 92,14 | 13 586 798 | 10 618 643 | 72,01 |
| 61_S26 | 13 822 985 | 13 749 308 | 93,36 | 12 904 543 | 10 282 872 | 74,39 |
| 62_S27 | 13 164 373 | 13 075 262 | 95,07 | 12 515 698 | 10 105 064 | 76,76 |
| 63_S28 | 13 016 480 | 12 895 035 | 94,54 | 12 305 752 | 10 208 008 | 78,42 |
| 64_S29 | 9 093 154 | 9 037 604 | 94,88 | 8 627 399 | 6 862 664 | 75,47 |
| 65_S30 | 16 848 317 | 16 732 560 | 93,89 | 15 818 384 | 12 520 297 | 74,31 |
| 66_S31 | 12 277 596 | 12 177 586 | 93,40 | 11 466 954 | 9 130 006 | 74,36 |
| 68_S32 | 12 136 515 | 12 062 141 | 93,38 | 11 333 254 | 8 949 317 | 73,74 |
| 69_S33 | 12 975 313 | 12 874 213 | 94,06 | 12 204 355 | 9 712 270 | 74,85 |
| 70_S34 | 13 278 078 | 13 201 935 | 93,50 | 12 415 317 | 9 711 055 | 73,14 |
| 71_S35 | 13 750 023 | 13 682 166 | 94,01 | 12 926 092 | 10 300 846 | 74,92 |
| 72_S36 | 12 817 356 | 12 739 836 | 94,04 | 12 054 075 | 9 621 997 | 75,07 |
