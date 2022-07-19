# NGS79

## Description

- Team : [Pharmacologie des dystrophies musculaires](https://www.istem.eu/randd/pharmacologie-des-dystrophies-musculaires/) 
- Sequencing type : Quantseq
- Samples : 36
- Organism : Mouse GRCm38.99

- Analysis : 
-- Raw & normalized counts, 
-- List of filtered DE genes (p-adj < 0.05 & |log2FC| >= 0.4 & BM >= 20). 

## Tools 

- FastQC ( v0.11.9 ) 
- cutadapt ( 1.18 ) 
- PRINSEQ-lite ( 0.20.4 ) 
- STAR ( 2.7.6a ) 
- samtools ( 0.1.18 )  
- HTSeq-count ( 0.12.4 ) 

- R Packages ( 4.1.2 ) 
 - RColorBrewer
 - gplots
 - DESeq2
 - pheatmap 
 - tidyverse 
 - gridExtra 
 - reshape 
 - ggforce 
 - cowplot 
 - reshape2 
 - viridis 

## Details samples 

| Names | Taille des reads | dates QC | QC | Encoding qualité | Nombre de reads  | Filtre | QC post-filtre | % reads post-filtre | Nbre de reads post-filtre | STAR : % read alignés  | Nbre de reads mappé | Nbre de reads mappé filtrés | STAR : % read alignés filtré |
| :----: | :----: | :--------: | :--------: | :--------: | :--------: | :--------: | :----: | :----: | :--------: | :----:  | :--------: | :--------: | :----: |
| 1_S1 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 327 255 | CutNextera – R20 – M20 | ok | 0,42 | 12 275 923 | 98,20 | 12 105 373 | 9 098 096 | 73,80 |
| 3_S2 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 248 834 | CutNextera – R20 – M20 | ok | 0,68 | 12 165 186 | 96,68 | 11 842 510 | 9 145 699 | 74,67 |
| 5_S3 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 135 852 | CutNextera – R20 – M20 | ok | 0,51 | 13 069 485 | 98,13 | 12 889 858 | 9 870 398 | 75,14 |
| 11_S4 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 10 787 406 | CutNextera – R20 – M20 | ok | 0,61 | 10 721 206 | 97,75 | 10 544 518 | 8 207 843 | 76,09 |
| 12_S5 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 402 330 | CutNextera – R20 – M20 | ok | 0,77 | 12 306 458 | 97,49 | 12 090 692 | 9 672 122 | 77,99 |
| 13_S6 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 14 490 263 | CutNextera – R20 – M20 | ok | 0,52 | 14 415 162 | 98,07 | 14 210 186 | 11 231 611 | 77,51 |
| 14_S7 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 305 003 | CutNextera – R20 – M20 | ok | 0,68 | 13 214 874 | 97,86 | 13 019 945 | 10 310 793 | 77,50 |
| 15_S8 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 571 362 | CutNextera – R20 – M20 | ok | 0,62 | 13 487 144 | 97,92 | 13 289 141 | 10 437 726 | 76,91 |
| 16_S9 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 11 112 174 | CutNextera – R20 – M20 | ok | 0,75 | 11 029 043 | 97,26 | 10 808 091 | 8 837 076 | 79,53 |
| 17_S10 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 867 347 | CutNextera – R20 – M20 | ok | 0,66 | 12 781 855 | 97,64 | 12 564 316 | 9 940 610 | 77,25 |
| 18_S11 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 442 825 | CutNextera – R20 – M20 | ok | 0,57 | 12 372 521 | 97,88 | 12 179 096 | 9 348 520 | 75,13 |
| 19_S12 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 637 899 | CutNextera – R20 – M20 | ok | 0,54 | 12 570 230 | 98,00 | 12 385 242 | 9 476 576 | 74,99 |
| 20_S13 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 581 473 | CutNextera – R20 – M20 | ok | 0,40 | 13 527 345 | 98,10 | 13 323 282 | 10 265 135 | 75,58 |
| 22_S14 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 439 393 | CutNextera – R20 – M20 | ok | 0,43 | 13 381 456 | 98,25 | 13 203 804 | 10 135 109 | 75,41 |
| 23_S15 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 14 145 293 | CutNextera – R20 – M20 | ok | 0,42 | 14 085 184 | 98,18 | 13 887 736 | 10 620 687 | 75,08 |
| 24_S16 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 14 067 552 | CutNextera – R20 – M20 | ok | 0,57 | 13 987 881 | 97,64 | 13 735 175 | 10 386 196 | 73,83 |
| 25_S17 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 15 424 319 | CutNextera – R20 – M20 | ok | 0,30 | 15 378 394 | 98,41 | 15 179 098 | 11 360 048 | 73,65 |
| 26_S18 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 782 515 | CutNextera – R20 – M20 | ok | 0,50 | 13 713 954 | 98,10 | 13 521 150 | 10 338 003 | 75,01 |
| 47_S19 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 290 336 | CutNextera – R20 – M20 | ok | 0,57 | 13 214 181 | 92,36 | 12 274 490 | 9 454 155 | 71,14 |
| 49_S20 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 805 545 | CutNextera – R20 – M20 | ok | 0,77 | 12 706 900 | 94,00 | 12 036 977 | 9 532 876 | 74,44 |
| 51_S21 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 637 086 | CutNextera – R20 – M20 | ok | 0,58 | 13 557 829 | 94,15 | 12 838 809 | 10 137 978 | 74,34 |
| 57_S22 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 799 433 | CutNextera – R20 – M20 | ok | 0,66 | 13 708 904 | 92,82 | 12 808 373 | 10 056 585 | 72,88 |
| 58_S23 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 15 154 044 | CutNextera – R20 – M20 | ok | 0,53 | 15 073 973 | 94,31 | 14 291 730 | 11 563 832 | 76,31 |
| 59_S24 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 988 141 | CutNextera – R20 – M20 | ok | 0,52 | 13 915 512 | 93,85 | 13 128 221 | 10 566 973 | 75,54 |
| 60_S25 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 14 745 430 | CutNextera – R20 – M20 | ok | 0,71 | 14 641 000 | 92,14 | 13 586 798 | 10 618 643 | 72,01 |
| 61_S26 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 822 985 | CutNextera – R20 – M20 | ok | 0,53 | 13 749 308 | 93,36 | 12 904 543 | 10 282 872 | 74,39 |
| 62_S27 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 164 373 | CutNextera – R20 – M20 | ok | 0,68 | 13 075 262 | 95,07 | 12 515 698 | 10 105 064 | 76,76 |
| 63_S28 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 016 480 | CutNextera – R20 – M20 | ok | 0,93 | 12 895 035 | 94,54 | 12 305 752 | 10 208 008 | 78,42 |
| 64_S29 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 9 093 154 | CutNextera – R20 – M20 | ok | 0,61 | 9 037 604 | 94,88 | 8 627 399 | 6 862 664 | 75,47 |
| 65_S30 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 16 848 317 | CutNextera – R20 – M20 | ok | 0,69 | 16 732 560 | 93,89 | 15 818 384 | 12 520 297 | 74,31 |
| 66_S31 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 277 596 | CutNextera – R20 – M20 | ok | 0,81 | 12 177 586 | 93,40 | 11 466 954 | 9 130 006 | 74,36 |
| 68_S32 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 136 515 | CutNextera – R20 – M20 | ok | 0,61 | 12 062 141 | 93,38 | 11 333 254 | 8 949 317 | 73,74 |
| 69_S33 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 975 313 | CutNextera – R20 – M20 | ok | 0,78 | 12 874 213 | 94,06 | 12 204 355 | 9 712 270 | 74,85 |
| 70_S34 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 278 078 | CutNextera – R20 – M20 | ok | 0,57 | 13 201 935 | 93,50 | 12 415 317 | 9 711 055 | 73,14 |
| 71_S35 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 13 750 023 | CutNextera – R20 – M20 | ok | 0,49 | 13 682 166 | 94,01 | 12 926 092 | 10 300 846 | 74,92 |
| 72_S36 | 76 | 06/07/22 | PolyA / ok | Sanger / Illumina 1.9 | 12 817 356 | CutNextera – R20 – M20 | ok | 0,60 | 12 739 836 | 94,04 | 12 054 075 | 9 621 997 | 75,07 |

