# Primaquine-Challenge

This repository provides data and code for the analysis of a phase 1 healthy volunteer study of primaquine in male thai/burmese G6PD deficient individuals. The study was pre-registered as two separate sub-studies (ascending dose primaquine and single high dose primaquine: Thai Clinical Trial Registry numbers TCTR20170830002 and TCTR20220317004).


The preprint with the clinical results (under review at elife): [Pharmacometric assessment of primaquine induced haemolysis in glucose-6-phosphate dehydrogenase deficiency](https://www.medrxiv.org/content/10.1101/2023.02.24.23286398v1).

We are still working on the paper for the red blood cell modelling results (update of a previous paper in elife: [Modelling primaquine-induced haemolysis in G6PD deficiency](https://elifesciences.org/articles/23061v1)).

## Data and SAP for clinical paper

The are data are provided in the folder *Data*. The statistical analysis plan is given in the pdf *PQ_SAP_v1.pdf*.


## Main result

The results in the clinical paper are generated in the Markdown *Main_analysis.qmd*.


The data from the ascending dose study:
![](Main_analysis_files/figure-html/Fig1-1.png) 



The data from the single high dose study:
![](Main_analysis_files/figure-html/Fig2-1.png) 

## Modelling work

This is still ongoing. The main models written in stan are in the folder *Stan_models*. 
