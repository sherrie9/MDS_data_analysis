# MDS_data_analysis
This repo showcases my work on a group project analyzing MDS data.

# Introduction
Myelodysplastic syndromes (MDS) are a heterogeneous group of hematopoietic stem cell (HSC) disorders characterized by defective haematopoiesis and cytopenia, associated with an increased risk of transformation to acute myeloid leukemia (AML). Various classification schemes exist to divide MDS patients into low-risk MDS and high-risk MDS for predicting the risk of transformation to AML, based on the percentage of blasts in the bone marrow or blood, karyotype, and number of cytopenia present [1, 2]. However, even low-risk MDS is associated with significant morbidity and mortality, with the median survival for low-risk MDS patients ranging from 5.3 to 8.8 years [3]. The median survival for high-risk MDS patients ranges from 0.8 to 1.6 years [3]. Thus, understanding molecular mechanisms involved specifically in low-risk MDS and high-risk MDS is urgent for development of better treatments for low-risk and high-risk MDS patients, respectively.

Recent evidence suggests that early MDS (low-risk MDS) is characterized by immunodeficiency, apoptosis, and chemokine signalling [4], whereas high-risk MDS showed upregulation of genes related to immature progenitor cells, with deregulated pathways involved in metabolism and DNA repair [5]. Pellagatti et al (2006) [6] demonstrated similarities of MDS gene expression profile to reported interferon-γ-induced gene expression in normal CD34+ cells, and the two most upregulated genes in MDS, are IFIT1 and IFITM1.

We used the raw reads (microarray data) generated by Pellagatti et al (2006) [6] to find differentially expressed genes (DEGs) and the resulting enriched functions in low-risk MDS and high-risk MDS using ErmineJ [7].

## My Work
Part of my contribution to the project was [data quality control](https://github.com/sherrie9/MDS_data_analysis/blob/master/Data_quality_control.md), and fitting a predictive model to distinguish lowrisk and high risk MDS bassed on differential gene expression, as shown in [predictiveModel.md](https://github.com/sherrie9/MDS_data_analysis/blob/master/predictiveModel.md).

## Poster
Poster is linked [here](https://github.com/sherrie9/Stat540_projects/blob/master/MDS_Poster%20(1).pdf)

Reference
---------

1.  Ades et al. Myelodysplastic syndromes, Lancet (2014) 383: 2239-2252. [link to the paper](https://doi.org/10.1016/S0140-6736(13)61901-7)

2.  Corey et al. Myelodysplastic syndromes: complexity of stem-cell diseases, Nature Reviews (2007) 7: 118-129. [link to the paper](https://www.nature.com/articles/nrc2047)

3.  MDS statistics from American Cancer Society. [link](https://www.cancer.org/cancer/myelodysplastic-syndrome/detection-diagnosis-staging/survival.html)

4.  Pellagatti et al. Deregulated gene expression pathways in myelodysplastic syndrome hematopoietic stem cells, Leukemia (2010) 24: 756-764. [link to the paper](https://www.nature.com/articles/leu201031)

5.  Shiozawa et al. Gene expression and risk of leukemic transformation in myelodysplsia, Blood (2017) 130: 2642-2653. [link to the paper](http://www.bloodjournal.org/content/130/24/2642.long?sso-checked=true)

6.  Pellagatti et al. Gene expression profiles of CD34+ cells in myelodysplastic syndromes: involvement of interferon-stimulated genes and correlation to FAB subtype and karyotype, Blood (2010) 108: 337-345. [link to the paper](http://www.bloodjournal.org/content/108/1/337.long)

7.  Lee et al. ErmineJ: Tool for functional analysis of gene expression datasets. BMC Bioinformatics (2005) 6: 269. [link to the paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-269)
