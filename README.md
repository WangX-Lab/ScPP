# ScPP

A simple and effective algorithm for recognizing cell subpopulations with specific phenotypes based on the expression profiles of phenotype-associated marker genes in bulks and single cells.

<img width="1118" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/e77f43bf-15ca-4f7e-ba5e-d55dac3cc1af">



&nbsp;
&nbsp;
# Description

To infer phenotypes of single cells from scRNA-seq data, ScPP requires three types of data as input, including single-cell transcriptomes, bulk transcriptomes and phenotypic features in bulk data. The phenotypic features can be categorical variables, continuous variables or clinical survival.

&nbsp;
&nbsp;
# Details

+ The function `singcell_Preprocessing()` is for single-cell expression data preprocessing. Its input is a count matrix with cells being column names and genes row names.

+ The function `marker_Binary()` is used for binary variables to generate marker genes or signatures correlated with phenotype+ or phenotype- groups, containing 4 parameters: bulk_data, features, ref_group and Log2FC_cutoff.
  + "bulk_data" is a log2-normalized bulk expression data with genes in row and samples in column.
  + "features" is the feature data of bulk samples, column1 are sample names (colname is "Sample") and column2 are phenotype labels (colname is "Feature") of samples.
  + "ref_group" is a character indicating which feature is the phenotype- group.
  + "Log2FC_cutoff" is the absolute cutoff value of log2(fold change) (default: 0.585).

+ The function `marker_Continuous()` is used for continuous variables to generate marker genes or signatures correlated with the phenotype, containing 4 parameters: bulk_data, features, method and estimate_cutoff.
  + "bulk_data" is a log2-normalized bulk expression data with genes in row and samples in column.
  + "features" is the feature data of bulk samples, such as copy number variation values of each sample.
  + "method" is the method in cor.test(), with default "spearman" and the alternative is "pearson."
  + "estimate_cutoff" is the absolute cutoff value of correlation coefficient (default: 0.2).

+ The function `marker_Survival()` is used for survival data to generate marker genes or signatures correlated with patients'prognosis, containing 2 parameters: bulk_data and survival_data.
  + "bulk_data" is a log2-normalized bulk expression data with genes in row and samples in column.
  + "survival_data" is the survival data with time in column1 and status in column2, and its row names are sample ID.

+ The function `ScPP()` is used for Single Cells’Phenotype Prediction, containing 3 parameters: sc_dataset, geneList and probs.
  + "sc_dataset" is a seurat object of scRNA-seq data, which can be the output of the function `singcell_Preprocessing()`.
  + "geneList" is a gene list correlated with the phenotypes of interest, which can be the output of functions `marker_Binary()`, `marker_Continuous()` and `marker_Survival()` for binary variables, continuous varaibles and survival data, respectively.
  + "probs" is the α value of ScPP ranging from 0.2 to 0.5 (default: 0.2), which is the cutoff for selecting top and bottom α proportion of single cells ranked based on AUCell scores.

&nbsp;
&nbsp;

# Installation

- The `Seurat` package (*version* >= 4.3.0) is used for loading data and preprocessing.

- To install `ScPP`, first install the `AUCell` package （*version* >= 1.20.2）, if it is not already installed:
&nbsp;
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AUCell")
```

- Finally, users can install the released version of `ScPP` with:
&nbsp;
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("WangX-Lab/ScPP")
```

&nbsp;
&nbsp;
# Examples

&nbsp;
&nbsp;

## **Apply ScPP with binary variables** 

### **Prepare data**

```R
library(ScPP)
load(system.file("data/binary.RData",package = "ScPP"))
```



**Bulk data**

```R
bulk[1:6,1:6]
```

|        | TCGA-CA-5256-01 | TCGA-AZ-6599-01 | TCGA-AA-3655-01 | TCGA-A6-6137-01 | TCGA-CK-4952-01 | TCGA-A6-5657-01 |
| ------ | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| HIF3A  | 3.7172          | 2.3437          | 2.0858          | 6.0759          | 1.9506          | 5.4777          |
| CAMK4  | 3.0698          | 4.9331          | 2.3709          | 4.1387          | 1.1557          | 4.1746          |
| RNF112 | 1.3702          | 2.4817          | 2.4947          | 3.5941          | 2.3486          | 4.9185          |
| SPN    | 5.5207          | 5.6704          | 6.8577          | 8.0598          | 5.0049          | 7.6076          |
| LRRTM1 | 3.2408          | 1.6031          | 0.9465          | 1.9142          | 0               | 3.2523          |
| GRIN1  | 3.0698          | 6.4944          | 4.3225          | 2.8073          | 7.346           | 4.5             |



**Binary data**

```R
head(binary)
```

|      | Sample          | Feature |
| ---- | --------------- | ------- |
| 1    | TCGA-CA-5256-01 | Tumor   |
| 2    | TCGA-AZ-6599-01 | Tumor   |
| 3    | TCGA-AA-3655-01 | Tumor   |
| 4    | TCGA-A6-6137-01 | Tumor   |
| 5    | TCGA-CK-4952-01 | Tumor   |
| 6    | TCGA-A6-5657-01 | Tumor   |



**scRNA-count**

```R
sc_count[1:6,1:6]
```

|          | SMC01.T_AAACCTGCATACGCCG | SMC01.T_AAACCTGGTCGCATAT | SMC01.T_AAACCTGTCCCTTGCA | SMC01.T_AAACGGGAGGGAAACA | SMC01.T_AAACGGGGTATAGGTA | SMC01.T_AAAGATGAGGCCGAAT |
| -------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| A1BG     | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1BG-AS1 | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1CF     | 0                        | 2                        | 0                        | 0                        | 3                        | 0                        |
| A2M      | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2M-AS1  | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2ML1    | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |


&nbsp;


### Execute ScPP to select the informative cells

```R
sc = sc_Preprocess(sc_count)
geneList = marker_Binary(bulk, binary, ref_group = "Normal")
res = ScPP(sc, geneList)
head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey","blue","red"))

```

<img width="642" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/1808015c-3790-4a07-ba15-411496d42d22">





&nbsp;
&nbsp;


## **Apply ScPP with continuous variables**

### Prepare data

```R
library(ScPP)
load(system.file("data/continuous.RData",package = "ScPP"))
```



**Bulk data**

```R
bulk[1:6,1:6]
```

|        | TCGA-AZ-6599-01 | TCGA-AA-3655-01 | TCGA-A6-6137-01 | TCGA-CK-4952-01 | TCGA-A6-5657-01 | TCGA-AD-6963-01 |
| ------ | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| HIF3A  | 2.3437          | 2.0858          | 6.0759          | 1.9506          | 5.4777          | 4.4634          |
| CAMK4  | 4.9331          | 2.3709          | 4.1387          | 1.1557          | 4.1746          | 3.2363          |
| RNF112 | 2.4817          | 2.4947          | 3.5941          | 2.3486          | 4.9185          | 1.4621          |
| SPN    | 5.6704          | 6.8577          | 8.0598          | 5.0049          | 7.6076          | 7.396           |
| LRRTM1 | 1.6031          | 0.9465          | 1.9142          | 0               | 3.2523          | 0               |
| GRIN1  | 6.4944          | 4.3225          | 2.8073          | 7.346           | 4.5             | 3.1816          |



**Continuous data**

```R
head(continuous)
```

|      | samp            | TMB_non_silent |
| ---- | --------------- | -------------- |
| 1    | TCGA-AZ-6599-01 | 178            |
| 2    | TCGA-AA-3655-01 | 65             |
| 3    | TCGA-A6-6137-01 | 91             |
| 4    | TCGA-CK-4952-01 | 206            |
| 5    | TCGA-A6-5657-01 | 63             |
| 6    | TCGA-AD-6963-01 | 67             |



**scRNA-count**

```R
sc_count[1:6,1:6]
```

|          | SMC01.T_AAACCTGCATACGCCG | SMC01.T_AAACCTGGTCGCATAT | SMC01.T_AAACCTGTCCCTTGCA | SMC01.T_AAACGGGAGGGAAACA | SMC01.T_AAACGGGGTATAGGTA | SMC01.T_AAAGATGAGGCCGAAT |
| -------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| A1BG     | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1BG-AS1 | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1CF     | 0                        | 2                        | 0                        | 0                        | 3                        | 0                        |
| A2M      | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2M-AS1  | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2ML1    | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |



&nbsp;

### Execute ScPP to select the informative cells

```R
sc = sc_Preprocess(sc_count)
geneList = marker_Continuous(bulk, continuous$TMB_non_silent)
res = ScPP(sc, geneList)
head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey","blue","red"))
```

<img width="642" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/6c7dac7b-138b-41e2-9ef1-af87baff06e3">



&nbsp;
&nbsp;

## **Apply ScPP with  survival data**

### Prepare data

```R
library(ScPP)
load(system.file("data/survival.RData",package = "ScPP"))
```



**Bulk data**

```R
bulk[1:6,1:6]
```

|         | TCGA-69-7978 | TCGA-62-8399 | TCGA-78-7539 | TCGA-73-4658 | TCGA-44-6775 | TCGA-44-2655 |
| ------- | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ |
| HIF3A   | 4.2598       | 11.6239      | 9.1362       | 5.0288       | 4.0573       | 5.5335       |
| RTN4RL2 | 8.2023       | 5.5819       | 3.5365       | 7.4156       | 7.7107       | 5.3257       |
| HMGCLL1 | 2.7476       | 5.8513       | 3.8334       | 3.6447       | 2.9188       | 4.882        |
| LRRTM1  | 0            | 0.4628       | 4.7506       | 6.8005       | 7.7819       | 2.2882       |
| GRIN1   | 6.6074       | 5.4257       | 4.9563       | 7.351        | 3.5361       | 3.3311       |
| LRRTM3  | 1.7458       | 2.0092       | 0            | 1.4468       | 0            | 0            |



**Survival data**

```R
head(survival)
```

|              | status | time  |
| ------------ | ------ | ----- |
| TCGA-69-7978 | 0      | 4.4   |
| TCGA-62-8399 | 0      | 88.57 |
| TCGA-78-7539 | 0      | 25.99 |
| TCGA-73-4658 | 1      | 52.56 |
| TCGA-44-6775 | 0      | 23.16 |
| TCGA-44-2655 | 0      | 43.5  |



**scRNA-count**

```R
sc_count[1:6,1:6]
```

|          | SMC01.T_AAACCTGCATACGCCG | SMC01.T_AAACCTGGTCGCATAT | SMC01.T_AAACCTGTCCCTTGCA | SMC01.T_AAACGGGAGGGAAACA | SMC01.T_AAACGGGGTATAGGTA | SMC01.T_AAAGATGAGGCCGAAT |
| -------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| A1BG     | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1BG-AS1 | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A1CF     | 0                        | 2                        | 0                        | 0                        | 3                        | 0                        |
| A2M      | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2M-AS1  | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |
| A2ML1    | 0                        | 0                        | 0                        | 0                        | 0                        | 0                        |


&nbsp;


### Execute ScPP to select the informative cells

```R
sc = sc_Preprocess(sc_count)
geneList = marker_Survival(bulk, survival)
res = ScPP(sc, geneList)
head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey","blue","red"))


```

<img width="642" alt="image" src="https://github.com/WangX-Lab/ScPP/assets/54932820/47404d94-abe4-485c-8657-0f5e47bc62c3">


&nbsp;
# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
