# Data preprocessing for TMOIA
Steps before running **TMOIA core**

## Preproccess datasets for KIPAN, BRCA, ROSMAP and LGG
```
# Load package
```
+ Pick samples' ID owning 3 omics including meth, rnaseqv2 and mirnaseq. Select columns by target features. Sort columns by ID (length = 16).
    ```
    filter_cols()
    ```
+ Set thresholds for var & NA and pick probes.
    ```
    filter_rows()
    ```
+ Download subtypes. (R lang)
+ Search labels by shared IDs. (R lang)
    ```
    get_label_func.R
    ```

## Preprocess multi-omic dataset of Arabidopsis thaliana
+ Get SNP matrix.
    ```
    python h5m2csv.py
    ```
+ Filter SNPs by a threshold variance.
    ```
    filter_rows_SNP()
    ```
+ Filter other omics.
    ```
    filter_ath()
    ```
+ Remove samples with missing omic.
    ```
    Rscript proc_dataset.R <trait(s)> <output dir> <path to omics> <...>
    ```
## Data splitting by Bootstrap
```
write_randomRepeats()
```
## Dimension reduction for SNP-gene & SNP-intergenic region using **local-connected** neural network
+ Search genes of SNPs.
    ```
    search_gene_of_SNP()
    ```
+ Train SNP dimension reduction models.
    ```
    run_DR() # on GPU
    ```
+ Calculate the "new omic" ***Genomic variation*** by trained models' weight.
    ```
    interpreter_inDir()
    ```