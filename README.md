### Model building process
Scripts here are used to build RAVmodel_536, annotated with MSigDB C2.all gene sets.

| Script | Processes | Main Output | 
|----------|:-------------:|:-------------:|
| 01_import.R | Import `_quant.sf` files using `tximport` | `{studyName}.rds` |
| 02_Filtering.R | Log transformation | `{studyName}_count.csv` | 
| 03_Top_Genes.R | Order genes in each study based on their variance | topGenesInTrainingData |
| 04_Common_Genes.R | Select top-varing common genes | topGenes |
| 05_PCA.R | Calculate mean and sd for all samples/ rowNormalization/ PCA | `_PCs_rowNorm.rds` |
| 06_Clustering.R | Clustering/ Build RAVindex/ Calculate avg.silhouette.width | `_PCclusters.rds` |
| 07_Final_Model.R | MeSH/GSEA annotation and build RAVmodel | `_RAVmodel_{annotGeneSets}.rds` |


### A list of RAVmodels
Here is the summary table of studies/samples/genes for different models.

| # of studies | # of samples | % top varying | # of genes |
|----------|:-------------:|:-------------:|------:|
| 536 | 44,890 | 90% | 13,934 |
| 735 | 51,605 | 90% | 11,685 |
| 1,399 | 75,433 | 90% | 7,951 |
