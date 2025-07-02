# scRNAseq-Workshop
Workshop for summer bioinformatics day IBL.

## Table of Contents
- [1. Preparation](#1-preparation)
- [2. Data structure and metadata](#2-data-structure-and-metadata)
- [3. Merging data](#3-merging-data)
- [4. Preprocessing - Cell Filtering](#4-preprocessing---cell-filtering)
- [5. Preprocessing: More downstream steps](#5-preprocessing-more-downstream-steps)
- [6. Clustering](#6-clustering)
- [7. Annotation](#7-annotation)
- [8. DEG analysis](#8-deg-analysis)

## <ins>1. Preparation<ins>
<ins>1.1. Install R (https://posit.co/download/rstudio-desktop/)<ins>
- Windows: [https://download1.rstudio.org/electron/windows/RStudio-2025.05.1-513.exe](https://cran.rstudio.com/bin/windows/base/R-4.5.1-win.exe)
- Mac Silicon: https://cran.rstudio.com/bin/macosx/big-sur-arm64/base/R-4.5.1-arm64.pkg
- Mac Intel: https://cran.rstudio.com/bin/macosx/big-sur-x86_64/base/R-4.5.1-x86_64.pkg

<ins>1.2. Install RStudio<ins>
- Windows: https://download1.rstudio.org/electron/windows/RStudio-2025.05.1-513.exe
- Mac: https://download1.rstudio.org/electron/macos/RStudio-2025.05.1-513.dmg

<ins>1.3. Data<ins>
- Download all objects from the following link. Make sure you save it on a location that you can easily access later on. For this workshop it is recommended to just save it in a new folder on your desktop. \
- Download link:
 
  https://leidenuniv1-my.sharepoint.com/:f:/g/personal/s4409728_vuw_leidenuniv_nl/EvmAL6jF6lVElRbVv3sjLe8BodEvQ8P_qTiYzufhwhimhQ?e=ronQEs


## <ins>2. Data structure and metadata<ins>
Now that we have everything set up we can start working in RStudio.

<ins>Useful links:<ins>
- Zebrafish references:
  - https://daniocell.nichd.nih.gov/
  - https://zfin.org/
- Seurat Website:
    - https://satijalab.org/seurat/
- Seurat commands cheat sheet:
  - https://satijalab.org/seurat/articles/seurat5_essential_commands


#### <ins>2.1. Install packages in RStudio<ins>
```bash
install.packages("Seurat")
install.packages("BiocManager")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
```

#### <ins>2.2. Load libraries in RStudio<ins>
```bash
library(Seurat)
library(hdf5r)
library(ggplot2)
library(dplyr)
```

#### <ins>2.3. Read in the data<ins>
The count matrices that are saved as .h5 files. These count matrices can be generated using different methods, but for his experiment I used the 10X genomics CellRanger count() algorithm. Generating these matrices requires a reference genome so that gene names or symbols are correctly assigned. As this generation can take 30 minutes or longer to process on an HPC, we will skip this and use the downloaded matrices.

- Count matrix names:
  - DMSO_filtered_feature_bc_matrix.h5
  - DMSO_Amp_filtered_feature_bc_matrix.h5

First of all the working directory in RStudio needs to be changed to the location on your computer where the downloaded data is. In your file explorer go to the location and copy the full path from on top.
  ```bash
  # Set the working directory
  setwd("<your_path>")
  ```
Second, load in the count matrices
```bash
  # Load in the count matrices
  DMSO_mtx <- Read10X_h5(filename = "DMSO_filtered_feature_bc_matrix.h5")
  DMSO_Amp_mtx <- Read10X_h5(filename = "DMSO_Amp_filtered_feature_bc_matrix.h5")
```

#### <ins>2.4. Check the output of both matrices. How are they structured? What information can you derive?<ins>
  ```bash
  DMSO_mtx
  DMSO_Amp_mtx
  ```

- Additional commands you can apply on any object:
  - class()
  - head()
  - dim()
  - str()
  - rownames()
  - colnames()
    
> [!NOTE]
> **Question 1:** What is the class of the DMSO matrix?\
> **Question 2:** How many rows and columns do DMSO_mtx and DMSO_Amp_mtx have?\
> **Question 3:** What output does rownames() and colnames() give you?

<details>
<summary>Answers</summary>
Question 1: A dgCMatrix is a specific type of sparse matrix class in R. It only stores non-zero values, making it memory-efficient for datasets with many zeros such as scRNAseq. Using dgCMatrix saves significant memory compared to dense matrices.<br />
<br />
Question 2: DMSO_mtx has 25432 rows and 27489 columns. DMSO_Amp has 25432 rows and 24514 columns.<br />
<br />
Question 3: rownames() give the gene names. colnames() gives the cellular barcodes. These refer to the barcodes of a gel bead droplet to distinguish between cells.<br />
</details>
----------------------------------------------------------------------------------------------------------------------------------------

 
#### <ins>2.5. Create Seurat objects<ins>  

A Seurat object is a specialized data structure in R used for storing and analyzing scRNA-seq data. It's the core component of the Seurat package and is needed in order to work with the Seurat scRNAseq package.

This step contains an initial filtering step. It is not required to do this filtering yet but it will need to happen at some point anyway. Several parameters require explanation:
- __counts__ = the count matrix name as which you loaded in your data
- __project__ = the identity name you want to give all of the cells inside your object. This can for example be the type of treatment.
- __min.cells__ = this paramater will define the minimum amount of cells in which genes need to be expressed. Genes expressed in less cells than this value will be filtered out.
- __min.percent__ = the minimum amount of features each cell should contain. Cells with fewer genes will be filtered out
  
```bash
DMSO_seur <- Seurat::CreateSeuratObject(object = DMSO_mtx,
                                        project = "DMSO", 
                                        min.cells = 3, 
                                        min.features = 200)
  
DMSO_Amp_seur <- Seurat::CreateSeuratObject(DMSO_Amp_mtx,
                                          project = "DMSO_Amp", 
                                          min.cells = 3, 
                                          min.features = 200)
```

Run the following commands and see if you can find the answers to the following questions:
  ```bash
  DMSO_seur
  head(DMSO_seur)
  View(DMSO_seur@meta.data)
  
  DMSO_Amp_seur
  head(DMSO_Amp_seur)
  View(DMSO_Amp_seur@meta.data)
  ```


  > [!NOTE]
  > **Question 1:** How many cells (samples) and how many features (genes) does the DMSO and the DMSO_Amp samples have?\
  > **Question 2:** What is shown in the nCount_RNA and the nFeature_RNA columns?


<details>
<summary>Answers</summary>
Question 1: DMSO has 27438 samples (cells) and 23793 features (genes). DMSO_Amp has 24484 samples (cells) and  23699 features genes).<br />
<br />
Question 2: nCount_RNA shows the total amount of sequenced transcripts per cell. nFeature_RNA shows the total amount of genes per cell.<br />
<br />
</details>
----------------------------------------------------------------------------------------------------------------------------------------




## <ins>3. Merging data<ins>
Merging of data is at some point needed in order to intgrate the data (removal of batch effects and have similar biological cells overlap on UMAP plots). While it also allows for a better visualization, it is not always recommended to do early on. Merging data can conflict with the preprocessing which is why it is a better approach to first preprocess every sample individually and merge them after. For the sake of this workshop we will merge the objects as we cannot run the preprocessing anwyay in this timespan.

```bash
fish_merged <- merge(DMSO_seur, y = DMSO_Amp_seur, add.cell.ids = TRUE)
```

If you check this merged object out you will see that both samples are present in different count layers.
```bash
fish_merged
```

## <ins>4. Preprocessing - Cell Filtering<ins>
An important step in the cell filtering process is the removal of low quality cells such as:
- cells with high mitochondrial gene content
- doublets
- emtpy droplets

### <ins>4.1 Empty droplets<ins>
Empty droplets are droplets in which no cells are present. This does not mean that they are completely negative for transcripts as they can contain ambient RNA. There are several packages such as SoupX to identify these empty droplets. Because for this data the CellRanger count algorithm was used, empty droplets are already filtered out and we can proceed to the next step.

### <ins>4.2 High mitochondrial gene content<ins> 
Cells with a high percentage of mitochondrial genes expressed can resemble dying cells, stressed cells and/or damaged cells.  
It is important to filter these cells out as they are of low quality and can interfere with downstream analysis. Depending on your experiment you might use a different cut-off than here. Everything depends on the biological context in which you are working.  

As you can see in the metadata using View(obj@meta.data), There is no column showing the percentage of mitochondrial genes per cell.    
Depending on which organism you have data from, the mitochondrial genes have a different naming.
For *Homo sapiens* these are identified by **MT-**  
For this dataset of *Danio rerio* they follow the pattern **mt-**    

We can calculate the percentage of mitochondrial gene content per cell using the following command
  ```bash
  fish_merged[["percent.mt"]] <- PercentageFeatureSet(object = fish_merged, pattern = "^mt-")
  ```
We have now added this data to the metadata of our Seurat Object.
 <img width="663" alt="image" src="https://github.com/user-attachments/assets/86bf3dad-9657-401c-95cc-61796c00d911" />

A good habit during cell filtering and data analysis is to use visualization methods. This allows you to get a clear idea on what you are doing and will result in better decisions.
Let's visualize the total mRNA counts, total gene counts, and the mitochondrial percentage per cell using Violin Plots.

```bash
VlnPlot(fish_merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
```
<img src="https://github.com/user-attachments/assets/6b38bebe-ca0d-4780-890f-319f275873ef" width="450" height="400">


Using the subset() command we will filter out cells with percent.mt > 10%
```bash
fish_merged_sub <- subset(fish_merged, subset = percent.mt < 10)
VlnPlot(fish_merged_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
```
Check if the filtering succeeded by plotting a violin plot.<br />
<img src="https://github.com/user-attachments/assets/102a23dc-f6d2-4123-8ccd-0fc3190b7fcb" width="450" height="400">


### <ins>4.3 Doublet removal<ins>
Doublets are droplets in which 2 cells are present. This results in the both of these cells being labelled with the same cellular barcode. Heterotypical doublets contain 2 different cell types and are easier to detect than homotypical doublets. A lot of different packages for doublet detection exist and are constantly released. 

Here we use the package scDblFinder (https://github.com/plger/scDblFinder). Because this package is made for SingleCellExperiment (sce) objects and not Seurat objects, our Seurat objects need to be converted. For the sake of time you can skip the processing part, and dive into the visualization after doublet detection instead. 

Load the following objects into RStudio. These are processed with scDblFinder:
```bash
DMSO_doublets <- readRDS("DMSO_for_doublets_workshop.rds")
DMSO_Amp_doublets <- readRDS("DMSO_Amp_for_doublets_workshop.rds")
```
Using the table() command it is possible to get the amount of cells in a certain non-numerical metadata group
```bash
table(DMSO_doublets$scDblFinder.class.30k)
table(DMSO_Amp_doublets$scDblFinder.class.30k)
```

<ins>Questions:</ins>
1. How many singlets and doublets are there in the DMSO and DMSO_Amp sample?
2. Visualize nCount_RNA and nFeature_RNA for these objects using violin plots. Is there a difference in amount of genes/ the spread between doublets and singlets?  
   <ins>Note:</ins> Plotting singlets vs doublets requires an identity change
   ```bash
   Idents(DMSO_doublets) <- "scDblFinder.class.30k"
   Idents(DMSO_Amp_doublets) <- "scDblFinder.class.30k"
   ```


<details>
<summary>Click to see the Violin plots for the singlets and doublets</summary>
<img src="https://github.com/user-attachments/assets/7aef7615-b7cb-4faa-b04a-b042ed585ec7" width="600" height="400">
<img src="https://github.com/user-attachments/assets/4117afa4-3342-46b1-bb5d-98484d0e507a" width="600" height="400">

</details>

## <ins>5. Preprocessing: More downstream steps<ins>
After cell filtering in which low quality cells are removed, the next part of preprocessing is count normalization, data scaling, pca determination, integration to remove batch effects, identifying neighbors for clustering, clustering and dimensionality reduction.  
Because most of these steps require a lot of computational resources and time, we will not be doing them during this workshop and we will move to clustering instead.  
A standard workflow would however look like this:
  - ScaleData()
  - RunPCA()
  - IntegrateLayers()
  - JoinLayers()
  - FindNeighbors()
  - FindClusters()
  - RunUMAP()

## <ins>6. Clustering<ins>
Clustering is required to group similar cell types together so that:
- cell types can be identified
- rare cell types can be found
- differential gene expression can be ran between different groups of cells

In Seurat, the FindClusters() command is used. In here the resolution parameter determined the amount of clusters the dataset will be split in. A higher resolution means more clusters, a lower resolution less clusters. Let's run this command for different resolutions to find one that makes biological sense and allows us to clearly distinguish groups from each other so that they can be annotated.

### 6.1 Load in the fish_workshop.rds object that has been fully preprocessed.
```bash
fish <- readRDS("fish_workshop.rds")
```
### 6.2 Run FindClusters for different resolutions
```bash
  fish <- FindClusters(object = fish, 
               resolution = 10,  
               algorithm = 4, 
               random.seed = 123, 
               verbose = T
               )
  
  fish <- RunUMAP(object = fish, 
                  seed.use = 123, 
                  dims = 1:36, 
                  reduction = "harmony", 
                  reduction.name = "umap.harmony", 
                  reduction.key = "UMAPHARMONY_")

  DimPlot(fish, label = T)
```
Let's look at the result:  
<img src="https://github.com/user-attachments/assets/813f279f-d7ba-4e2d-867f-fef09cf0f33a" width="600" height="400">



With clustering resolution 10, we get 146 clusters. This, depending on the dataset not only does not make biological sense, it also does not make it easy to annotate the clusters.
Try to find a resolution that makes sense.  


> [!TIP]
> You only need to run the RunUMAP() command once.


> [!TIP]
> It is not needed to identify every possible cell type yet. One approach can be an initial clustering to distinguish tissues on a general level, and later on subcluster each tissue cluster again with another resolution to identify cell types per tissue.

<details>
<summary>My clustering resolution</summary>
I chose to run the first clustering with resolution 0.1. This so I could first annotate all the tissues. I afterwards subclustered each tissue with a different resolution that made the most biological sense to find the subtypes.<br />
<br />
<img src="https://github.com/user-attachments/assets/9171fee8-f601-4761-bebb-4a8bfc1273c0" width="600" height="400">

<img src="https://github.com/user-attachments/assets/bf28bddb-934b-45d7-b1e3-2a4d5fab114f" width="600" height="400">
</details>

## <ins>7. Annotation<ins>
Normally we would first annotate all our clusters before going into actual DEG analysis. While both steps use the same commands, the actual DEG analysis goes more in depth and is done with more parameters set in place.  

### 7.1 Annotation
Let's try to identify which tissue type a certain cluster is. Do the following:
- Set the Idents() of the fish object to "RNA_snn_res.0.1"
- Run the FindMarkers command and replace the ident.1 parameter with the cluster number
- Use online references such as DanioCell or Zfin to identify the tissue based on the genes that show up

  ```bash
  Idents(fish) <- "RNA_snn_res.0.1"
  markers <- FindMarkers(fish, ident.1 = 10, group.by = "RNA_snn_res.0.1")
  head(markers)
  ```
> [!NOTE]
> **Question 1:** What tissue type is cluster 6 related to?\
> **Question 2:** Can you find on Daniocell which cell type cluster 6 probably is?\
>
> **Question 3:** What tissue type is cluster 8 related to?\
> **Question 4:** Can you find on Daniocell which cell type cluster 8 probably is?\
> 
> **Question 5:** What tissue type is cluster 9 related to?\
> **Question 6:** Can you find on Daniocell which cell type cluster 9 probably is?

<img src="https://github.com/user-attachments/assets/91164f44-06e4-470e-a473-d326eb34e7d4" width="600" height="400">

<details>
<summary>Click to see the tissue annotations</summary>

<img src="https://github.com/user-attachments/assets/1068302f-b18a-4ddc-9ecb-6287703ee6d4" width="600" height="400">  
 
</details>


   
> [!IMPORTANT]
> If you end up getting an error saying that your memory limit was reached, perform the following steps<br/>
> 1. Run: library(usethis)
> 2. Run: usethis::edit_r_environ()
> 3. Add this to the first line of the R.environ tab that opened (replace 16Gb with the amount of RAM your laptop has): R_MAX_VSIZE=16Gb
> 4. Save and restart RStudio


## <ins>8. DEG analysis<ins>
Lets now do some DEG analysis. DEG or differential expressed gene analysis in single-cell RNA sequencing identifies genes that are expressed differently between groups of cells or conditions.

In order to run the DEG analysis between a tissue of the DMSO_Amp sample vs DMSO, we need to create an identifier in our data. Let's create a new metadata column that contains the sample name and the tissue name.
```bash
fish$sample.tissue <- paste0(fish$orig.ident, ".", fish$tissue_updated)
```
We can now check all groups in this column
```bash
table(fish$sample.tissue)
```

Lets run the following code for our DEG:
```bash
Idents(fish) <- "sample.tissue"
markers <- FindMarkers(fish, ident.1 = "DMSO.Muscle", ident.2 = "DMSO_Amp.Muscle", min.pct = 0.25, logfc.threshold = 0.25)
View(markers)
```

As you can see there are also genes present with high adjusted_p_values. These are not significant and we want to filter them out:
```bash
markers <- markers %>% filter(p_val_adj < 0.01)
View(markers)
```

Last, lets make a heatmap using 2 different packages.
First instead of plotting all genes, only take the top and bottom 20 differentially expressed genes. 
```bash
gene_order <- markers %>% arrange(desc(-avg_log2FC)) %>% pull(genes) %>% unique()
markers$genes <- factor(markers$genes, levels = gene_order)

top_genes <- top_n(markers, 20, avg_log2FC)
bottom_genes <- top_n(markers, -20, avg_log2FC)

markers_selection <- rbind(top_genes, bottom_genes)
```

#### ggplot2:
```bash
DEG_name <- paste0("DMSO_Amp vs DMSO")
p <- ggplot(markers_selection, aes(x = DEG_name, y = genes, fill = avg_log2FC)) +
  geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  ggtitle(paste0("DEGs for DMSO_Amp vs DMSO \nLog10 Harmony MAST") 
```
<img width="437" alt="image" src="https://github.com/user-attachments/assets/454cfd40-96e0-4a15-9c62-97de5bffb871" />


#### Complexheatmap:
This package does not make plots based on dataframes but on matrices that can only hold numeric values. We therefore need to tweak out dataframe a bit.
```bash
markers_mtx <- markers_selection[, c("genes", "avg_log2FC")]
markers_mtx
markers_mtx$genes <- NULL
markers_mtx <- as.matrix(markers_mtx)
```

Lets plot:
```bash
ComplexHeatmap::Heatmap(markers_mtx, 
                        width = unit(4, "cm"),
                        height = unit(10, "cm"),
                        row_names_gp = gpar(fontsize = 5),
                        column_names_gp = gpar(fontsize = 10),
                        column_names_rot = 45,
                        row_dend_width = unit(2, "cm"),
                        row_dend_gp = gpar(lwd = 0.3),
                        heatmap_legend_param = list(title = "avg_FoldChange"))
```

<img width="438" alt="image" src="https://github.com/user-attachments/assets/38a9eb23-45df-4775-99f8-7feb7fcbe297" />

# Thank you for your participation!
--------------------------------------------------------------------------------------------------
Thomas Schrauwen MSc. Molecular Genetics & Biotechnology
--------------------------------------------------------------------------------------------------






















