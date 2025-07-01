# scRNAseq-Workshop
Workshop for summer bioinformatics day IBL

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
  ```bash
  https://leidenuniv1-my.sharepoint.com/:f:/g/personal/s4409728_vuw_leidenuniv_nl/EvmAL6jF6lVElRbVv3sjLe8BodEvQ8P_qTiYzufhwhimhQ?e=ronQEs
  ```

## 2. <ins>Data structure and metadata<ins>
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
install.packages("tidyverse")
```

#### <ins>2.2. Load libraries in RStudio<ins>
```bash
library(Seurat)
library(hdf5r)
library(ggplot2)
```

#### <ins>2.3. Read in the data<ins>
The count matrices that are saved as .h5 files. These count matrices can be generated using different methods, but for his experiment I used the 10X genomics CellRanger count() algorithm. Generating these matrices requires a reference genome so that gene names or symbols are correctly assigned. As this generation can take 30 minutes or longer to process on an HPC, I already prepared them.

- Count matrix names:
  - DMSO_filtered_feature_bc_matrix.h5
  - DMSO_Amp_filtered_feature_bc_matrix.h5
 
  ```bash
  # Set the working directory location to the folder where you stored the downloaded data.
  # Replace [Username] with your username on your pc.
  setwd("C:\Users\[Username]\Desktop")
  
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

- Questions:
  1. What is the class of the DMSO matrix?
  2. How many rows and columns do DMSO_mtx and DMSO_Amp_mtx have?
  3. What output does rownames() and colnames() give you?
 
#### <ins>2.5. Create a seurat objects. These are needed in order to work with the Seurat scRNAseq package<ins>  

A Seurat object is a specialized data structure in R used for storing and analyzing scRNA-seq data. It's the core component of the Seurat package.

This step contains an initial filtering step. It is not required to do this filtering yet but it will need to happen at some point anyway. Several parameters require explanation:
- __object__ = the object name as which you loaded in your data
- __project__ = the identity name you want to give all of the cells inside your object. This can for example be the type of treatment.
- __min.cells__ = this paramater will define the minimum amount of cells in which genes need to be expressed. Genes expressed in less cells than this value will be filtered out.
- __min.percent__ = the minimum amount of features each cell should contain. Cells with fewer genes will be filtered out
  
```bash
DMSO_seur <- Seurat::CreateSeuratObject(object = DMSO,
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

  Questions:
  1. How many cells (samples) and how many features (genes) does the DMSO and the DMSO_Amp samples have?
  2. What is shown in the nCount_RNA and the nFeature_RNA columns?


      
## <ins>3.Preprocessing: Cell Filtering<ins>
An important step in the cell filtering process is the removal of low quality cells such as :
- cells with high mitochondrial gene content
- doublets
- emtpy droplets

### <ins>3.1 Empty droplets<ins>
Empty droplets are droplets in which no cell is present. This does not mean that they are completely negative for mRNA as they can contain ambient RNA. There are several packages such as SoupX to identify these empty droplets. Because for this data the CellRanger count algorithm was used, empty droplets are already filtered out and we can proceed to the next step.

### <ins>3.2 High mitochondrial gene content<ins> 
Cells with a high percentage of mitochondrial genes expressed can resemble dying cells, stressed cells and/or damaged cells.  
It is important to filter these cells out as they are of low quality and can interfere with downstream analysis. Depending on your experiment you might use a different cut-off than here. Everything depends on the biological context in which you are working.  

As you can see in the metadata using View(obj@meta.data), There is no column showing the percentage of mitochondrial genes per cell.    
Depending on which organism you have data from, the mitochondrial genes have a different naming.
For * *Homo sapiens* * these are identified by **MT-**  
For this dataset of * *Danio rerio* * they follow the pattern **mt-**    

We can calculate the percentage of mitochondrial gene content per cell using the following command
  ```bash
  DMSO_seur[["percent.mt"]] <- PercentageFeatureSet(object = DMSO_seur, pattern = "^mt-")
  ```
We have now added this data to the metadata of our Seurat Object.
 <img width="663" alt="image" src="https://github.com/user-attachments/assets/86bf3dad-9657-401c-95cc-61796c00d911" />

A good habit during cell filtering and data analysis is to use visualization methods. This allows you to get a clear idea on what you are doing and will result in better decisions.
Let's visualize the total mRNA counts, total gene counts, and the mitochondrial percentage per cell using Violin Plots.

```bash
VlnPlot(DMSO_seur, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
```
<img src="https://github.com/user-attachments/assets/4155eb82-c892-412a-ad45-87e80d3ce341" width="450" height="400">

Using the subset() command we will filter out cells with mt.percent > 10%
```bash
DMSO_seur_sub <- subset(DMSO_seur, subset = percent.mt < 10)
VlnPlot(DMSO_seur_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
```
Check if the filtering succeeded  
<img src="https://github.com/user-attachments/assets/c0626fdf-e87d-48c9-8ad6-5f455d4659ab" width="450" height="400">

### <ins>3.3 Doublet removal<ins>
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

## <ins>3.Preprocessing: More downstream steps<ins>
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

## <ins>4.Clustering<ins>
Clustering is required to group similar cell types together so that:
- cell types can be identified
- rare cell types can be found
- differential gene expression can be ran between different groups of cells

In Seurat, the FindClusters() command is used. In here the resolution parameter determined the amount of clusters the dataset will be split in. A higher resolution means more clusters, a lower resolution less clusters. Let's run this command for different resolutions to find one that makes biological sense and allows us to clearly distinguish groups from each other so that they can be annotated.

### 4.1 Load in the fish_workshop.rds object that has been fully preprocessed.
```bash
fish <- readRDS("fish_workshop.rds")
```
### 4.2 Run FindClusters for different resolutions
```bash
  fish <- FindClusters(object = fish, 
               resolution = 5,  
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
<img src="https://github.com/user-attachments/assets/8aa81931-c5f0-47b6-8dde-7e2d99372444" width="600" height="400">


With clustering resolution 5, we get 96 clusters. This ,depending on the dataset not only does not make biological sense, it also does not make it easy to annotate the clusters.
Try to find a resolution that makes sense.  

> [!TIP]
> It is not needed to identify every possible cell type yet. One approach can be an initial clustering to distinguish tissues on a general level, and later on subcluster each tissue cluster again with another resolution to identify cell types per tissue.

## <ins>5.Annotation/DEG analysis<ins>
Normally we would first annotate all our clusters before going into actual DEG analysis. While both steps use the same commands, the actual DEG analysis goes more in depth and is done with more parameters set in place.  

### 5.1 Annotation
























