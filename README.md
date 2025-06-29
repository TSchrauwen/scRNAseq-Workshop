# scRNAseq-Workshop
Workshop for summer bioinformatics day IBL

## 1. Install R and Rstudio
- R:
  - Windows: https://cran.rstudio.com/bin/windows/base/R-4.5.1-win.exe
  - Mac Silicon: https://cran.rstudio.com/bin/macosx/big-sur-arm64/base/R-4.5.1-arm64.pkg
  - Mac Intel: https://cran.rstudio.com/bin/macosx/big-sur-x86_64/base/R-4.5.1-x86_64.pkg
 
- Rstudio:
  - Windows: https://download1.rstudio.org/electron/windows/RStudio-2025.05.1-513.exe
  - Mac: https://download1.rstudio.org/electron/macos/RStudio-2025.05.1-513.dmg
 
Check if R is properly installed:
1. Open Rstudio
2. In the terminal tab:
   ```bash
   R --version
   ```

## 2. Install packages in Rstudio
### 2.1. The Seurat package. An R toolkit for single cell genomics developed by Satija Lab
```bash
# R toolkit for single cell genomics developed by Satija Lab
install.packages("Seurat")

# allows users to install and manage packages from the Bioconductor project wich contains open source software
install.packages("BiocManager")

# One of the many packages to detected doublets
BiocManager::install("scDblFinder")
```


## 3. Load packages in Rstudio
```bash
library(Seurat)
library(scDblFinder)
```
