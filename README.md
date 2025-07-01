# scRNAseq-Workshop
Workshop for summer bioinformatics day IBL

## 1. Preparation
### 1.1. Install R (https://posit.co/download/rstudio-desktop/)
- Windows: [https://download1.rstudio.org/electron/windows/RStudio-2025.05.1-513.exe](https://cran.rstudio.com/bin/windows/base/R-4.5.1-win.exe)
- Mac Silicon: https://cran.rstudio.com/bin/macosx/big-sur-arm64/base/R-4.5.1-arm64.pkg
- Mac Intel: https://cran.rstudio.com/bin/macosx/big-sur-x86_64/base/R-4.5.1-x86_64.pkg

### 1.2. Install RStudio
- Windows: https://download1.rstudio.org/electron/windows/RStudio-2025.05.1-513.exe
- Mac: https://download1.rstudio.org/electron/macos/RStudio-2025.05.1-513.dmg

### 1.3. Data
- Download all objects from the following link: \
  https://leidenuniv1-my.sharepoint.com/:f:/g/personal/s4409728_vuw_leidenuniv_nl/EvmAL6jF6lVElRbVv3sjLe8BodEvQ8P_qTiYzufhwhimhQ?e=ronQEs


## 2. Data analysis
Now that we have everything set up we can start working in RStudio.

### Useful links:
- Zebrafish references:
  - https://daniocell.nichd.nih.gov/
  - https://zfin.org/
- Seurat Website:
    - https://satijalab.org/seurat/
- Seurat commands cheat sheet:
  - https://satijalab.org/seurat/articles/seurat5_essential_commands



### 2.1. Install packages in RStudio
```bash
install.packages("Seurat")
install.packages("BiocManager")
install.packages("tidyverse")
```

### 2.2. Load Libraries in RStudio
```bash
library(Seurat)
library(hdf5r)
library(ggplot2)
```

