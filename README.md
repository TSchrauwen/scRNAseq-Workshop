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
- Download both .h5 files from github and save in a folder on your desktop




## 2. Packages to install in Rstudio
```bash
# Toolkit for single cell genomics developed by Satija Lab
install.packages("Seurat")

# allows users to install and manage packages from the Bioconductor project wich contains open source software
install.packages("BiocManager")

# One of the many packages to detected doublets
BiocManager::install("scDblFinder")
```


## 3. Data analysis
Now that we have everything set up and the packages installed that we will use during this workshop, we can start looking into the data.






## 1. Install Python
- Windows: https://www.python.org/ftp/python/3.13.5/python-3.13.5-amd64.exe
  - (IMPORTANT: Check the box **Add python.exe to PATH** at the bottom in installation window)
  - In Powershell (press Windows key + R, type Powershell, press Enter):
    ```bash
    python --version # to verify installation
    ```
- Mac: https://www.python.org/ftp/python/3.13.5/python-3.13.5-macos11.pkg
    - In Terminal (press cmd + spacebar, type terminal, press Enter):
      ```bash
      python3 --version # to verify installation
      ```

## 2. Install Jupyter Lab
- Windows Powershell:
  ```bash
  pip install jupyterlab
  ```
- Mac Terminal:
  ```bash
  pip3 install jupyterlab
  ```

## 3. Install R
- Windows Powershell:
  ```bash
  winget show RProject.R
  ```


## 3. Load packages in Rstudio
```bash
library(Seurat)
library(scDblFinder)
```
