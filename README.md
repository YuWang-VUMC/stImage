# **stImage** 

StImage is a R package for integrated analysis of spatial transcriptomics and the corresponding modality images. StImage is freely available at https://github.com/YuWang-VUMC/stImage.



<img align="top" src="https://github.com/YuWang-VUMC/stImage/blob/master/man/figures/stImage_workflow.png" alt="drawing" width="600"/>
  
  
## Installation


Before installing stImage, dependencies should be installed first:

```r
# SPARK/SPARKX and SpatialPCA
library(devtools)
install_github('xzhoulab/SPARK')
install_github("shangll123/SpatialPCA")

#BioConductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rhdf5", "SingleCellExperiment", "BayesSpace", "omicade4"))

```

Users should also make sure to successfully install the *tensorflow* and *keras* following instruction [Tensorflow for R](https://tensorflow.rstudio.com/install/).

```r
#installation of tensorflow
install.packages("tensorflow")

#create python virtualenv. If already installed python, replace 'install_python()' with the path of excutable python.
library(reticulate)
path_to_python <- install_python()
virtualenv_create("r-reticulate", python = path_to_python)

library(tensorflow)
install_tensorflow(envname = "r-reticulate")

#installation of keras
install.packages("keras")
library(keras)
install_keras(envname = "r-reticulate")

#testing the installation
library(tensorflow)
tf$constant("Hello Tensorflow!")
```
Once Tensorflow and keras were successfully installed, you can then install the latest version of stImage from GitHub with:

```r
library(devtools)
install_github("YuWang-VUMC/stImage")
```

## Tutorial

-   [Adult Mouse Olfactory Bulb Visium Dataset processed by Space Ranger 2.0.0](https://htmlpreview.github.io/?https://github.com/YuWang-VUMC/stImage/blob/master/vignettes/MOB_workflow.html)

-   [HER2 breast cancer data by ST platform](https://htmlpreview.github.io/?https://github.com/YuWang-VUMC/stImage/blob/master/vignettes/HER2_workflow.html)

The tutorial includes main example codes for multiple spatial transcriptomics datasets (e.g. Adult Mouse Olfactory Bulb and Human breast tumor)


## License

stImage is licensed under the MIT License.

## Citation

stImage [xxx](xxx).

doi: xxx

