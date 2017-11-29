# Analysis of chromatin loop predictions

This repository contains the analysis of chromatin loopsing preidctions using 
the R package [chromloop](https://github.com/ibn-salem/chromloop). 
It aims to reproducibly document all analyis in the associated publicaton. 
More details on the method itself can be found here on the chromloop 
[website](https://ibn-salem.github.io/chromloop/).


# Order of analysis scripts
1. [download.sh](download.sh)
1. [data_prepare_candidates.R](R/data_prepare_candidates.R)
1. [analyse_motif_significance.R](R/analyse_motif_significance.R)
1. [screen_TFslfc.R](R/screen_TFslfc.R)
1. [analyse_selected_models.R](R/analyse_selected_models.R)
1. [loop_freaters_EDA.R](R/loop_freaters_EDA.R)
1. [analyse_HeLa.R](R/analyse_HeLa.R)


## Download external data

To download all external data for this analysis run the [download.sh](download.sh) script
from within the `data` directory.
```
cd data
sh download.sh
cd ..
```
This will create sub-folders, download external data sets and formats them.

## Predictions using selected models

### Prepare candidates
First we need to read the CTCF motif data and prepare candidate loops as motif pairs. 
Also candidates are labled whether they represent true loops according to Hi-C and ChIA-PET.
```
Rscript R/data_prepare_candidates.R
```

To run moitf significance threshold analysis run:
```
Rscript R/analyse_motif_significance.R
```

## Screen TFs
In this step all TF ChIP-seq experiments from ENCODE in human GM12878 cells are analysed and compared for ther perfrmance in predicting chromatin interaction loops. 
```
Rscript R/screen_TFslfc.R
```

To run the predictions and performance evaluation using six selected TF ChIP-seq data sets, run the follwoing R script. 
```
Rscript R/analyse_selected_models.R
```

## Exploratory analysis
Now we can run the exploratory analysis of loop featers using the following script.
```
Rscript R/loop_freaters_EDA.R
```

## Predictions in HeLa cells
Here, we use a prediction model that was trained in GM12878 cells and evaluate its performance in human HeLa cells.
```
Rscript R/analyse_HeLa.R
```

