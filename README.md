# Analysis of chromatin loop predictions

This repository contains the analysis of chromatin looping interaction 
predictions using the R package [sevenC](https://github.com/ibn-salem/sevenC). 

This repository documents all analysis in the associated publication:

***Computational Chromosome Conformation Capture by Correlation of ChIP-seq at CTCF motifs***
by Jonas Ibn-Salem and Miguel Andrade-Navarro

More details on using the sevenC R package itself can be found on the [sevenC website](https://ibn-salem.github.io/sevenC/).

# Order of analysis scripts

1. [download.sh](data/download.sh)
1. [data_prepare_candidates.R](R/data_prepare_candidates.R)
    1. [analyse_motif_significance.R](R/analyse_motif_significance.R)
    1. [analyse_motif_overlap.R](R/analyse_motif_overlap.R)
1. [analyse_selected_models.R](R/analyse_selected_models.R)
1. [analyse_selected_models_on_loop_subset.R](R/analyse_selected_models_on_loop_subset.R)
1. [screen_TFs_lfc.R](R/screen_TFs_lfc.R)
1. [binary_predictions.R](R/binary_predictions.R)
1. [compare_other_tools.R](R/compare_other_tools.R)
1. [loop_freaters_EDA.R](R/loop_freaters_EDA.R)
1. [analyse_HeLa.R](R/analyse_HeLa.R)
1. [analyse_input_data_types.R](R/analyse_input_data_types.R)
    - [analyse_input_data_types_performance.R](R/analyse_input_data_types_performance.R)
    - [analyse_input_data_types_features.R](R/analyse_input_data_types_features.R)

## Download external data

To download and process all external data for this analysis run the 
[download.sh](data/download.sh) script from within the `data` directory.
```
cd data
sh download.sh
cd ..
```

This will create sub-folders, download external data sets and formats them.
Note that the scripts downloads a lot of data, installs tools, and executes many 
external software, it might be easier to execute the commands in the script 
interactively step by step.

## Predictions using selected models

### Prepare candidates
First we need to read the CTCF motif data and prepare potential candidate loops as motif pairs. 
Motif pairs will be labeled whether they represent true loops according to Hi-C and ChIA-PET.
```
Rscript R/data_prepare_candidates.R
```

To run motif significance threshold analysis run:
```
Rscript R/analyse_motif_significance.R
```

## Screen TFs
In this step all TF ChIP-seq experiments from ENCODE in human GM12878 cells are analysed and compared by their performance in predicting chromatin interaction loops. 
```
Rscript R/screen_TFs_lfc.R
```

To run the predictions and performance evaluation using six selected TF ChIP-seq data sets, run the following R script. 
```
Rscript R/analyse_selected_models.R
```

## Loop feature analysis
Now we can run the analysis of loop features using the following script.
```
Rscript R/loop_freaters_EDA.R
```

## Predictions in HeLa cells
Here, we use a prediction model that was trained in GM12878 cells and evaluate its performance in human HeLa cells.
```
Rscript R/analyse_HeLa.R
```

## Input types for loop prediction
We analyse different input data types for ChIP-seq and other genomic assays for there performance in loop prediction.
```
Rscript R/analyse_input_data_types.R
Rscript R/analyse_input_data_types_performance.R
Rscript R/analyse_input_data_types_features.R
```


