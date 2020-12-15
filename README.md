pTADS
# Introduction
pTADS(prediction of TAD boundary and strength) is a method to evaluate the performance of histone and transcription factor information in predicting the TAD boundaries and boundary strength across multiple cell lines. The pTADS method is consisted of random forest models to predict the TAD boundary and lasso based boundary score to characterize the TAD boundary strength, which is independent of the interaction matrix information of Hi-C.

# How to run it ?
pTADS is developed in R and can be downloaded from https://github.com/YunlongWang-ylw/pTADS. This repository contains scripts,examples and required packages for pTADS.

Scripts:
  run_pTADS.ori.R ;
  varSelRF.R ;

packages:
  Rscript, ggplot2, randomForest, caret, PRROC, pROC;
  
When you run the program, please Follow the README in the ./test directory.

# Required data
To run pTADS, the following data should be prepared:
-i1  the model of Random forest have been trained, Stored in an *.RData file 

-i2  Matrix data containing sample features (warning: The input matrix data, feature ID name and order shall be consistent with the matrix data in the example). 

# Example use: 
cd ./test

Rscript run_pTADS.ori.R -i1 ./model/GM12878.Pred.tP.score.2020.1.10.RData -i2 test.chr1.40M_60M.matrix.txt 10 1 0.5 100000 test_outfile


Contacct us
If you have any questions or suggestions, please send an email to Yunlong Wang(yunlong@webmail.hzau.edu.cn).
