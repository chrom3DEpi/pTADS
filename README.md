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

-m  the model of Random forest have been trained, Stored in an *.RData file 

-d  Matrix data containing sample features (warning: The input matrix data, feature ID name and order shall be consistent with the matrix data in the example). 

# Example use: 
input file：

./model/GM12878.Pred.tP.score.2020.1.10.RData  ： The RF model of GM12878 training, as well as the LASSO coefficients of key features

./example/test.chr1.40M_60M.matrix.txt    ：Matrix data containing key features

parameter：

Rscript ./Scripts/run_pTADS.ori.R -h

 -m : The RF model of GM12878 training, as well as the LASSO coefficients of key features
 
 -d : Matrix data containing key features
 
 -w : Defines the size of the sliding window, the number of bin.(example: -win 10 ,represent the 10 bins)
 
 -s  ：Defines the window slide distance.(example: -slide 1, represent the 1 bin)
 
 -p： The smooth parameters of the curve are between 0 and 1
 
 -r ： the size of bin, Match the size of the input matrix sample.(example: -res 100000, represent the 1 bin)  
 
 -o : Output directory
 


use:

Rscript ./Scripts/run_pTADS.ori.R -m ./model/GM12878.Pred.tP.score.2020.1.10.RData -d ./example/test.chr1.40M_60M.matrix.txt -w 10 -s 1 -p 0.5 -r 100000 -o Results

two result files:

*.predicted.TAD_boundary.* :  TAD boundaries are predicted

*.RF_BSSM.*    ： inclued the model prediction results of each sample,TAD boundary strength,Optimized Boundary Score judgment results,pTADS predicted results



Contacct us
If you have any questions or suggestions, please send an email to Yunlong Wang(yunlong@webmail.hzau.edu.cn).
