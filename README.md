pTADS
# Introduction
The pTADS(prediction of TAD boundary and strength) method can simultaneously predict the TAD boundary and characterize the boundary strength across multiple cell lines by integrating sequence and epigenetic profile information such as histone and transcription factor binding information. The pTADS method is consisted of random forest models to predict the TAD boundary and lasso based boundary score to characterize the TAD boundary strength, which is independent of the contact matrix-based interaction matrix information of Hi-C.

# How to run it ?
The pTADS is developed in R and can be downloaded from https://github.com/YunlongWang-ylw/pTADS. This repository contains scripts,examples and required packages for pTADS.

Scripts:

  run_pTADS.ori.R ;


packages:

  Rscript, ggplot2, randomForest, caret, PRROC, pROC, getopt;
  
  
When you run the program, please Follow the README.txt in the ./pTADS directory.

# Required data
To run pTADS, the following data should be prepared:

-m  the pre-trained model of Random forest in  *.RData file 

-c  the coefficients for importance features in LASSO function in *.RData file

-d  Matrix data containing sample features (warning: The input matrix data, feature ID name and order shall be consistent with the matrix data in the example). 

# Example use: 
input file：

./model/GM12878_model.RData  ： The pre-trained random forest model for GM12878 cell line

./model/GM12878_coeff1.RData :  The coefficients for features in LASSO function for GM12878 cell line

./example/test.chr1.40M_60M.matrix.txt    ：The pre-calculated features of sequence and epigenetic profile information such as histone and 
  transcription factor binding information for the pre-trained model of Random forest and LASSO function

parameter：

Rscript ./Scripts/run_pTADS.ori.R -h

 -m : the pre-trained model of Random forest
  
 -c : the coefficients for importance features in LASSO function
  
 -d : Matrix data containing importance features
 
 -w : Defines the size of the sliding window .(example: -win 10 ,represent the 10 bins)
 
 -s  ：Defines the window slide distance.(example: -slide 1, represent the 1 bin)
 
 -p： The smooth parameters of the curve are between 0 and 1
 
 -r ： the size of bin, Match the size of the input matrix sample.(example: -res 100000, represent the 100kbp,equal 1 bin here)  
 
 -o : Output directory
 


use:

Rscript ./Scripts/run_pTADS.ori.R -m ./model/GM12878_model.RData -c ./model/GM12878_coeff1.RData -d ./example/test.chr1.40M_60M.matrix.txt -w 10 -s 1 -p 0.5 -r 100000 -o Results

two result files:

*.predicted.TAD_boundary.* :  TAD boundaries are predicted

*.RF_BSSM.*    ： inclued the model prediction results of each sample,TAD boundary strength,Optimized Boundary Score judgment results,pTADS predicted results



Contacct us
If you have any questions or suggestions, please send an email to Yunlong Wang(yunlong@webmail.hzau.edu.cn).
