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

-d: the pre-calculated features of sequence and epigenetic profile information for each sample in matrix format. (It should be noted that the input matrix sample must be in the same feature order with the given example).

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
  
 -d : the pre-calculated features for each sample in matrix format.
 
 -w : the size of the sliding window.(for example: -win 10 ,represent the 10 bins)
 
 -s ：the step of sliding windows across the whole genome.(for example: -slide 5, represent the 5 bins interval for adjacent regions)
 
 -p ：the smooth parameters for the optimized TAD boundary scores which it is between 0 and 1, default value is 0.5
 
 -r ： the resultion of bins which is consistent with the resultion of the input samples.(example: -res 100000, represent the 100kbp,equal 1 bin here)  
 
 -o : Output directory
 


use:

Rscript ./Scripts/run_pTADS.ori.R -m ./model/GM12878_model.RData -c ./model/GM12878_coeff1.RData -d ./example/test.chr1.40M_60M.matrix.txt -w 10 -s 1 -p 0.5 -r 100000 -o Results

two result files:

*.predicted.TAD_boundary.* :  pTADS predicted results

*.RF_BSSM.*    ： inclued the model prediction results of each sample,TAD boundary strength,Optimized Boundary Score judgment results,pTADS predicted results



Contacct us
If you have any questions or suggestions, please send an email to Yunlong Wang(yunlong@webmail.hzau.edu.cn).
