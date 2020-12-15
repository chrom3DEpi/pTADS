library(ggplot2)
library('getopt')
library(glmnet)
library(randomForest)
#source("varSelRF.R")
library(caret)
library(PRROC)
library(pROC)


command=matrix(c( 
  'help', 'h', 0,'logical', 'Help document',
  'file1', 'm', 1,'character', 'The RF model of cell training',
  'file2', 'c', 1,'character', 'LASSO coefficients of key features',
  'file3', 'd', 1, 'character', 'Matrix data containing key features',
  'win', 'w', 1, 'numeric', 'Defines the size of the sliding window, the number of bin.(example: -win 10 ,represent the 10 bins)',
  'slide', 's', 1, 'numeric', 'Defines the window slide distance.(example: -slide 1, represent the 1 bin)',
  'spar', 'p', 1, 'numeric', 'The smooth parameters of the curve are between 0 and 1',
  'res', 'r', 1, 'numeric', 'the size of bin, Match the size of the input matrix sample.(example: -res 100000, represent the 1 bin)',
  'output', 'o', 1, 'numeric', 'Output directory'),
  byrow=T,ncol=5)
  
args=getopt(command)
 
if (!is.null(args$help) || is.null(args$file1) || is.null(args$file2) || is.null(args$file3) ||is.null(args$win) || is.null(args$slide) || is.null(args$spar) || is.null(args$res) || is.null(args$output)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}



## Rscript run_pTADS.R test.chr1.40M_60M.matrix.txt 10 1 0.5 100000 test_outfile
#args   <- commandArgs(TRUE);
#inputfile<-args[1]; # test.chr1.40M_60M.matrix.txt
#winn <- args[2];
#sidee<- args[3];
#sparr<-args[4];
#ress<- args[5];
#outfilee<-args[6]; ## 

inputfile<-args$file3; # test.chr1.40M_60M.matrix.txt
winn <- args$win;
sidee<- args$slide;
sparr<-args$spar;
ress<- args$res;
outfilee<-args$output; ## 

winn<-as.numeric(as.character(winn))
sidee<-as.numeric(as.character(sidee))
sparr<-as.numeric(as.character(sparr))
ress<-as.numeric(as.character(ress))

########### step1: 读取测试文件
load(args$file1)
model<-model 
load(args$file2)
coeff1<-coeff1

# 提供 RF预测模型，  model
# 筛选特征拟合的lasso系数： coeff1+ nam 
# 
testX<-read.table(inputfile,head=T)


########### step2: 用训练集模型预测testX
rfPred <- predict(model, newdata = testX)
tP<-as.numeric(rfPred)-1

########### step3: lasso 系数计算全基因组等bin样本
coeff_15<-coeff1
name<-colnames(testX)
coe<-NULL;
List<-NULL;
for(i in 1:length(name)){
  p<-grep(name[i],name)
  coe<-c(coe,coeff_15[[p]])
  apart<-cbind(name[i],name[p],coeff_15[[p]])
  List<-rbind(List,apart) 
}

da<-testX

score<-NULL;
for(i in 1:dim(testX)[1])
 {
  for(j in 1:dim(testX)[2])
   {
     da[i,j]<-testX[i,j]*coe[j]
   }
   score<-c(score,sum(da[i,]))
}

########### step4: combine lasso 系数表征的边界强度+ 模型预测的边界
Pred<-cbind(rownames(testX),tP,score)
Pred<-as.data.frame(Pred)
colnames(Pred)<-c("sample","tP","score")
Pred[,1]<-as.character(Pred[,1])
Pred[,3]<-as.numeric(as.character(Pred[,3]))


########### step5:  窗口滑动，在每条染色体中范围内预测TAD边界

# head(Pred)
c<-unlist(strsplit(as.vector(Pred$sample),split=":"))
mat<-matrix(c,ncol=2,byrow=T)
name<-setdiff(unique(matrix(c,ncol=2,byrow=T)[,1]),c("chrX","chrY"))

########### step0: 参数

win=winn # 10bin=1M 以win滑动，每次滑动0.1M( bin),最后统计所有边界，去重
side=sidee # 1bin=0.1M 

for(ij in 1:length(name)){
  chr_Pred<-Pred[which(mat[,1]==name[ij]),]
  matt<-mat[which(mat[,1]==name[ij]),]
  
  bb<-chr_Pred
  ddd<-cbind(c(1:dim(bb)[1]),bb[,3])
  ddd<-as.data.frame(ddd)
  colnames(ddd)<-c("x","y")

  ## 去掉异常值
  x<-ddd[,1]
  y<-ddd[,2]
  y[which(y==0)]<-min(y[which(y>0)])
  yy<-1/y
  sort_yy<-sort(yy,decreasing = TRUE)
  outliers<-boxplot.stats(yy)$out  ## 异常值
  up<-outliers[which(outliers>median(yy))]
  down<-outliers[which(outliers<median(yy))]
  yy[which(yy>=sort_yy[which(sort_yy==min(up))+1])]<-sort_yy[which(sort_yy==min(up))+1]
  yy[which(yy<=sort_yy[which(sort_yy==max(down))-1])]<-sort_yy[which(sort_yy==max(down))-1]
  
  res=ress # bin=100k
  n<-ceiling((length(yy)-win)/side)+1
  List_pos<-NULL;
  List_bou<-NULL;
  for(ii in 1:n){
     e<-win+(ii-1)*side;
     label<- length(yy)-e
     if(label>=0){
       s<-e-win+1
       le<-e-s+1
       x_pos<-c(1:le)
       y_score<-yy[s:e]
      }
     if(label<0){
       s<-s+1
	   e<-length(yy)
       le<-e-s+1
       x_pos<-c(1:le)
       y_score<-yy[s:e]
      }
   
     ss<- smooth.spline(x_pos, y_score,spar=sparr)
     len<- length(ss$y)-1 
     slope<-NULL;
     for(i in 1:len){
        j<- i+1
        apart<- ss$y[j]-ss$y[i]
        slope<-c(slope,apart)
      }

     lenn<-length(slope)-1
     bou<-NULL;
     position<-NULL;
     start<-s # bin
     for(i in 1:lenn){
         j<-i+1
         if(slope[i]<0 & slope[j]>0 ){
	       position<-c(position,c(i,i+1,i+2))
		   #position<-position+start-1
		   if(start==1){
		     bin_e<- (i+2)*res
			 bin_s<- bin_e-res*3
		   }
		   if(start>1){
		     bin_e<- (start+i+2-1)*res
	         bin_s<- bin_e-res*3
		   }
           bou<-rbind(bou,c(name[ij],bin_s,bin_e))
         }  
     }
	 position<-position+start-1
	 predict_model<-sum(as.numeric(as.character(chr_Pred[position,2])))
	 if(predict_model>0)
     {List_pos<-c(List_pos,position)
      List_bou<-rbind(List_bou,bou)
	 }
  }
  
  
  # all select BSSM position
  out_pos<-unique(List_pos)
  resultt<-cbind(chr_Pred,rep("N",dim(chr_Pred)[1]))
  colnames(resultt)[4]<-c("BSSM")
  resultt[,4]<-as.character(resultt[,4])
  resultt[out_pos,4]<-"B"
  
  # filter
  #out_pos<-out_pos[which(resultt[out_pos,2]==1)]
  nn<-length(out_pos)-1
  arr<-NULL;
  FF<-NULL;
  for(i in 1:nn){
    j<-i+1
    p<-out_pos[j]-out_pos[i]
    if(p==1){
      arr<-c(arr,out_pos[i],out_pos[j])
    }
   if(p>1){
      arr<-c(arr,out_pos[i])
      arr<-unique(arr)
	  if(length(arr)>=4){
	      if(length(arr)%%2==0){
		      ne<-length(arr)/2
			  FF<-c(FF,arr[c(ne,ne+1)])
		  }
		  else{
		      ne<-floor(length(arr)/2)+1
			  FF<-c(FF,arr[ne])
		  }
	     #FF<-c(FF,arr[which(resultt[arr,3]==max(resultt[arr,3]))])
	   }
	  if(length(arr)<=3){
	     FF<-c(FF,arr)
	   }
	  arr<-NULL;
     }
  }
  
  arr<-c(arr,out_pos[i])
  arr<-unique(arr)
  if(length(arr)>=4){
	  if(length(arr)%%2==0){
		   ne<-length(arr)/2
		   FF<-c(FF,arr[c(ne,ne+1)])
	  }
	 else{
		   ne<-floor(length(arr)/2)+1
		   FF<-c(FF,arr[ne])
	  }
	}
  if(length(arr)<=3){
	     FF<-c(FF,arr)
	   }

 resultt<-cbind(resultt,rep("N",dim(chr_Pred)[1]))
 colnames(resultt)[5]<-c("pTADS")
 resultt[,5]<-as.character(resultt[,5])
 resultt[FF,5]<-"TAD boundary"
 
 options(scipen = 200)
 outfile1<-paste("./",outfilee,"/",outfilee,".",name[ij],".RF_BSSM.",win,".",side,".",sparr,".bed",sep="")
 write.table(resultt,file=outfile1,row.names=F,col.names=T,sep='\t',quote=F)

   ## TAD boundary region
   #out_bou<-matt[FF,]
   #c1<-unlist(strsplit(as.vector(out_bou[,2]),split="-"))
   #mat1<-matrix(c1,ncol=2,byrow=T)
   #out_bou<-cbind(out_bou[,1],mat1)
   out_bou<-List_bou[!duplicated(List_bou[,2]), ]

  for(iii in 1:dim(out_bou)[1]){
    out_bou[iii,2]<-format(as.numeric(as.character(out_bou[iii,2])),scientific=FALSE)
    out_bou[iii,3]<-format(as.numeric(as.character(out_bou[iii,3])),scientific=FALSE)
  }
  out_bou<-as.data.frame(out_bou)
  out_bou[,2]<-as.numeric(as.character(out_bou[,2]))
  out_bou[,3]<-as.numeric(as.character(out_bou[,3]))

  out_sort<-out_bou[order(out_bou[,2]),]
  sss<-as.numeric(as.character(unlist(strsplit(as.vector(matt[,2]),split="-"))[1]))
  out_sort[,2]<-out_sort[,2]+sss
  out_sort[,3]<-out_sort[,3]+sss
 
 predict_boundary_com<-NULL;
 label<-0
 for(jj in 1:c(dim(out_sort)[1]-1)){
 
   if(out_sort[jj,3]<out_sort[jj+1,2] & label==0){
     predict_boundary_com<-rbind(predict_boundary_com,out_sort[jj,])
	 label<-0
   }
   if(out_sort[jj,3]>=out_sort[jj+1,2] & label==0){
     s<-min(out_sort[jj,2],out_sort[jj+1,2])
	 e<-max(out_sort[jj,3],out_sort[jj+1,3])
	 label<-label+1
   }
   if(out_sort[jj,3]>=out_sort[jj+1,2] & label>0){
     s<-s 
	 e<-max(e,out_sort[jj+1,3])
	 label<-label+1
   }
   if(out_sort[jj,3]<out_sort[jj+1,2] & label>0){
     apart<-as.data.frame(cbind(as.character(out_sort[jj,1]),s,e))
	 colnames(apart)<-c("V1","V2","V3")
	 apart[,2]<-as.numeric(as.character(apart[,2]))
	 apart[,3]<-as.numeric(as.character(apart[,3]))
	 predict_boundary_com<-rbind(predict_boundary_com,apart)
	 label<-0
   }
   
 }

# last 1
if(out_sort[jj+1,2]>out_sort[jj,3] & label==0){
    predict_boundary_com<-rbind(predict_boundary_com,out_sort[jj+1,])
	label<-0
}
# last more 
if(out_sort[jj,3]>=out_sort[jj+1,2] & label>0){
  apart<-as.data.frame(cbind(as.character(out_sort[jj,1]),s,e))
  colnames(apart)<-c("V1","V2","V3")
  apart[,2]<-as.numeric(as.character(apart[,2]))
  apart[,3]<-as.numeric(as.character(apart[,3]))
  predict_boundary_com<-rbind(predict_boundary_com,apart)
  label<-0
}
  
predict_boundary_com[,2]<-as.numeric(as.character(predict_boundary_com[,2]))
predict_boundary_com[,3]<-as.numeric(as.character(predict_boundary_com[,3]))
  
   ## all predicted region as 300kbp
  bigger_position<-which(predict_boundary_com[,3]-predict_boundary_com[,2]>res*3)
  predict_boundary_com[bigger_position,2]<-(predict_boundary_com[bigger_position,2]+predict_boundary_com[bigger_position,3])/2 -(res*3/2)
  predict_boundary_com[bigger_position,3]<-(predict_boundary_com[bigger_position,2]+predict_boundary_com[bigger_position,3])/2 +(res*3/2)

  
 ## all predicted region as 200kbp
  #bigger_position<-which(predict_boundary_com[,3]-predict_boundary_com[,2] < res*2)
  #predict_boundary_com[bigger_position,2]<-predict_boundary_com[bigger_position,2] -(res/2)
 #predict_boundary_com[bigger_position,3]<-predict_boundary_com[bigger_position,3] +(res/2)
 
 #bigger_position<-which(predict_boundary_com[,3]-predict_boundary_com[,2] < res*3)
  #predict_boundary_com[bigger_position,2]<-(predict_boundary_com[bigger_position,2]+predict_boundary_com[bigger_position,3])/2 -(res*3/2)
  #predict_boundary_com[bigger_position,3]<-(predict_boundary_com[bigger_position,2]+predict_boundary_com[bigger_position,3])/2 +(res*3/2)
  
  outfile2<-paste("./",outfilee,"/",outfilee,".",name[ij],".predicted.TAD_boundary.",win,".",side,".",sparr,".bed",sep="")
  options(scipen = 200)
  write.table(predict_boundary_com,file=outfile2,row.names=F,col.names=F,sep='\t',quote=F)

}



