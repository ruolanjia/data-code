#psych包用于计算相关性、p值等信息

library(psych)

#pheatmap包用于绘制相关性热图

library(pheatmap)

#reshape2包用于输出数据的整合处理
library(ggplot2)
library(reshape2)

#读取微生物代谢物数据
phy <-read.csv(file = "D:/论文/1实验数据/牡蛎/代谢组/spearman系数热图/微生物.csv", row.names = 1)
met <-read.csv(file = "D:/论文/1实验数据/牡蛎/代谢组/spearman系数热图/代谢物.csv", row.names = 1)

#计算相关性矩阵（可选：”pearson”、”spearman”、”kendall”相关系数）、p值矩阵

cor <-corr.test(phy, met, method = "spearman",adjust= "none")

#提取相关性、p值

cmt <-cor$r

pmt <- cor$p

head(cmt)

head(pmt)

#输出相关系数表格,第一行为代谢物信息，第一列为物种信息

cmt.out<-cbind(rownames(cmt),cmt)

write.table(cmt.out,file= "D:/论文/1实验数据/牡蛎/代谢组/spearman系数热图/result/cor.txt",sep= "t",row.names=F)

#输出p值表格，第一行为代谢物信息，第一列为物种信息

pmt.out<-cbind(rownames(pmt),pmt)

write.table(pmt.out,file= "D:/论文/1实验数据/牡蛎/代谢组/spearman系数热图/result/pvalue.txt",sep= "t",row.names=F)

#以关系对的形式输出表格

#第一列为物种名，第二列为代谢物名，第三、第四列对应显示相关系数与p值

df <-melt(cmt,value.name= "cor")

df$pvalue <- as.vector(pmt)

head(df)

write.table(df,file= "D:/论文/1实验数据/牡蛎/代谢组/spearman系数热图/result/cor-p.txt",sep= "t")

#对所有p值进行判断，p< 0.01的以“**”标注，p值 0.01<p< 0.05的以“*”标注

if(!is.null(pmt)){
  
  ssmt <- pmt< 0.01
  
  pmt[ssmt] <- '**'
  
  smt <- pmt > 0.01& pmt < 0.05
  
  pmt[smt] <- '*'
  
  pmt[!ssmt&!smt]<- ''
  
} else{
  
  pmt <- F
  
}

#自定义颜色范围

mycol<-colorRampPalette(c("#68A9CF","white","#EC976A"))(800)
library(RColorBrewer)
#绘制热图,可根据个人需求调整对应参数

#scale=”none” 不对数据进行均一化处理 可选 "row", "column"对行、列数据进行均一化

#cluster_row/col=T 对行或列数据进行聚类处理，可选F为不聚类

#border=NA 各自边框是否显示、颜色，可选“white”等增加边框颜色

#number_color=”white” 格子填入的显著性标记颜色

#cellwidth/height=12 格子宽度、高度信息

p <- pheatmap(cmt,scale = "none",cluster_row = T, cluster_col = T, border= "white",
         
         display_numbers = pmt,fontsize_number = 12, number_color = "black",
         
         cellwidth = 20, cellheight =20,color=mycol)
