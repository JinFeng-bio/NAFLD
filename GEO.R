BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("RColorBrewer")
BiocManager::install("impute")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")

setwd("C:/J's/liver")
#����R��
library(affyPLM)
#��������
Data<-ReadAffy()
#�����ݼ����лع����
Pset<-fitPLM (Data)
#�������ƣ��鿴�Ҷ�ͼ
image(Data[,1])
#���ݼ���������Ȩ��ͼ
image(Pset,type="weights",which=1,main="Weights")
#���ݼ����������в�ͼ
image(Pset,type="resids",which=1,main="Residuals")
#���ݼ����������в����ͼ
image(Pset,type="sign.resids",which=1,main="Residuals.sign")

#��������:��Զ������RLE��
library(RColorBrewer)
#������ɫ
colors<-brewer.pal(12,"Set3")
#����RLE����ͼ
par(mar=c(10.1,4.1,4.1,2.1))
Mbox(Pset,col=colors,main="RLE",las=3,cex.axis=0.8)

#����NUSE����ͼ
boxplot(Pset,col=colors,main="NUSE",las=3)

#�������ƣ�RNA����ͼ
library(affy)
#��ȡ��������
data.deg<-AffyRNAdeg(Data)
#����RNA����ͼ
plotAffyRNAdeg(data.deg,col=colors)
#�����ϲ�λ����ͼע
legend("topleft",sampleNames(Data),col=colors,lwd=1,inset=0.05,cex=0.3)


RMA��Ԥ����normal����
library(affyPLM)
library(affy)
Data<-ReadAffy()
#sampleNames(Data)
#N=length(Data)
#��RMAԤ��������
eset.rma<-rma(Data)
#��ȡ�������ݲ����������
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="../normal.expres.txt",sep='\t',quote=F,row.names=F)

RMA��Ԥ����tumor����
library(affyPLM)
library(affy)
Data<-ReadAffy()
#sampleNames(Data)
#N=length(Data)
#��RMAԤ��������
eset.rma<-rma(Data)
#��ȡ�������ݲ����������
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="../tumor.expres.txt",sep='\t',quote=F,row.names=F)

�ϲ�N��T������
#setwd(" ")
normal_exprs<-read.table("normal.expres.txt",header=T,sep="\t")
tumor_exprs<-read.table("tumor.expres.txt",header=T,sep="\t")
#��T��N�ϲ�
probe_exprs<-merge(normal_exprs,tumor_exprs,by="probeid")
write.table(probe_exprs,file="cancer.probeid.exprs.txt",sep='\t',quote=F,row.names=F)


Probe IDת��ΪGene symbol����ƽ̨�ļ�����������
#setwd("")
#��ȡ��������ļ�
probe_exp<-read.table("cancer.probeid.exprs.txt",header=T,sep="\t",row.names=1)
#��ȡ̽���ļ�
probeid_geneid<-read.table("../GPL570-55999.txt",header=T,sep="\t")
#probeid_geneid<-read.table("../GPL6244-17930.txt",header=T,sep="\t",colClasses=c("ID"="character"))
probe_name<-rownames(probe_exp)
#probe����ƥ��
loc<-match(probeid_geneid[,1],probe_name)
#ȷ����ƥ���ϵ�probe����ֵ
probe_exp<-probe_exp[loc,]
#ÿ��probeidӦ�Ե�geneid
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
#raw_geneid<-as.numeric(as.matrix(probeid_geneid[,1]))
#�ҳ���geneid��probeid����������
index<-which(!is.na(raw_geneid))
#��ȡ��geneid��probe
geneid<-raw_geneid[index]
#�ҵ�ÿ��geneid�ı���ֵ
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
#���̽���Ӧ1������������ȡƽ��ֵ
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
#geneid��Ϊ����
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="Gastric.cancer.geneid.exprs.txt",sep='\t',quote=F,row.names=F)
#��gene id ת��Ϊgene symbol
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
#loc<-match(rownames(gene_exp_matrix),probeid_geneid[,1])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="Gastric.cancer.genesyb.exprs.txt",sep='\t',quote=F,row.names=F)


����ȱʧֵ                         
(��Ҫ��genesyb����ļ���Ҫ����)����ھӷ���KNN��k-Nearest Neighbor�������˷�����Ѱ�Һ���ȱʧֵ�Ļ���ı��������Ƶ���������ͨ����Щ����ı���ֵ�����ձ����������Լ�Ȩ�������ȱʧֵ
#���ú���impute
library(impute)
#��ȡ����ֵ
gene_exp_matrix<-read.table("Gastric.cancer.genesyb.exprs.txt",header=T,sep="\t",row.names=1)
gene_exp_matrix<-as.matrix(gene_exp_matrix)
#KNN������ȱʧֵ
imputed_gene_exp<-impute.knn(gene_exp_matrix,k=10,rowmax=0.5,colmax=0.8,maxp=3000,rng.seed=362436069)
#��������ȱʧֵ����������
GeneExp<-imputed_gene_exp$data
#д�����
genesymbol<-rownames(GeneExp)
GeneExp<-cbind(genesymbol,GeneExp)
write.table(GeneExp,file="Gastric.cancer.gene.exprs.txt",sep='\t',quote=F,row.names=F)


library(limma)
rt<-read.table("Gastric.cancer.gene.exprs.txt",header=T,sep="\t",row.names="genesymbol")
#differential
class<-c(rep("Normal",40),rep("Tumor",32))
design<-model.matrix(~factor(class))
colnames(design)<-c("Normal","Tumor")
fit<-lmFit(rt,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)
#write table
diffLab<-allDiff[with(allDiff, ((logFC>1 |logFC<(-1)) & adj.P.Val<0.05)),]
write.table(diffLab,file="diffEXp.xls",sep="\t",quote=F)
#������
diffExpLevel<-rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLevel.xls",sep="\t",quote=F)

#heatmap
hmExp=log10(diffExpLevel+0.00001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",height=8,width=8)
par(oma=c(3,3,3,5))
plot_color = ifelse(class=='Mild','grey','orange')
heatmap.2(hmMat,col='greenred',trace="none",cexCol=.8,density.info="histogram",ColSideColors = plot_color,labRow = "",key.par=list(mar=c(1,0,3,1)))
legend(y=1.1654, x=.8, xpd=TRUE,     
       legend = c('Mild','Advanced'),
       col = c('grey','orange'), 
       border = 'red',
       lty= 5,             
       lwd = 5,           
       cex= 1
)
dev.off()

#volcano
pdf(file="vol.pdf")
xMax=max(abs(allDiff$logFC))
yMax=max(-log10(allDiff$adj.P.Val))
plot(allDiff$logFC,-log10(allDiff$adj.P.Val),xlab="Log fold change",ylab="-Log10 adj.P-Val",main="Volcano",xlim=c(-xMax,xMax),ylim=c(0,yMax),pch=20,cex=0.4)
diffSub1=subset(allDiff,allDiff$adj.P.Val<0.05 & allDiff$logFC > 1)
diffSub2=subset(allDiff,allDiff$adj.P.Val<0.05 & allDiff$logFC < -1)
points(diffSub1$logFC,-log10(diffSub1$adj.P.Val),pch=20,col="red",cex=0.4)
points(diffSub2$logFC,-log10(diffSub2$adj.P.Val),pch=20,col="blue",cex=0.4)
legend('topright', xpd=TRUE,     
       legend = c('DOWN','NOT','UP'),
       col = c('blue','black','red'), 
       #box.lwd = 0,box.col = "white",bg = "white",
       pch= 16,          
       cex= 1
)
legend("bottomright", paste(nrow(diffSub1)), bty="n",text.col = 'red',text.font = 4,pt.cex = 1,cex=2)
legend("bottomleft", paste(nrow(diffSub2)), bty="n",text.col = 'blue',text.font = 4,pt.cex = 1,cex=2)
#abline(v=0,lty=2,lwd=3)
dev.off()

#DEG 
library(clusterProfiler)
rt=read.table("X.txt",sep="\t",head=T,check.names=F)
geneFC=rt$logFC
gene<-rt$ENTREZ_GENE_ID
names(geneFC)=gene
#kegg
kk<-enrichKEGG(gene=gene,organism="human",qvalueCutoff=0.05)#,readable=TRUE)
write.table(summary(kk),file="KEGG.xls",sep="\t",quote=F,row.names=F)
pdf(file="KEGG.barplot.pdf")
barplot(kk,drop=TRUE,showCategory=12)
dev.off()
pdf(file="KEGG.cnetplot.pdf")
cnetplot(kk,categorySize="geneNum",foldChange=geneFC)
dev.off()
library("pathview")
keggxls=read.table("KEGG.xls",sep="\t",header=T)
for(i in keggxls$ID){
  pv.out<-pathview(gene.data=geneFC,pathway.id=i,species="hsa",out.suffix="pathview")}

library("ggplot2")
#1800*600
enrich <- read.delim("C:/J's/liver/chart.txt")
#enrich2 = subset(enrich,enrich$FDR<=0.05)
enrich$FDR=-log10(enrich$FDR)
plot<-ggplot()+ylab('-Log10 (FDR)')+xlab('')
plot<-plot+theme_classic()+theme(panel.background=element_rect(fill='transparent',color='black'),title=element_text(size=20,vjust=1.4,hjust=0.5,face='bold',color='black'),legend.text=element_text(size=20,face='bold'),legend.key.height = unit(1, "cm"),axis.text.x=element_text(size=15,face='bold',color='black'),axis.title.x=element_text(hjust=0.5,face='bold'),axis.text.y=element_text(size=20,face='bold',color='black'))
plot<-plot+geom_bar(data=enrich,aes(x=reorder(Term, FDR), y=FDR, fill=Category),stat="identity")+coord_flip()
plot
