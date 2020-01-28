#setwd('/Users/panzhang/Desktop/RNA-seq_pig/Pig_Hypo116')
library("edgeR")
# read gene raw count data, each row represents a gene, each column is a sample
data <-  read.delim("gene_counts.txt", row.names=1, stringsAsFactors=FALSE)
dim(data)
colnames(data)
targets <- data.frame(Lane = c(1:112),Sex = c(rep("F",31),rep("M",31),rep("F",21),rep("M",29)),
                      PigletTrt = c(rep("con",9),rep("fast",10),rep("poly",12),rep("con",12),rep("fast",9),rep("poly",10),rep("con",6),rep("fast",6),rep("poly",9),rep("con",11),rep("fast",9),rep("poly",9)),
                      GiltTrt = c(rep("Control", 62),rep("Prrs", 50)))
Group <- factor(paste(targets$Sex, targets$PigletTrt, targets$GiltTrt, sep="."))
cbind(targets,Group=Group)
y <- DGEList(counts=data[,1:112],group=Group)
colnames(y) <- colnames(data)
dim(y)
keep <- rowSums(cpm(y)>1) >= 9
y <- y[keep,]
dim(y)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y
plotMDS(y,col=c(rep("black",9),rep("yellow",10), rep("blue",12),
                rep("red",12),rep("purple",9),rep("orange",10),
                rep("green",6),rep("dark blue",6),rep("dark red",9),
                rep("cyan",11),rep("grey",9),rep("chocolate",9)))


#only genewise dispersion
#PigletTrt
design <- model.matrix(~0+GiltTrt*Sex*PigletTrt,data=targets)
colnames(design)
x <- estimateGLMTrendedDisp(y, design)
x <- estimateGLMTagwiseDisp(x, design)
fit <- glmFit(x,design)

#GiltTrt main effect
qlf <- glmLRT(fit,contrast =c(1,-1,0,0,0,0,0,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="GLMTrend_GiltTrt_main.txt",sep = "\t",quote = FALSE)

#Sex main effect
qlf <- glmLRT(fit,coef = 3)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="GLMTrend_Sex_main.txt",sep = "\t",quote = FALSE)

#PigletTrt main effect
qlf <- glmLRT(fit,coef = 4:5)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="GLMTrend_PigletTrt_main.txt",sep = "\t",quote = FALSE)

#GiltTrt-by-Sex interaction
qlf <- glmLRT(fit,coef = 6)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="GLMTrend_GiltTrt-by-Sex.txt",sep = "\t",quote = FALSE)

#PigletTrt-by-sex interaction
qlf <- glmLRT(fit,coef = 9:10)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="GLMTrend_PigletTrt-by-Sex.txt",sep = "\t",quote = FALSE)

#PigletTrt-by-GiltTrt interaction
qlf <- glmLRT(fit,coef = 7:8)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="GLMTrend_PigletTrt-by-GiltTrt.txt",sep = "\t",quote = FALSE)


#pairwise comparison
design <- model.matrix(~0+Group,data=targets)
colnames(design)
x <- estimateGLMTrendedDisp(y, design)
#x <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(x,design)


qlf <- glmLRT(fit,contrast =c(1,0,0,0,0,0,-1,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="F.con.C_vs_M.con.C.txt",sep = "\t",quote = FALSE)

qlf <- glmLRT(fit,contrast =c(0,1,0,0,0,0,-1,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="F.con.P_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,1,0,0,0,-1,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="F.fast.C_vs_M.con.C.txt",sep = "\t",quote = FALSE)

qlf <- glmLRT(fit,contrast =c(0,0,0,1,0,0,-1,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="F.fast.P_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,1,0,-1,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="F.poly.C_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,0,1,-1,0,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="F.poly.P_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,0,0,-1,1,0,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="M.con.P_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,0,0,-1,0,1,0,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="M.fast.C_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,0,0,-1,0,0,1,0,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="M.fast.P_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,0,0,-1,0,0,0,1,0)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="M.poly.C_vs_M.con.C.txt",sep = "\t",quote = FALSE)


qlf <- glmLRT(fit,contrast =c(0,0,0,0,0,0,-1,0,0,0,0,1)) #contrast = c(1,-1,0,0)
top <- topTags(qlf, n=Inf)
summary(de <- decideTestsDGE(qlf))
write.table(top,file="M.poly.P_vs_M.con.C.txt",sep = "\t",quote = FALSE)
