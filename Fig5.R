library(vegan)
library(ggplot2)
library(psych)
library(igraph)
library(tibble)
library(readxl)

####Fig4a####
data <- KO <- read_excel("Supplcompst.xlsx", sheet = "KO_table") %>%
  column_to_rownames(var = colnames(.)[1])
group <- read_excel("Supplenv.xlsx", sheet = "group")
ko <- t(data)
ko1 <- ko[,colSums(ko) > 0 & specnumber(t(ko)) > 20 ] 


spman.r0  <-  corr.test(ko1, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)


Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

r.cutoff = 0.6
p.cutoff = 0.05
df <- data.frame(Cor.df,  P.df)
da.tmp<-df.sig<- df[ df$Cor > r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)

data <- read_excel("Supplcompst.xlsx", sheet = "KO_table")
tax <- read_excel("Supplcompst.xlsx", sheet = "ko_tax")
tax <- tax[tax$KO %in% data$KEGG_ko, ] 
ID.ko <- tax
IDsub1<-ID.ko[ID.ko$KO %in% V1$v1, ]; IDsub2<-ID.ko[ID.ko$KO %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "KO", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "KO", all.x = T); M2<-M2[order(M2$id), ]
da<-data.frame(da.tmp, M1, M2)

g <- graph_from_data_frame(da, directed=FALSE)  # 创建无向网络
length(V(g)) 
length(E(g))
fun.fc <- cluster_fast_greedy(g)
print(modularity(fun.fc))

ko.modu.memb <- data.frame(c(membership(fun.fc) ))
ko.modu.memb <- data.frame("KOID"= row.names(ko.modu.memb), 
                           "Modular"= ko.modu.memb$c.membership.fun.fc..)
ko.sub <- ID.ko [ID.ko$KO %in% ko.modu.memb$KOID,]
ko.modu.memb.ID <- merge(ko.modu.memb, ko.sub, by.x= "KOID", by.y = "KO",all = T)


print(sizes(fun.fc))###

fun.comps <- membership(fun.fc)
colbar <- c("#0000FF","#FF3030")
V(g)$color <- colbar[fun.comps]
set.seed(123)
pdf("konetwk.pdf", width = 6, height = 6)
plot(g, layout = layout_with_kk, edge.width=0.07,edge.color="grey", vertex.frame.color=NA,vertex.label=NA,edge.lty=1,
     edge.curved=T,vertex.size=1,margin=c(0, 0,0,0))
dev.off()

####Fig4b####
library(vegan)
library(psych)
library(igraph)
library(ggplot2)
library(psych)
library(igraph)
library(tidyr)
library(nlme) 
library(MuMIn)

data <- KO <- read_excel("Supplcompst.xlsx", sheet = "KO_table") %>%
  column_to_rownames(var = colnames(.)[1])
env <- read_excel("Supplenv.xlsx", sheet = "group")
ko <- t(data)
ko1 <- ko[,colSums(ko) > 0 & specnumber(t(ko)) > 20 ]


ko.modu.memb.ID.mod1 <- ko.modu.memb.ID[ko.modu.memb.ID$Modular==1,]
ko.mgm.mod1 <- ko1[,colnames(ko1) %in% ko.modu.memb.ID.mod1$KOID]
ko.mgm.mod1 <- ko.mgm.mod1[order(rownames(ko.mgm.mod1)),]
env <- env[order(env$sampleid),]
row.names(ko.mgm.mod1) == env$sampleid

ko.mgm.mod1.pH<- data.frame(ph=env[,"pH"], ko.mgm.mod1)

ko.mgm.mod1.pH.long <- gather(ko.mgm.mod1.pH, KOID, abu, K00001:K22504, factor_key=TRUE)
ko.mgm.mod1.pH.mean <- ko.mgm.mod1.pH
env$mod1<-rowSums(ko.mgm.mod1) / ncol(ko.mgm.mod1)
#write.csv(ko.modu.memb.ID.mod1, "ko_M1.csv")

ggplot() + 
  geom_smooth(data=ko.mgm.mod1.pH.long, aes(x= ph, y= log(abu+1), group = KOID), col= "grey", method = "lm" , linewidth = 0.05,  fill=NA)+
  geom_smooth(data=env, aes(x= pH, y= log(mod1+1) ), col= "#0000FF", method = "lm", linewidth = 5)+
  labs(x = "Soil pH",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"))

lme0<-lme(abu~ph,random=~1|KOID,data=ko.mgm.mod1.pH.long, control=lmeControl(opt = "optim"))
summary(lme0)
anova(lme0)
r.squaredGLMM(lme0)
####Fig4c####
ko.modu.memb.ID.mod2 <- ko.modu.memb.ID[ko.modu.memb.ID$Modular==2,]
ko.mgm.mod2 <- ko1[,colnames(ko1) %in% ko.modu.memb.ID.mod2$KOID]
ko.mgm.mod2 <- ko.mgm.mod2[order(rownames(ko.mgm.mod2)),]

row.names(ko.mgm.mod2) == env$sampleid

ko.mgm.mod2.pH<- data.frame(ph=env[,"pH"], ko.mgm.mod2)

ko.mgm.mod2.pH.long <- gather(ko.mgm.mod2.pH, KOID, abu, K00002:K22508, factor_key=TRUE)
ko.mgm.mod2.pH.mean <- ko.mgm.mod2.pH
env$mod2<-rowSums(ko.mgm.mod2) / ncol(ko.mgm.mod2)
#write.csv(ko.modu.memb.ID.mod2,"ko_M2.csv")

ggplot() + 
  geom_smooth(data=ko.mgm.mod2.pH.long, aes(x= ph, y= log(abu+1), group = KOID), col= "grey", method = "lm" , linewidth = 0.05,  fill=NA)+
  geom_smooth(data=env, aes(x= pH, y= log(mod2+1) ), col= "#FF3030", method = "lm", linewidth = 5,  fill=NA)+
  labs(x = "Soil pH",y = "Log (abundance + 1)")+
  ylim(-1,10)+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"))

lme0<-lme(abu~ph,random=~1|KOID,data=ko.mgm.mod2.pH.long, control=lmeControl(opt = "optim"))
summary(lme0)
anova(lme0)
r.squaredGLMM(lme0)

####Fig4d####
library(ggrepel)
M1.ko.level <- ko.modu.memb.ID.mod1
M2.ko.level <- ko.modu.memb.ID.mod2
data <- read_excel("Supplcompst.xlsx", sheet = "eggnog.KO.TPM_ant")

M1.ko.level <- merge(M1.ko.level, data[,c(1,43:45)], by.x = "KOID", by.y = "KEGG_ko")
M2.ko.level <- merge(M2.ko.level, data[,c(1,43:45)], by.x = "KOID", by.y = "KEGG_ko")

# 计算每个 level2 有多少不同的 KOID
M1.ko.level2 <- aggregate(KOID ~ level2, data = M1.ko.level, FUN = function(x) length(unique(x)))
colnames(M1.ko.level2)[2] <- "M1_count"



M2.ko.level2 <- aggregate(KOID ~ level2, data = M2.ko.level, FUN = function(x) length(unique(x)))
colnames(M2.ko.level2)[2] <- "M2_count"


M1M2_level2<-merge(M1.ko.level2,M2.ko.level2,by="level2")


M1M2_level2<-tibble::column_to_rownames(M1M2_level2,var = "level2")
group_list <- factor(c(rep("M1_count",1),rep("M2_count",1)))
exprSet <- DGEList(counts = M1M2_level2, group = group_list)
exprSet
bcv = 0.01
et <- exactTest(exprSet, dispersion=bcv^2)

gene<-data.frame(topTags(et, n = nrow(exprSet$counts)))

gene[which(gene$FDR < 0.05 & gene$logFC <= -0.5),'sig'] <- 'Down'
gene[which(gene$FDR < 0.05 & gene$logFC >=0.5),'sig'] <- 'Up'
gene[which(gene$FDR >= 0.05 | abs(gene$logFC) < 0.5),'sig'] <- 'None'
sum(gene$sig=="Down")
sum(gene$sig=="Up")
sum(gene$sig=="None")

#write.csv(gene,"pathwy2_result.csv")
dataset <- read.csv("pathwy2_result.csv", header = T)


ggplot(dataset, aes(x = logFC,y = -log10(PValue),colour=sig)) +
  geom_point(alpha=0.5, size=6) +
  xlim(-6,6)+ ylim (0,10)+
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.02),lty=4,col="grey",lwd=0.5) +
  scale_y_continuous(breaks=seq(0,10,by=2))+
  geom_text_repel(data=subset(dataset, FDR < 0.05 & sig != "None"),vjust="inward",hjust="inward",
                  aes(label=X, color = sig), angle = 0, size=5, fontface= "bold")+
  labs(x="log2(fold change)", y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"))

####Fig4e####
rm(list = ls())

data <- read_excel("Supplcompst.xlsx", sheet = "Feng_kotable") %>%
  column_to_rownames(var = colnames(.)[1])
ko <- t(data)
ko1 <- ko[,colSums(ko) > 0 & specnumber(t(ko)) > 18 ] 
spman.r0  <-  corr.test(ko1, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)


Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

r.cutoff = 0.6
p.cutoff = 0.05
df <- data.frame(Cor.df,  P.df)
da.tmp<-df.sig<- df[ df$Cor > r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)

data <- read_excel("Supplcompst.xlsx", sheet = "Feng_kotable")
tax <- read_excel("Supplcompst.xlsx", sheet = "ko_tax")
tax <- tax[tax$KO %in% data$KO, ] 
ID.ko <- tax
IDsub1<-ID.ko[ID.ko$KO %in% V1$v1, ]; IDsub2<-ID.ko[ID.ko$KO %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "KO", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "KO", all.x = T); M2<-M2[order(M2$id), ]
da<-data.frame(da.tmp, M1, M2)


g <- graph_from_data_frame(da, directed=FALSE)
length(V(g)) 
length(E(g))
fun.fc <- cluster_fast_greedy(g)
print(modularity(fun.fc))

ko.modu.memb <- data.frame(c(membership(fun.fc) ))
ko.modu.memb <- data.frame("KOID"= row.names(ko.modu.memb), 
                           "Modular"= ko.modu.memb$c.membership.fun.fc..)
ko.sub <- ID.ko [ID.ko$KO %in% ko.modu.memb$KOID,]
ko.modu.memb.ID <- merge(ko.modu.memb, ko.sub, by.x= "KOID", by.y = "KO",all = T)


print(sizes(fun.fc))
fun.comps <- membership(fun.fc)
colbar <- c("#FF3030","#0000FF")
V(g)$color <- colbar[fun.comps]
set.seed(123)
pdf("konetwk.pdf", width = 6, height = 6)
plot(g, layout = layout_with_kk, edge.width=0.07,edge.color="grey", vertex.frame.color=NA,vertex.label=NA,edge.lty=1,
     edge.curved=T,vertex.size=1,margin=c(0, 0,0,0))
dev.off()




####Fig4f####
data <- read_excel("Supplcompst.xlsx", sheet = "Feng_kotable") %>%
  column_to_rownames(var = colnames(.)[1])
ko <- t(data)
ko1 <- ko[,colSums(ko) > 0 & specnumber(t(ko)) > 18 ] 
env <- read_excel("Supplenv.xlsx", sheet = "Feng_group")

ko.modu.memb.ID.mod2 <- ko.modu.memb.ID[ko.modu.memb.ID$Modular==2,]
ko.mgm.mod2 <- ko1[,colnames(ko1) %in% ko.modu.memb.ID.mod2$KOID]
ko.mgm.mod2 <- ko.mgm.mod2[order(rownames(ko.mgm.mod2)),]
env <- env[order(env$site),]

row.names(ko.mgm.mod2) == env$site
env$log2_EC <- log2(env$Ec)

ko.mgm.mod2.EC<- data.frame(EC=env[,"log2_EC"], ko.mgm.mod2)

ko.mgm.mod2.EC.long <- gather(ko.mgm.mod2.EC, KOID, abu, K12132:K21479, factor_key=TRUE)
ko.mgm.mod2.EC.mean <- ko.mgm.mod2.EC
env$mod2<-rowSums(ko.mgm.mod2) / ncol(ko.mgm.mod2)
#write.csv(ko.modu.memb.ID.mod2, "ko_M2.csv")

ggplot() + 
  geom_smooth(data=ko.mgm.mod2.EC.long, aes(x= EC, y= log(abu+1), group = KOID), col= "grey", method = "lm" , linewidth = 0.05,  fill=NA)+
  geom_smooth(data=env, aes(x= log2_EC, y= log(mod2+1) ), col= "#0000FF", method = "lm", linewidth = 5)+
  labs(x = "log2(EC)",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"))
lme0<-lme(abu~EC,random=~1|KOID,data=ko.mgm.mod2.EC.long, control=lmeControl(opt = "optim"))
summary(lme0)
anova(lme0)
r.squaredGLMM(lme0)

####Fig4g####

ko.modu.memb.ID.mod1 <- ko.modu.memb.ID[ko.modu.memb.ID$Modular==1,]
ko.mgm.mod1 <- ko1[,colnames(ko1) %in% ko.modu.memb.ID.mod1$KOID]
ko.mgm.mod1 <- ko.mgm.mod1[order(rownames(ko.mgm.mod1)),]
row.names(ko.mgm.mod1) == env$site


ko.mgm.mod1.EC<- data.frame(EC=env[,"log2_EC"], ko.mgm.mod1)

ko.mgm.mod1.EC.long <- gather(ko.mgm.mod1.EC, KOID, abu, K01768:K19820, factor_key=TRUE)
ko.mgm.mod1.EC.mean <- ko.mgm.mod1.EC
env$mod1<-rowSums(ko.mgm.mod1) / ncol(ko.mgm.mod1)
#write.csv(ko.modu.memb.ID.mod1, "ko_M1.csv")

ggplot() + 
  geom_smooth(data=ko.mgm.mod1.EC.long, aes(x= EC, y= log(abu+1), group = KOID), col= "grey", method = "lm" , linewidth = 0.05,  fill=NA)+
  geom_smooth(data=env, aes(x= log2_EC, y= log(mod1+1) ), col= "#FF3030" , method = "lm", linewidth = 5)+
  labs(x = "log2(EC)",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"))

lme0<-lme(abu~EC,random=~1|KOID,data=ko.mgm.mod1.EC.long, control=lmeControl(opt = "optim"))
summary(lme0)
anova(lme0)
r.squaredGLMM(lme0)

####Fig4h####
M1.ko.level <- ko.modu.memb.ID.mod2
M2.ko.level <- ko.modu.memb.ID.mod1
data <- read_excel("Supplcompst.xlsx", sheet = "Feng_KO_annotation")

M1.ko.level <- merge(M1.ko.level, data[,c(1,39:41)], by.x = "KOID", by.y = "KO")
M2.ko.level <- merge(M2.ko.level, data[,c(1,39:41)], by.x = "KOID", by.y = "KO")


M1.ko.level2 <- aggregate(KOID ~ level2, data = M1.ko.level, FUN = function(x) length(unique(x)))
colnames(M1.ko.level2)[2] <- "M1_count"



M2.ko.level2 <- aggregate(KOID ~ level2, data = M2.ko.level, FUN = function(x) length(unique(x)))
colnames(M2.ko.level2)[2] <- "M2_count"


M1M2_level2<-merge(M1.ko.level2,M2.ko.level2,by="level2")


M1M2_level2<-tibble::column_to_rownames(M1M2_level2,var = "level2")
group_list <- factor(c(rep("M1_count",1),rep("M2_count",1)))
exprSet <- DGEList(counts = M1M2_level2, group = group_list)
exprSet
bcv = 0.01
et <- exactTest(exprSet, dispersion=bcv^2)

gene<-data.frame(topTags(et, n = nrow(exprSet$counts)))

gene[which(gene$FDR < 0.05 & gene$logFC <= -0.5),'sig'] <- 'Down'
gene[which(gene$FDR < 0.05 & gene$logFC >=0.5),'sig'] <- 'Up'
gene[which(gene$FDR >= 0.05 | abs(gene$logFC) < 0.5),'sig'] <- 'None'
sum(gene$sig=="Down")
sum(gene$sig=="Up")
sum(gene$sig=="None")

#write.csv(gene,"Feng_pathwy2_result.csv")
dataset <- read.csv("Feng_pathwy2_result.csv", header = T)

ggplot(dataset, aes(x = logFC,y = -log10(PValue),colour=sig)) +
  geom_point(alpha=0.5, size=6) +
  xlim(-6,6)+ ylim (0,10)+
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.02),lty=4,col="grey",lwd=0.5) +
  scale_y_continuous(breaks=seq(0,10,by=2))+
  geom_text_repel(data=subset(dataset, FDR < 0.05 & sig != "None"),vjust="inward",hjust="inward",
                  aes(label=X, color = sig), angle = 0, size=5, fontface= "bold",max.overlaps = Inf,
                  force = 3,force_pull = 2,
                  segment.size = 0.5)+
  labs(x="log2(fold change)", y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"))

