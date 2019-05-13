---
title: "Telomere candidate gene search"
author: "JT Lovell, EV Shakirov"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
    toc: true
    theme: united
    fig_width: 8
    fig_height: 8
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### 1. Set up the environment
#### Load required packages, install updates from github
```{r loadPackages, warnings = FALSE, results = 'hide'}
rm(list = ls())

pkg<-c("qtlTools","limmaDE2","knitr","data.table","plyr",
       "ggplot2","reshape2","mpMap")
suppressMessages(sapply(pkg, require, character.only = TRUE))
```

### 2. Load data
#### The Arabidopsis magic population data
```{r loadMAGIC,results = 'hide'}
setwd("/Users/John/Desktop/magic")
load('/Users/John/Desktop/magic/ArabidopsisExample.RData')
tel<-read.delim("~/Desktop/magic/telomere.all.txt")

### Add the phenotype data to the MAGIC R object
phe<-armpc$pheno
phe$SUBJECT.NAME<-rownames(phe)
phe<-merge(phe, tel, by = "SUBJECT.NAME", all.x = T)
rownames(phe)<-phe$SUBJECT.NAME
phe$SUBJECT.NAME<-NULL
phe<-phe[match(rownames(armpc$pheno),rownames(phe)),]
phe$mean.trf[is.na(phe$mean.trf)]<-mean(phe$mean.trf)
armpc$pheno<-phe
```

#### The Arabidopsis tair 10 gene annotation
```{r loadTAIR, results = 'hide'}
gff<-read.delim("/Users/John/Desktop/juengerLabProjects/telomere_data/TAIR10_GFF3_genes.gff", header=F) 
colnames(gff)<-c("chr","tair.version","type","start","end","unk","strand","unk2","info")
gff<-gff[gff$type == "gene",]
gff<-gff[gff$chr=="Chr5",]

gff<-gff[!duplicated(gff$info),]
gff$geneID<-gsub("Parent=","",gff$info)
gff$geneID<-gsub("ID=","",gff$geneID)
gff$geneID<-substring(gff$geneID,1,9)
gff$info<-as.character(gff$info)
```

#### The parental alleles of the 19 genomes
```{r loadAlleles, results = 'hide'}
snps<-fread("/Users/John/Desktop/magic/chr5.alleles.txt") 
snps<-data.frame(snps)
```

#### The parental expression counts of the 19 genomes
```{r loadCounts, results = 'hide'}
counts<-read.delim("/Users/John/Desktop/juengerLabProjects/telomere_data/counts19genomes.tab", header=T)
colnames(counts)<-gsub("X.fml.ag.raetsch.nobackup.projects.sequencing_runs.A_thaliana_magic.results.alignments.new2_splice.","",colnames(counts))
colnames(counts)<-gsub("R1.4.nodup.","",colnames(counts))
counts<-counts[,c(1, grep("nodup.RPKM",colnames(counts)))]
colnames(counts)<-gsub(".R2.4.nodup.RPKM","",colnames(counts))
```

#### The southern blotting candidate gene results
```{r}
tdna<-read.csv("/Users/John/Desktop/juengerLabProjects/telomere_figures/southernResults.csv",
               stringsAsFactors=F)
tdna$Gene<-toupper(tdna$Gene)
```


### 3. Run mpMap composite interval mapping
```{r runMpMap, results = 'hide'}
s1 <- mpIM(object=armpc, ncov=1, responsename="mean.trf", dwindow=50, window = 50)
```

```{r}
tmp<-armpc$prob$Chr5
tmp<-ifelse(tmp<0.5, FALSE, TRUE)
tmp[tmp<0.5]<-FALSE
tmp[tmp>=0.5]<-TRUE
tmp<-data.frame(tmp)

wh <- apply(tmp[,grep("MN5_22714627", colnames(tmp))],1, which.max)
als <- as.numeric(sapply(colnames(tmp[,grep("MN5_22714627", colnames(tmp))])[wh], function(x) strsplit(x,".", fixed = T)[[1]][4]))
aln <- rownames(armpc$founders)
names(als)<-names(wh)
alo <- aln[als]
names(alo)<-names(wh)

al.dt <- data.table(id = names(alo), haplotype = alo)
p.dt <- data.table(id = tel$SUBJECT.NAME, phe = tel$mean.trf)
p <- merge(al.dt, p.dt, by = "id")
p$is.sf2 <- p$haplotype == "Sf"

library(lme4)
library(MuMIn)

#Fit Model
m <- lmer(phe ~ 1 + (1|haplotype), data = p)

#Determine R2:
r.squaredGLMM(m) 

m
mnames<-names(armpc$map$Chr5) 
out<-lapply(mnames, function(x){
  dat<-tmp[,grep(x,colnames(tmp))]
  rownames(armpc$founders)[as.numeric(unlist(apply(dat, 1, which)))]
})
out<-do.call(cbind, out)
out<-data.frame(out)
colnames(out)<-mnames
rownames(out)<-rownames(tmp)
pdat<-data.frame(armpc$pheno,allele = out[,"MN5_22714627"])
```


#### Print results
```{r printSummary}
kable(summary(s1))
```

#### Plot the results
```{r plotS1}
par(mfrow=c(2,1))
plot(s1)
plot(s1,chr = 5, lodsupport=5)
```

### 4. Analyze the effect of the focal QTL
#### Get % variance explained by each QTL
```{r}
dat<-fillmiss(armpc)
map<-dat$map
grid.mar<-unlist(lapply(map, function(x) {
  grid<-seq(from = 0, to = max(x), by = 5)
  out<-sapply(grid, function(y)  names(x)[which.min(abs(x-y))])
  return(out)
  }))
  
tomat<-apply(dat$finals, 2, function(x) ifelse(x == min(x), -1, ifelse(x == max(x), 1,0)))
amat<-A.mat(tomat[,grid.mar])
y = dat$pheno$mean.trf  


tmp<-do.call(cbind, armpc$prob)
tmp<-ifelse(tmp<0.5, FALSE, TRUE)
tmp[tmp<0.5]<-FALSE
tmp[tmp>=0.5]<-TRUE
tmp<-data.frame(tmp)
qtlposs<-unlist(attr(res$qtl,"index"))
qtlmars<-names(unlist(map))[qtlposs]
qtlmars<-sapply(qtlmars, function(y) strsplit(y, ".", fixed = T)[[1]][2])
out<-lapply(qtlmars, function(x){
  dat<-as.matrix(tmp[,grep(x,colnames(tmp))])
  line = rownames(dat)
  wh = apply(dat, 1, which)
  wh[sapply(wh,length)==0]<-NA
  wh<-as.numeric(unlist(wh))
  rownames(armpc$founders)[wh]
})
out<-do.call(cbind, out)
out<-data.frame(out)
colnames(out)<-qtlmars
rownames(out)<-rownames(tmp)

qtlMar<-out
identical(rownames(tomat), rownames(qtlMar))

y = dat$pheno$mean.trf
adat<-cbind(y, qtlMar)
fit<-aov(y~MN5_22714627, data = adat)
af <- anova(fit)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
 
aovnull<-aov(y~1, data = adat)
fixed.matrix<-model.matrix(form, data = qtlMar2)
test.q5<-rrBLUP::mixed.solve(y=y, K = amat,
                              X = fixed.matrix)
print(test.q5$Vu); print(test.q5$Ve)
print(test.q5$Vu/(test.q5$Vu+test.q5$Ve))

form<-as.formula("~ MASC04211")
fixed.matrix<-model.matrix(form, data = qtlMar2)
test.null<-rrBLUP::mixed.solve(y=y, K = amat)
print(test.null$Vu); print(test.null$Ve)
print(test.null$Vu/(test.null$Vu+test.null$Ve))
```


#### Pull out means and ses
```{r getMeans}
res<-s1$QTLresults
est<-do.call(cbind,res$fndrfx)
se<-do.call(cbind,res$se)
qtl.est<-est[,attr(res$qtl,"index")$Chr5]
founder.id<-rownames(armpc$founders)
qtl.se<-se[,attr(res$qtl,"index")$Chr5]
qtl.estimate<-data.frame(founder.id,qtl.est,qtl.se, stringsAsFactors=F)
```

#### Plot means and SE
```{r plotEffect}
ggplot(qtl.estimate, aes(x = founder.id, y = qtl.est, fill = qtl.est))+
  geom_bar(stat="identity", width = .5)+ 
  scale_fill_gradient(low = "grey", high = "red", guide = F)+
  theme_jtlbar()+
  geom_errorbar(aes(ymin = qtl.est-qtl.se, ymax = qtl.est+qtl.se), width = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  labs(x = "MAGIC Founder Parent", y = "QTL effect estimate (+/- SE)",
       title = "Chr 5 QTL Effect")
```

#### Pull out all genes in the interval
```{r getInterval}
wh<-which.max(-log10(summary(s1)$pvalue))
for(i in c(2.5,5,7.5,10)){
 interval<-as.numeric(supportinterval(s1, lodsupport=i)$support[,wh])
 print(names(s1$map$Chr5[s1$map$Chr5 %in% interval]))
}

intervals<-data.frame(level = c(2.5,7.5),
                      lower = c(22415000,22414731),
                      upper = c(23115566,23115566), stringsAsFactors=F)
intervals$width = with(intervals, upper-lower)
interval.bp<-c(22397775, 23246248) 
```

### 5. Get statistics for expression and mutations of all genes in the interval
```{r compareFunction, echo = F}
# this function takes the gff, snps from parents and expression counts, plus a window
# it parses the files based on where the genes are and calculates various statistics
#########################################################
#########################################################
compareSF2<-function(gff, snps, counts, lowCI, highCI){
  #prepare annotation - subset to region
  gff_drop<-gff[gff$start>lowCI & gff$end<highCI,]
  cat("comparing SF-2 to the other 19 ecotypes for",nrow(gff_drop),"genes\n")
  genos<-colnames(snps)[6:24]
  narrow<-snps[snps$bp>lowCI & snps$bp<highCI,]
  gff_drop$info<-as.character(gff_drop$info)
  cat("calculating the number of SF-2 private alleles in each gene\n")
  out<-apply(narrow[,genos],1,function(x) {
    sfa<-x[14]
    if(sfa==0){
      return(NA)
    }else{
      s<-sum(x[-14]==sfa)
      return(s)
    }
  })
  narrow$propsf2<-out
  narrow$onlysf<-ifelse(narrow$propsf2==0 & !is.na(narrow$propsf2),1,0)
  #count the number of private alleles to sf2
  n.private_sf2<-sapply(unique(gff_drop$info), function(i){
    wind<-c(gff_drop$start[gff_drop$info==i], gff_drop$end[gff_drop$info==i])
    dat<-narrow[narrow$bp>wind[1] & narrow$bp<wind[2],]
    
    sum(dat$onlysf)
  })
  gff_drop$n.privateSf2<-n.private_sf2
  
  #pull out the genes in the region
  genes_drop<-gff_drop
  narrowGenes<-genes_drop$geneID
  genes_drop$gene<-genes_drop$geneID
  
  #subset the expression counts data by the genes
  counts_5<-counts[counts$gene %in% narrowGenes,]
  
  #calculate statistics about the distribution of counts
  cat("calculating expression differences between SF-2 and the other ecotypes\n")
  exp.test<-apply(counts_5[,-1], 1, function(x) {
    sf2<-x[14]
    other<-x[-14]
    if(sd(x)<0.000001){
      out.t<-NA
    }else{
      out.t<-t.test(other, mu=sf2)$p.value # one sample t-test where there is some variance
    }
    out.percentile<-sum(other<sf2) #n of lines with expression less than sf2
    c(out.t,out.percentile)
  })
  exp.test<-data.frame(t(exp.test))
  colnames(exp.test)<-c("t.test_pvalue","n_ecotypes_belowSF2")
  exp.test$gene<-counts_5$gene
  exp.test$t.test_qvalue<-p.adjust(exp.test$t.test_pvalue, method="bonferroni")
  exp.test<-exp.test[complete.cases(exp.test),]
  out<-merge(genes_drop, exp.test, by="gene")
  return(out)
}
#########################################################
#########################################################
```

#### Calculate stats for all genes in the window, comparing SF-2 to all other parents
```{r}
dropstats<-compareSF2(gff=gff, snps=snps, counts=counts, 
                      lowCI=min(interval.bp)-100, 
                      highCI=max(interval.bp)+100)
dropstats$priv.cat<-with(dropstats,ifelse(n.privateSf2==0,"same","sf2private"))
```

```{r}
ggplot(dropstats, aes(x=n_ecotypes_belowSF2, y=(-log10(t.test_pvalue)), col=priv.cat))+
  geom_point()+
  scale_x_continuous("n ecotype with expression less than SF")+
  scale_y_continuous("-log10 Pvalue of 1-sample t-test for mean 18 ecotypes != SF-2")+
  geom_abline(intercept=-log10(.05/nrow(dropstats)), slope=0, col="red", lty=2)+
  ggtitle(paste("drop7.5 LOD window contains",nrow(dropstats),"genes"))+
  theme_bw()
candidateGenes<-dropstats$gene[dropstats$n.privateSf2>0 & dropstats$t.test_qvalue<0.05]
print(candidateGenes)
```

### 6. Examine genes in the interval
```{r}

boxplot(tdna$mean.TRF ~ tdna$Gene, xaxt = "n", range=0, ylab = "% Increase from Wt (Col-0)")


wt<-tdna[tdna$Gene == "WT-COL-0",]
tdna<-tdna[tdna$Gene != "WT-COL-0",]
mean.wt<-mean(wt$mean.TRF)
se.wt<-sd((wt$mean.TRF - mean.wt)/mean.wt)/sqrt(length(wt$mean.TRF))
tdna$prop.wt<-with(tdna, (mean.TRF-mean.wt)/mean.wt)
tdna.mean.se<-ddply(tdna, .(Gene), summarize,
                    mean = mean(prop.wt),
                    sd = sd(prop.wt),
                    se = sd(prop.wt)/sqrt(length(prop.wt)))

tdna2<-merge(tdna.mean.se, with(gff[gff$geneID %in% tdna.mean.se$Gene,], 
                      data.frame(Gene = geneID, pos = start, stringsAsFactors=F)),
                      by = "Gene")
tdna2$order<-order(tdna2$pos)/nrow(tdna2)

p.min<-min(dropstats$start)
p.max<-max(dropstats$start)
d.diff<-p.max-p.min
tdna2$x<-(tdna2$order*d.diff)+p.min


abline(h = -se.wt, lwd = 3, col = rgb(1,0,0,.5))
with(tdna2, plot(x,mean))
with(tdna2, segments(x0 = x, x1=x, y1 = mean-sd, y0 = mean+sd))
tdna2$mean<-tdna2$mean*100
tdna2$sd<-tdna2$sd*100
tdna2$se<-tdna2$se*100

tdna2$mean<-tdna2$mean*100
tdna2$sd<-tdna2$sd*100
tdna$prop.wt<-tdna$prop.wt*100
se.wt<-se.wt*100
```


```{r}
cand = unique(tdna$Gene)
```


```{r}
par(mar = c(0, 12, 0, 2) + 0.1, mfrow = c(2,1))
boxplot(tdna$prop.wt ~ tdna$Gene, xaxt = "n", range=0, ylab = "% Increase from Wt (Col-0)")
rect(xleft = -1, xright = 22, ytop = se.wt, ybottom = -se.wt, col = rgb(1,0,0,.5),
     border = NA)
abline(h = 0, lwd = 3, col = rgb(1,0,0,1))

with(dropstats, plot(c(min(start), max(start)),c(0,0),type = "l", ylim = c(-.2,1), axes=F,
                     xlab = "Chr 5 physical position (Mbp)",
                     ylab = NA))
with(dropstats, segments(x0 = start,x1=start,y0 = -0.05, y1 = 0.05))
with(dropstats[dropstats$gene %in% unique(tdna$Gene),], 
     points(x = start, y = rep(.6, length(start)), pch = 8, col = "green"))

with(dropstats[dropstats$n.privateSf2>0,], 
     points(x = start, y = rep(.2, length(start)), pch = 1, col = "darkblue", cex = .5))

with(dropstats[dropstats$n.privateSf2>0 & dropstats$t.test_qvalue<=0.05,], 
     points(x = start, y = rep(.6, length(start)), pch = 8, col = "green"))

with(dropstats[dropstats$n.privateSf2==0 & dropstats$t.test_qvalue<0.05,], 
     points(x = start, y = rep(.4, length(start)), pch = 1, col = "darkred", cex = .5))
with(tdna2, segments(x0 = x, x1 = pos, y0 = .7,y1 = .6, lty = 3))

with(tdna2, text(x = x,y = .72, labels = Gene, srt = 90, adj = 0, cex = .8))

axis(1, at = c(22420000,22760000,23100000), labels = c(22.42000,22.76000,23.10000))
axis(2, at = c(0,.2,.4,.6), 
     labels = c("genes",
                ">=1 SF2 Private SNP", 
                "SF2 Expression P <= 0.05",
                "Private SNP and Expression"),
     las = 2, line=NULL)



par(mar = c(5, 4, 4, 2) + 0.1)
```


```{r}
oli2<-counts[counts$gene == "AT5G55920",-1]
oli2<-data.frame(parents = names(oli2), exp = as.numeric(oli2), stringsAsFactors = F)
oli2$rank<-rank(oli2$exp)
with(oli2, plot(rank, exp))
with(oli2, text(x = rank, y = exp, label = parents, adj = c(1,1)))
```

```{r}
x=sqrt(mean.trf/100)
head(tel)
hist(tel$mean.trf/1000, breaks = 40)
with(tel, hist(sqrt(mean.trf)*mean.trf,mean.trf))
hist(log2(tel$mean.trf/1000), breaks = 40, xaxt = "n", xlab = "telo.length")
axis(1,at = log2(2:8), labels = 2:8)
```

