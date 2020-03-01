library(GEOquery)
library(oligo)

### Downloading files
gse <- getGEO(filename='./GSE124483_family.soft')
gsm_list <- names(GSMList(gse))
filePaths <- getGEOSuppFiles('GSE124483')

celfiles <- list.files("GSE124483_gz")
rawData <- read.celfiles(paste0('GSE124483_gz/',celfiles))

gsm <- getGEO("GSE124483",GSEMatrix=FALSE)


### Attaching pData to raw data
names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))]
Meta(gsm)[!(names(Meta(gsm)) %in% names(Meta(gse)))]


genotype <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}

age <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][3]
}

pd <- data.frame(age=as.factor(sapply(GSMList(gse),age)), genotype=as.factor(sapply(GSMList(gse),genotype)))
pd$age <- as.factor(pd$age)
levels(pd$age) <- c("old","young")
levels(pd$genotype) <- c("KO", "wild")
pData(rawData)$genotype <- as.factor(pd$genotype)
pData(rawData)$age <- as.factor(pd$age)

### Data normalisation
normData <- rma(rawData)

par(mfrow=c(1,2))
boxplot(rawData, target="core", main="Raw Data")
boxplot(normData, target="core", main="Normalised Data")


library(limma)

group <- with(pData(normData), paste(genotype, age, sep = "_"))
group <- factor(group)

design <- model.matrix(~ 0 + normData$genotype*normData$age)
colnames(design) <- levels(group)
cm <- makeContrasts(genotype_old = wild_old - KO_old,
                    genotype_young = wild_young - KO_young,
                    age_KO = KO_old - KO_young,
                    age_wild = wild_old - wild_young,
                    interaction = (KO_old - KO_young) - (wild_old - wild_young),
                    levels = design)

# Fitting model
fit <- lmFit(normData, design)

#Fitting contrasts
fit2 <- eBayes(contrasts.fit(fit,contrasts = cm))

ps <- rownames(topTable(fit2))

#Summarise results
results <- decideTests(fit2)
summary(results)

library(mta10transcriptcluster.db)
library(RColorBrewer)
head(keys(mta10transcriptcluster.db,keytype="PROBEID"))
AnnotationDbi::select(mta10transcriptcluster.db,ps,c("SYMBOL","ONTOLOGY", "ENTREZID","GENENAME"),keytype="PROBEID")

par(mfrow=c(1,1))
volcanoplot(fit2)

interesting_genes <- topTable(fit2,number=Inf,p.value = 0.05,lfc=2)

eset_of_interest <- normData[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest),labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

### MDS plot for troubleshooting
plotMDS(normData,col=c(rep("red",3),rep("blue",3),rep("orange",3), rep("black",3)))
