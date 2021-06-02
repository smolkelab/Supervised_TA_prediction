logFPKMs <- read.csv(file="D:/CME 250/aba.matrix.FPKM.vf.082511.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
logFPKMs$X <- NULL
colnames(logFPKMs)[3:13] <- unlist(lapply(logFPKMs[1,c(3:13)], as.character), use.names=FALSE)
logFPKMs <- na.omit(logFPKMs[-1,])
logFPKMs[,3:13] <- lapply(logFPKMs[,3:13], function(x) as.numeric(as.character(x)))
logFPKMvals <- data.matrix(logFPKMs[,3:13])
rownames(logFPKMvals) <- paste(logFPKMs$Locus.ID, logFPKMs$Unified.Functional.Annotation, logFPKMs$PFAM, sep=">")
candidates <- logFPKMvals[grep("efflux|transporter|PF00854|PF16974|PF01554|PF00664|PF00005|ABC|littorine|hyoscyamine|putrescine|tropinone|locus_4635", 
                               rownames(logFPKMvals), ignore.case=TRUE),]

#Using CYP80F1 littorine monooxygenase as bait gene
baitdata <- colMeans(candidates[grep("littorine", rownames(candidates), ignore.case = TRUE),]) 
model <- apply(candidates, 1, function(x) summary(lm(baitdata~x))$coefficients[8])
CYP80F1_pvals <- data.frame(p=sapply(model, function(x) x[1]))
CYP80F1_pvals <- na.omit(CYP80F1_pvals[order(CYP80F1_pvals$p), , drop=FALSE])

#Using hyoscyamine 6b-hydroxylase/oxygenase as bait gene
baitdata <- colMeans(candidates[grep("hyoscyamine", rownames(candidates), ignore.case = TRUE),]) 
model <- apply(candidates, 1, function(x) summary(lm(baitdata~x))$coefficients[8])
H6H_pvals <- data.frame(p=sapply(model, function(x) x[1]))
H6H_pvals <- na.omit(H6H_pvals[order(H6H_pvals$p), , drop=FALSE])

# Take hits for each of the bait genes and trim duplicates, then compute product of the p-values
top_hits <- data.frame(ID=unique(c(rownames(CYP80F1_pvals), rownames(H6H_pvals)))) #compile lists
top_hits$ID <- sub('>.*', '', top_hits$ID) #remove everything but locus IDs from rownames
top_hits$combined_p <- sapply(seq(1:length(top_hits$ID)), function(x) 
  log10(CYP80F1_pvals$p[grep(top_hits$ID[x], rownames(CYP80F1_pvals), ignore.case = TRUE)])
  + log10(H6H_pvals$p[grep(top_hits$ID[x], rownames(H6H_pvals), ignore.case = TRUE)]))
top_hits <- top_hits[order(top_hits$combined_p), , drop=FALSE] #order by combined p value
top_hits <- top_hits[top_hits$combined_p < -2,] #drop any with combined P < 0.01

# Generate subset of original log2FPKM values with the top_hits candidates
top_hits_FPKMs <- data.matrix(logFPKMs[sapply(top_hits$ID, function(x)
  grep(x, logFPKMs$Locus.ID, ignore.case = TRUE)), 3:13, drop=FALSE])
rownames(top_hits_FPKMs) <- top_hits$ID
#top_hits_FPKMs <- 2^top_hits_FPKMs

#Generate list of known TA bait genes and experimentally validates TA transporters
TA_gene_IDs <- read.table(file="D:/CME 250/TA_gene_locus_IDs.txt", header=FALSE, stringsAsFactors = FALSE)
rownames(top_hits_FPKMs) <- sapply(rownames(top_hits_FPKMs), function(x) 
  { if( length(grep(x, TA_gene_IDs$V1)) > 0 ) {TA_gene_IDs$V2[grep(x, TA_gene_IDs$V1)]} else {x} })

# Generate heatmap for top hits, non-normalized
library(gplots)
library(RColorBrewer)

row_labels <- sapply(rownames(top_hits_FPKMs), function(x) 
{ if( length(grep("locus", x)) > 0 ) {""} else {x} })

heatmap.2(top_hits_FPKMs, key=TRUE, col=colorRampPalette(c('red', 'black', 'green')), Colv=FALSE, dendrogram="row",
          scale="row", trace="none", sepwidth=c(0,0), density.info = 'none', 
          key.title = NA, colsep=1:ncol(top_hits_FPKMs), rowsep=1:nrow(top_hits_FPKMs),
          labRow = row_labels, cexRow = 0.3, lhei = c(1,7))
