#Analysing Illumina-microarrays

##Requirements
Install the limma package from Bioconductor.
```r
source(“http://bioconductor.org/biocLite.R”)
biocLite(“limma”)
```

##Quality control
Here we load the raw data coming fro the Illumina Beadstudio software and crete several plots from the raw and normalised data.

Set the working directory where the data is located and load the limma pack
```r
setwd("~/projects/012345/data")
library(limma)
```

Load the Illumina data and control files. 
```r
x <- read.ilmn("Group_Probe_Profile.txt", ctrlfiles=c("TableControl_lps data.txt", "TableControl.txt"), other.columns="Detection")
```

Check type of probes on the microarray-chip
```r
table(x$genes$Status)

# Example output:
#    BIOTIN            CY3_HYB       HOUSEKEEPING           LABELING LOW_STRINGENCY_HYB
#         2                  6                 14                  8                  8
#  NEGATIVE            regular
#       936              45281
```

Estimate the proportion of expressed probes for each of the arrays.
```r
propexpr(x, labels=c("NEGATIVE", "regular"))

# Example output:
# unstimulated_rep1  unstimulated_rep2  unstimulated_rep3           LPS_rep2 Batf2 KD_ LPS_rep2
#         0.4122342          0.4392318          0.3782704          0.3745475          0.3565128
#          LPS_rep1  Batf2 KD_LPS_rep1           LPS_rep3  Batf2 KD_LPS_rep3  unstimulated_rep4
#         0.3648191          0.3709147          0.3422734          0.3539387          0.3943150
#          LPS_rep4  Batf2 KD_LPS_rep4  unstimulated_rep5           LPS_rep5  Batf2 KD_LPS_rep5
#         0.3754312          0.3988355          0.3968508          0.3656047          0.4068437
```

Create some plots for the raw data.
```r
# Boxplot of regular probes and control
png(file="raw_boxplot.png", width=1200, height=1200, units="px")
par(mar=c(10, 5, 5, 5))
boxplot(log2(x$E[x$genes$Status=="regular",]),range=0,las=2,ylab="log2 raw intensities", main="Regular probes")
dev.off()

png(file="raw_boxplot_negativeControls.png", width=1200, height=1200, units="px")
par(mar=c(10, 5, 5, 5))
boxplot(log2(x$E[x$genes$Status=="NEGATIVE",]),range=0,las=2,ylab="log2 raw intensities", main="Negative control probes")
dev.off()

png(file="raw_mds.png", width=1200, height=1200, units="px")
par(mar=c(10, 5, 5, 5))
plotMDS(x,labels=x$targets$SampleNames, main="Raw probes (all)")
dev.off()
```

Normalise the data. Perform normexp-by-control background correction, quantile normalization and log2-transformation to the raw data:
```r
y <- neqc(x, negctrl="NEGATIVE")
```
Filtering out non-detected probes based on the p-Value of detection. Filter out probes which were not detected in at least one sample with detection pV=0.01.

**Attention! This depends on the version of the Beadstudio software.**

```r
# After Version 2 of Beadstudio software
detected <- apply(y$other$Detection < 0.01, 1, any)
# Version 2 or earlier of Beadstudio uses 1-Th in file
detected <- apply(y$other$Detection > 0.99, 1, any)

# Subselect setected probes form the normalised set.
y <- y[detected,]
```

Create plots for the normalised data.
```r
png(file="norm_boxplot.png", width=1200, height=1200)
par(mar=c(10, 5, 5, 5))
boxplot(y$E,range=0,las=2,ylab="normalised intensities", main="Regular probes", units="px")
dev.off()

par(mar=c(10, 5, 5, 5))
png("norm_mds_detected_regular_probes.png", width=1200, height=1200)
plotMDS(y,labels=y$targets$SampleNames, main="Normalised probes (detected)")
dev.off()

#
# END
#
```

Now you can investigate the plotsin the working directory for strange artifacts.

##Calling differential expressed genes
```r
# 
# Set working directory with the data
#
setwd("~/projects/012345/data")

#
# Load limma package
#
library(limma)

#
# LOAD DATA AND CONTROL PROBE
x <- read.ilmn("Group_Probe_Profile.txt", ctrlfiles=c("TableControl_lps data.txt", "TableControl.txt"), other.columns="Detection")
```

Sometimes QC let us to remove a replicate. If you know the column numbers of the ones that are ok, this can be achieved like:
```r
# all columns of replicates that we want to include
E <- subset(x$E, select=c(1,3,6,7,8,9,10,11,12,13,14,15))  # subset on expression
D <- subset(x$other$Detection, select=c(1,3,6,7,8,9,10,11,12,13,14,15)) # subset on detection values
x$E <- E
x$other$Detection <- D
```

Normalise the expression.Perform normexp background correction using negative control probes and quantile normalization using negative and positive control probes. After normalization, the intensities are log2 transformed and the control probes are removed.
```r
y <- neqc(x, negctrl="NEGATIVE")
```

Filtering out non-detected probes based on the p-Value of detection. Filter out probes which were not detected in at least one sample with detection pV=0.01.

**Attention! This depends on the version of the Beadstudio software.**

```r
# Keep probe if at leat 1 replicate has the probe as detected

# After Version 2 of Beadstudio software
#detected <- apply(y$other$Detection < 0.01, 1, any)
# Version 2 or earlier of Beadstudio uses 1-Th in file
detected <- apply(y$other$Detection > 0.99, 1, any)

# Get some stats.
length(detected[detected==TRUE]) # print number detected
length(detected[detected==TRUE])*100/length(y$E[,1]) # print percentage detected
yD <- y[detected,] # subselect detected

#
# WRITE NORMALISED EXPRESSION OF DETECTED PROBES TO FILE
#
write.table(yD, file="RESULTS_EXPR_NORM_LOG2_DETECTED.txt", sep="\t", quote=FALSE, row.names=FALSE)
```

Now create groups for the differntial expression analysis. These are the ones that will be compared to each other.

```r
#
# CREATE GROUPS FOR DIFF EXPR ANALYSIS
#
groups <- c("unstimulated",  # unstimulated_rep1
            "unstimulated",  # unstimulated_rep3
            "LPS",           # LPS_rep1
            "Batf2_KD_LPS",  # Batf KD_ LPS_rep1
            "LPS",           # LPS_rep3
            "Batf2_KD_LPS",  # Batf KD_ LPS_rep3
            "unstimulated",  # unstimulated_rep4
            "LPS",           # LPS_rep4
            "Batf2_KD_LPS",  # Batf KD_ LPS_rep4
            "unstimulated",  # unstimulated_rep5
            "LPS",           # LPS_rep5
            "Batf2_KD_LPS")  # Batf KD_ LPS_rep5
```

Create the design matrix and contrasts based on experimental setup and experiments.
```r
ct <- factor(groups) # factors
design <- model.matrix(~0+ct) # design matrix
colnames(design) <- levels(ct) # set column names of design matrix
# Here we define what we want to compare
contrasts <- makeContrasts(LPS-unstimulated, Batf2_KD_LPS-LPS, levels=design) # create contrast matrix
```

Now fit the data using limma package
```r
fit <- lmFit(yD$E,design) # fit
contrasts.fit <- eBayes(contrasts.fit(fit, contrasts))
```

Get the results for all probes and write them into files.
```r
t1.res   <- topTable(contrasts.fit, coef=1, n=1000000) # get all
# associate targetIDs to probes
t1.final <- merge(t1.res, y$genes, by.x="ID", by.y="ProbeID")
# sort based on adjusted p-Value
t1.final <- t1.final[order(t1.final$adj.P.Val),] 

t2.res   <- topTable(contrasts.fit, coef=2, n=1000000) # get all
# associate targetIDs to probes
t2.final <- merge(t2.res, y$genes, by.x="ID", by.y="ProbeID") 
# sort based on adjusted p-Value
t2.final <- t2.final[order(t2.final$adj.P.Val),]

# WRITE THE RESULTS OF THE TWO COMPARISIONS INTO A FILES
write.table(t1.final, file="RESULTS_LPS_VS_UNSTIMULATED.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(t2.final, file="RESULTS_BATF2KD-LPS_VS_LPS.txt", sep="\t", quote=FALSE, row.names=FALSE)

#
# END
#

```
