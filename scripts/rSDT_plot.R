#-------  plot results from SDT_Mac64.py
# load packages
library(plyr)
library(ggplot2)
library(reshape2)
library(psych)
#
#
setwd("/Applications/SDT_Mac64/output/")
#--------------------------------- functions ---------------------------------
# read.SDT file
read.SDT <- function(x){ 
              m = list(NULL)
                for (i in 1:nrow(x)) {
                    m[i] <- strsplit(x[i,], split = ',')
                    }
              m <- ldply(m, rbind)
            rownames(m) <- as.character(m[,1])
            colnames(m) <- c("names", as.character(m[,1]))
            x <- as.matrix(m[-1])
            x[is.na(x)] <- 0
            x <- apply(as.matrix(x),2,as.numeric)
            rownames(x) <- colnames(x)
            (x)
}

#------------------- plotSDT
plot.SDT <- function(y, size = ""){
  z <- melt(as.matrix(y), na.rm= TRUE)
  z$value <- as.numeric(as.character(z$value))
  z[z==0] <- NA
  ggplot(data = z, aes(Var2, Var1, fill = value*100))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                         na.value="white",
                         #midpoint = mean(melt_halfdd$value), space = "Lab", 
                         name="Pairwise\nSimilarity") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = size  , hjust = 1),
          axis.text.y = element_text(angle = 0, vjust = 1, size = size, hjust = 1),
          panel.grid.major = element_blank(), #look for a way to turn off or on the grids 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    coord_fixed()
}
#-------------------------- Load data ---------------------------
# --------------------------- RNA1 ---------------------------
ORF.up <- read.table(file='Comovirus_ORF_RNA1_ord.aln_mat.txt', as.is = T, stringsAsFactors = F)
pol.do <- read.table(file='Comovirus_polyprotein_RNA1_ord.aln_mat.txt', as.is = T, stringsAsFactors = F)

# SDT
(mUP <- read.SDT(ORF.up))
(mDO <- read.SDT(pol.do))

RNA1 <- lowerUpper(t(mUP), mDO)
diag(RNA1) <- 1
plot.SDT(RNA1*1, size = 12)

# --------------------------- RNA1 ---------------------------
ORF2.up <- read.table(file='Comovirus_ORF_RNA2_reord.aln_mat.txt', as.is = T, stringsAsFactors = F)
pol2.do <- read.table(file='Comovirus_polyprotein_RNA2_reord.aln_mat.txt', as.is = T, stringsAsFactors = F)

# SDT
(mUP2 <- read.SDT(ORF2.up))
(mDO2 <- read.SDT(pol2.do))

RNA2 <- lowerUpper(t(mUP2), mDO2)
diag(RNA2) <- 1
plot.SDT(RNA2*1, size = 12)

#------------------------- Species Demarcation ---------------------
# ---- POL ---- 
POL <- read.table(file='Comovirus_pol.fasta_mat.txt', as.is = T, stringsAsFactors = F)
mPol <- read.SDT(POL)
plot.SDT(mPol*1, size = 12)

# ---- Pro-Pol ---- 
DOM <- read.table(file='Comovirus_pro-pol_domain.fasta_mat.txt', as.is = T, stringsAsFactors = F)
mDOM <- read.SDT(DOM)
plot.SDT(mDOM*1, size = 12)

# ---- CP ----
CP <- read.table(file='Comovirus_CP.fasta_mat.txt', as.is = T, stringsAsFactors = F)
mCP <- read.SDT(CP)
plot.SDT(mCP*1, size = 12)
