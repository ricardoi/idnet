#rSDT
#Is a species demarcation tool implemente in R (rSDT) 
#following a similar framework as Species Demarcation Tool (SDT)  
# but adapted for R users

#Packages needed from Bioconductor and CRAN
#Multiple Sequence Alignment for R (https://bioconductor.org/packages/devel/bioc/vignettes/msa/inst/doc/msa.pdf)
source("http://www.bioconductor.org/biocLite.R")
biocLite("msa")
biocLite("Biostrings")
install.packages("heatmaply")
library(msa)
# require(Biosntrings)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(heatmaply)
library(ape)
library(sm)
setwd("/Users/RicardoI/++Manuscripts/++2017/rSDT/")
setwd("/Users/ricardoialcala/Documents/++Manuscripts/++2017/rSDT/")

#It is automatically performed at the beginning of an R session and may be used to initialize the environment
.First <- function() {
  options(digits=2, length=999)
  library(msa)
  require(Biostrings)
  library(Biostrings)
  library(seqinr)
  library(ggplot2)
  library(heatmaply)
  library(ape)
  library(sm)
}


#requires Biostrings

mySequences <- readDNAStringSet("SCMV_rSDT/SCMV_2017.fasta")
mySequences

#aligning with msa
#myFirstAlignment <- msa(mySequences, "Muscle")

#aligning with msaMuscle FASTER THAN msa
myFirstAlignment <- msaMuscle(mySequences, gapOpening = 420, verbose=TRUE)

#printing
print(myFirstAlignment)
print(myFirstAlignment, showConsensus=FALSE, halfNrow=5)


#I need to convert my data to a seqin file 
mySecondAlignment <- msaConvert(myFirstAlignment, type="seqinr::alignment")

#Similarity EL BUENO
dd <- dist.dna(as.DNAbin(mySecondAlignment, model = "K80", as.matrix= TRUE))
dd
save(dd, file="SCMV_similarity.Rdata")

per2plot <- function(x){
  x= round(as.matrix(-(x-(1))*100), digits=2)
  x[ x == 0] <- 100
  (x)
}
dd2plot= per2plot(dd)
# This two line were combined in a function to produce a percentage matrix
#dd2plot = round(as.matrix(-(dd-(1))*100), digits=2)
#dd2plot[ dd2plot == 0] <- 100

#Density plot 
plot(density(dd2plot[1,], cut=0.22), col="red", border="blue", main= "Density plot")


#Density plot 
# Base plot
d1= density(dd2plot[1,])
plot(d1,  main="pairwise similairity distribution")
polygon(d1, col="tomato", border="blue")

### prepraring data for ggplot
dd2plot.d <- as.data.frame(dd2plot[1,])
names(dd2plot.d)[1]<- "similarity"     #borrar###dd2plot.d$group <- list

# ggplot
ggplot(dd2plot.d, aes(x = similarity))+ 
  geom_density(fill="tomato", alpha = 0.5)+
  ggtitle("Pairwise Similairity Distribution")+
  theme(plot.title = element_text(hjust=0.5, color="black", size=16, face="bold"))


## Este saca sequence genetic distance otra cosa
da <- as.matrix(dist.alignment(mySecondAlignment, "similarity"))
da

#heat map 
heatmaply(as.matrix(dd2plot), margins = c(175, 200), colors= c("blue", "green", "yellow", "red"), fontsize_row = 5, fontsize_col = 5)


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#half heat map
library(reshape2)

halfdd<- get_upper_tri(dd2plot)
melt_halfdd <- melt(halfdd, na.rm= TRUE)

#ggplot for rSDT
ggplot(data = melt_halfdd, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                       #midpoint = mean(melt_halfdd$value), space = "Lab", 
                       name="Pairwise\nSimilarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()


 rSDT <- function(x){
   x <- get_upper_tri(x)
   y <- melt(as.matrix(x), na.rm= TRUE)
   ggplot(data = y, aes(Var2, Var1, fill = value))+
     geom_tile(color = "white")+
     scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                          #midpoint = mean(melt_halfdd$value), space = "Lab", 
                          name="Pairwise\nSimilarity") +
     theme_minimal()+ 
     theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                      size = 8, hjust = 1))+
     coord_fixed()
 }

rSDT(dd2plot)

#aqui estoy trratando de meter los comandos como me dijo Masata para poder manipularlos 
#en la función rSDT
x= mySecondAlignment
function.rSDT <- function(data, dist.dna = TRUE, model = "", size =""){
  if(dist.dna == TRUE) { 
    x <- dist.dna(as.DNAbin(data, 'model', model, as.matrix= TRUE)) 
    x= round(as.matrix(-(x-(1))*100), digits=2)
    #x[ x == 0] <- 100
    diag(x) <-100
    (x)  
    y <- get_upper_tri(x)
    z <- melt(y, na.rm= TRUE)
  } 
  else {
    y <- as.matrix(get_upper_tri(data))
    z <- melt(y, na.rm= TRUE)
  }
  ggplot(data = z, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                         #midpoint = mean(melt_halfdd$value), space = "Lab", 
                         name="Pairwise\nSimilarity") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 5 , hjust = 1),
          panel.grid.major = element_blank(), #look for a way to turn off or on the grids 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    coord_fixed()
}
#diag(DsMV.dis[ ,3:129])
function.rSDT(DsMV.dis[ ,3:129], dist.dna = FALSE, mode = 'K80')
dev.off()
function.rSDT(mySecondAlignment, dist.dna = TRUE, mode = 'K80')
nrow(DsMV.dis)
#aln = alignment = mySecondAlignment
str(mySecondAlignment)

library(ggdendro)
?dendro_data

dendro_data(SCMVtree, type="rectangle")


plot(density(melt_halfdd$value))

#N-J tree
SCMVtree <- nj(dd)
plot(SCMVtree, main="Phylogenetic Tree of SCMV WG Sequences")
#function que hizo Masato
Density.plot <- function(data){
  data2 <- get(data)
  plot(density(data2[,1], cut = 0), 
       main = paste("Desnit plot ", data, sep = ""))
}




#Genetic distances
gd <- dist.gene(as.matrix(mySecondAlignment), "pairwise")
class(gd)

#plotting different groups
list <- as.numeric(c("1","2","2","2","2","2","2","2","1","1","1","3","1","3","3","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","3","1","1","1","2","2","2"))
dim(list)
class(list)
dd.all = (ddd[1,])
dd.sc= subset(ddd[1,], rownames(ddd)  %in% sugarcane)
sm.density.compare(dd.all, list, model= "equal", xlab="percentage of similarity")


DsMVsdt <- read.csv("DsMV_pairwisematrix_codon.csv", row.names= 1, header = TRUE)
DsMV_halfdd<- get_upper_tri(DsMVsdt)
DsMVmelt_halfdd <- melt(as.matrix(DsMV_halfdd), na.rm= TRUE)


rSDT(DsMVsdt)
