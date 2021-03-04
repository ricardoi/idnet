# rSDT_package source
# Set of functions to plot SDTv1.8.py outputs using ggplot in R.
rm(list=ls(all=TRUE)) 

# -------------------------------------------------------------------
# load libraries
.First <- function() {
  options(digits=2, length=999)
  library(msa)
  require(Biostrings)
  library(Biostrings)
  library(seqinr)
  library(ggplot2)
  library(heatmaply)
  library(reshape2)
  library(ape)
  library(sm)
}
# -------------------------------------------------------------------
# SDT output from python, is a melted kind of format which requires
# to be formatted to make it compatible with ggplot
# -------------------------------------------------------------------
# Functions
#
# this function is to round matrices to two digits in percentage 
per2plot <- function(x){
  x= round(as.matrix(-(x-(1))*100), digits=2)
  x[ x == 0] <- 100
  (x)
}
# -------------------------------------------------------------------
# Functions for triangles
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
# -------------------------------------------------------------------
# This funtion is trasnformer 
sdt2plot <- function(x){
  n <- x$V1
  ij <- x[,-1] 
  i <- as.matrix(ij) 
  j <- as.matrix(t(ij))
  k <- matrix(NA, nrow(ij), ncol(ij)) #empty matrix to fill in
  k[lower.tri(k)] <- i[lower.tri(i)]
  k[upper.tri(k)] <- j[upper.tri(j)]
  diag(k) <- 1
  colnames(k) <- n
  rownames(k) <- n
  (k)
}
# -------------------------------------------------------------------
# This function is to plot matrices with pairwise, similarity or identity data 
# rSDT a function to plot matrices from SDT
rSDT <- function(x, size="", title="", angle.x=""){
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  x <- get_upper_tri(x)
  y <- melt(as.matrix(x), na.rm= TRUE)
  ggplot(data = y, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                         #midpoint = mean(melt_halfdd$value), space = "Lab", 
                         name="Pairwise\nSimilarity") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = angle.x, vjust = 1, size = size, hjust = 1),
          axis.text.y = element_text(vjust = 1, size = size, hjust = 1),
          panel.grid.major = element_blank(), #look for a way to turn off or on the grids 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    labs(title=title, 
         #subtitle="",
         caption="Identity plot: SDT v1.0 ",
         x="",
         y="")+
    coord_fixed()
}

