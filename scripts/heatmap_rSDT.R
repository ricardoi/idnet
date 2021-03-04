#########################################################
### A) Installing and loading required packages
#########################################################

if(!require(source("http://www.bioconductor.org/biocLite.R"))){
  biocLite("msa")
  biocLite("Biostrings") 
}
if (!require("ape")) {
  install.packages("ape", dependencies = TRUE)
  library(ape)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}
if (!require("reshape")) {
  install.packages("reshape", dependencies = TRUE)
  library(reshape)
}
if (!require("ggdendro")) {
  install.packages("ggdendro", dependencies = TRUE)
  library(ggdendro)
}
if (!require("plotly")) {
  install.packages("plotly", dependencies = TRUE)
  library(plotly)
}
}
if (!require("gridBase")) {
  install.packages("gridBase", dependencies = TRUE)
  library(grid)
}

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

setwd("/Users/RicardoI/++Manuscripts/++2017/rSDT/")
DsMVsdt <- read.csv("DsMV_pairwisematrix_codon.csv", stringsAsFactors = F, header = TRUE)
# Remove metadata to plot: heat map 
x <- DsMVsdt[,3:ncol(DsMVsdt)]
row.names(x)  <- colnames(x)
y <- x 
head(y) #aqui tengo que cambiar la y por la x en el script, esto es para afinar el dendograma

#Identity score
DsMVid <- read.csv("DsMV_FlplusWorld_2016_Final00_ordenadas_NoGap.fasta_mat.csv", row.names = 1, header = T)

head(DsMVid)
x = DsMVid
#here it works with headers
i <- as.matrix(DsMVid) 
j <- as.matrix(t(DsMVid))
k <- matrix(NA, nrow(DsMVid), ncol(DsMVid)) #empty matrix to fill in
k[lower.tri(k)] <- i[lower.tri(i)]
k[upper.tri(k)] <- j[upper.tri(j)]
diag(k) <- 1
colnames(k) <- colnames(DsMVid)
rownames(k) <- rownames(DsMVid)
tail(k)

#here TRY to make it work without mods
DsMVid2 <- read.csv("DsMV_FlplusWorld_2016_Final00_ordenadas_NoGap.fasta_mat_orig.csv", as.is = T, header = F)
x = DsMVid2
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
DsMV.sdt <- sdt2plot(DsMVid2)


#to put in the function
names <- DsMVid$V1 #getting the matrix names
i <- as.matrix(DsMVid[,-1]) #lower triangle
j = t(i) #upper triangle
k <- matrix(NA, nrow=106, ncol=106) #empty matrix to fill in
k[lower.tri(k)] <- i[lower.tri(i)]
k[upper.tri(k)] <- j[upper.tri(j)]
diag(k)<-1
rownames(k) <- names
colnames(k) <- names 
head(k)


#########################################################
### C) Customizing and plotting the heat map
#########################################################

function.rSDT <- function(x, dist.dna = TRUE, dendro = TRUE, model = "", size =""){{
  if(dist.dna == TRUE) { 
    x <- dist.dna(as.DNAbin(x, 'model', model, as.matrix = TRUE)) 
    x= round(as.matrix(-(x-(1))*100), digits=2)
    #x[ x == 0] <- 100
    diag(x) <-100
    #(x) #esto se puede borrar o mostrar un avance con print
    x[lower.tri(x)]<- NA
    #  y <- get_upper_tri(x)
    x <- melt(x, na.rm= TRUE)
  }
  else {
    x[lower.tri(x)]<- NA
    x <- melt(as.matrix(x), na.rm= TRUE)
  }
  #p1 <- 
    ggplot(data = x, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                         #midpoint = mean(melt_halfdd$value), space = "Lab", 
                         name="ID %") +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1))+
    #theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 5 , hjust = 1),
          axis.text.y = element_text(angle = 0, vjust = 1, size = 5 , hjust = 1),
          panel.grid.major = element_blank(), #look for a way to turn off or on the grids 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    coord_fixed()
  } 
}
function.rSDT(DsMV.sdt, dist.dna = FALSE, dendro = FALSE)


# I need to figure out how to put it in my fucntion
  #if(dendro == TRUE ){
    dd.col <- as.dendrogram(hclust(dist(x)))
    dd.row <- as.dendrogram(hclust(dist(x)))
    col.ord <- order.dendrogram(dd.col)
    row.ord <- order.dendrogram(dd.row)
    ddata_x <- dendro_data(dd.col)
    ddata_y <- dendro_data(dd.row)
    h <- x[col.ord, row.ord]
    h[lower.tri(h)]<- NA
    h <- melt(as.matrix(h), na.rm= TRUE)
p2 <- ggplot(segment(ddata_x)) + 
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
      #theme_none + 
      theme(axis.title.x= element_blank(),
            axis.title.y= element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(), 
            panel.background = element_blank()
      )
p3 <- ggplot(segment(ddata_y)) + 
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
      coord_flip() + 
      theme(axis.title.x= element_blank(),
            axis.title.y= element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(), 
            panel.background = element_blank()
      )
p4 <- ggplot(data = h, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradientn(colours= c("blue", "green", "yellow", "red"),
                           #midpoint = mean(melt_halfdd$value), space = "Lab", 
                           name="ID %") +
      theme(legend.position = c(0, 1), legend.justification = c(0, 1))+
      #theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 5 , hjust = 1),
            axis.text.y = element_text(angle = 0, vjust = 1, size = 5 , hjust = 1),
            panel.grid.major = element_blank(), #look for a way to turn off or on the grids 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      coord_fixed()
      grid.newpage()
      print(p4, vp=viewport(0.8, 0.8, x=0.4, y=0.4))
      print(p2, vp=viewport(0.503, 0.14, x=0.432, y=0.87)) # I need to automitize this
      print(p3, vp=viewport(0.14, 0.76, x=0.735, y=0.45)) # I need to automitize this
  #}
#}

      
      

