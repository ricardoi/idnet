---
title:   "Network ananlysis for Demarcation Species of Viruses"
output: 
word_document: default
pdf_document: default
html_document: default
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
#setwd("/Users/RicardoI/Documents/Documents-ehecatl/git_db/id_net/")
library(reshape2)
library(ggplot2)
library(ape)
library(igraph)
```

R. I. Alcala-Briseno, Plant Pathology Department, University of Florida.

## A network analysis of viral identities for demarcation of species

Viruses are classified using sequence iditities, several methods had been used to assign viral species. 
This package is designed to visualize the pairwise similarity matrix in a network visualization, taking advantage of the statistical framewokr theory of netwokr analysis testing likelihoods, hierarchichal clustering, max flow, etc... 

The package can manage output formats from a different of alignment programs, such as muscle, which can be run in R or locally in your computer; outputs exported from gui's like geneious or CLC bench; or genereted with the program Species Demarcation Tool (SDT, Varsani et al. 201?). Basically this matrices are loaded to R and converted to a melt format to manipulated it with igraph. Each of the nodes is a virus and the links are the ideintities shared among them, besides taking advantage of statistical network theory.    
Also, a function to generate the classical heat map for the Species Demarcation windows Tool has been developed to be plotted in R.

The functions introduced here are part of an R package under development, "iNet", for Viral Identities Network Analysis (VINA).  VINA is the analysis of taxonomic characteristics used as information for the characterization of viral species (Alcala-Briseno and Garret 2018).

## Generating an output file from pairwiese similarities.



# SDT
The sdt calculates the identity scores are as 1-(M/N) where M is the number of mismatching nucleotides and N the total number of positions along the alignment at which neither sequence has a gap character.

Running SDT on bash

<font size="1"> py2 SDT_Mac64.py test.fas muscle </font>  

Three outputs will be generated, but we are going to load only number 3

1- test.fas.sdt	#SDT file 

2- test.fas_identity_scores.txt	

3- test.fas_mat.txt # load the matrix with identity scores to R

```{r echo = TRUE}
DsMVmat <- read.csv("DsMV_2017_final.fasta_mat.csv", stringsAsFactors = F, header = F)
dim(DsMVmat)
``` 

The output 3 can be utilized for iNet, however, output 2 needs to be transformed or output 3 can be loaded directly to iNet. 

# Loading or converting the SDT outpot to a plotting format with sdt2plot( )
To load directly ot iNet use the file from SDT ending with "_mat".

```{r echo = TRUE}
sdt2plot <- function(x, diag = " "){
  n = DsMVmat
  n <- x$V1
  ij <- x[,-1] 
  i <- as.matrix(ij) 
  j <- as.matrix(t(ij))
  k <- matrix(NA, nrow(ij), ncol(ij)) #empty matrix to fill in
  k[lower.tri(k)] <- i[lower.tri(i)]
  k[upper.tri(k)] <- j[upper.tri(j)]
  diag(k) <- diag
  colnames(k) <- n
  rownames(k) <- n
  round((k)*100, 2)
}
```

```{r echo = TRUE}
DsMVmat1 <- sdt2plot(DsMVmat, diag = 0)
dim(DsMVmat1)
DsMVmat1[1:6,1:3]
```
Note that the diagonal contains zeros, because we are not interested in the diagonal, later on this values will be removed before load it into igraph.

## iNet graph 

The identity matrix is loaded as an iGraph object, the iNet function produce a network of viral species, represented in the number of nodes. Demarcation thresholds per each group of organism can be set up. We can add metadata as column for location, plant host, vector, etc.

Adding metada as csv, this metadata was collected from the NCBI or manuscript information. It is very important to have names columns as your pairwise input in order link both informations, metadata and pairwise similarity.
```{r echo=TRUE}
DsMVmd <- read.csv("DsMV_metadata_2017.csv", stringsAsFactors = F, header = TRUE)
dim(DsMVmd)
head(DsMVmd)
```

The function iNet follows this code:

```{r echo= TRUE}
DsMVnet <- melt(DsMVmat1) 
colnames(DsMVnet) <- c("From", "To", "value")
head(DsMVnet)
```
Formating the matrix; we are going to remove zeros that corresponds to the diagonal and we are gonna set a treshold, the demarcation species trheshold for potyvirus is 78 at nucleotide level and 79.6 at aminoacid identity. However, to see relationship among viral isolate, the threeshold can be modified for visual inspection. 
```{r echo = TRUE}
DsMVnet0 <- DsMVnet[!rowSums(DsMVnet[-c(1:2)] == 0) >= 1,] #removing 0's
DsMVnett <- DsMVnet0[DsMVnet0[-c(1:2)] > 85,] # setting up the treshold
head(DsMVnett)
```

```{r echo = FALSE}
#Quality check to check if the IDs for the pairwise simialrity matrix and the metadate match
#match(unique(DsMVnett$To), unique(DsMVmd$ID))
```

The pairwise identity matrix and the metada will be loaded in igraph, we don't need to specify the direction of the link, and for adding the metadata, we are going to specify in the vertices information. We are going to add a ramp palette to colour the identities (links) and the viral isolates (nodes). 
# Specifying the countries
```{r echo = TRUE}
#iNet <- function(x){
x = DsMVnett
igraph<-graph_from_data_frame(x, directed = FALSE, vertices=DsMVmd)
rbPal <- colorRampPalette(c("green", "yellow", "red")) 
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
V(igraph)$size=5
V(igraph)$xx <- as.numeric(as.factor(V(igraph)$Country)) # make the categories of x into numeric values for color ramp
V(igraph)$color <- counPal(10)[cut(as.numeric(V(igraph)$xx),breaks = 10)]
CounColor <- unique(cbind(V(igraph)$Country, V(igraph)$color))
E(igraph)$xx <- as.numeric(unlist(E(igraph))) # make the categories of x into numeric values for color ramp
E(igraph)$color <- rbPal(10)[cut(as.numeric(E(igraph)$xx),breaks = 10)]
#}
```
Working in the function but is not done yet. I need to figure out how to make a function loading igraph and the information to plot it.
```{r echo = TRUE}
#igraph<- iNet(DsMVnett) # it doesn't work yet
```
Plotting the igraph vector using the countries information
```{r echo = TRUE}
plot(igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)
```
# Specifying the host 
We are using the host information to plot the network, I can simplify this ... 
```{r echo = TRUE}
igraph<-graph_from_data_frame(x, directed = FALSE, vertices=DsMVmd)
V(igraph)$size=5
V(igraph)$xx <- as.numeric(as.factor(V(igraph)$Host)) # make the categories of x into numeric values for color ramp
V(igraph)$color <- counPal(10)[cut(as.numeric(V(igraph)$xx),breaks = 10)]
CounColor <- unique(cbind(V(igraph)$Host, V(igraph)$color))
E(igraph)$xx <- as.numeric(unlist(E(igraph))) # make the categories of x into numeric values for color ramp
E(igraph)$color <- rbPal(10)[cut(as.numeric(E(igraph)$xx),breaks = 10)]
```
Plotting the igraph vector using the host information... However, I can combine the country and host information in a single plot. 
```{r echo = TRUE}
plot(igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_mds)
legend(x=-0.9, y=-0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2, cex=.8, bty="n", ncol=1)
```



```{r echo = FALSE}

# all vertex shapes, minus "raster", that might not be available
shapes <- setdiff(shapes(), "")
#g <- make_ring(length(shapes))
set.seed(42)
#################################################################
# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}
# no clipping, edges will be below the vertices anyway
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=3))
```

# Combining countries and host 
We are using the host information to plot the network, I can simplify this ... 
```{r echo = TRUE}
igraph3<-graph_from_data_frame(x, directed = FALSE, vertices=DsMVmd)
V(igraph3)$size=5
V(igraph3)$xx <- as.numeric(as.factor(V(igraph3)$Country)) # make the categories of x into numeric values for color ramp
#shapes <- shapes(shape = NULL)
shapes <- paste("star", 2:13)
V(igraph3)$shape <- as.matrix(shapes)[1:length(unique(V(igraph3)$Host)),] 
### f(V(igraph3)$Host ==  1, #lo saque de la linea anterior
#V(igraph3)$Host <- as.numeric(as.factor(V(igraph3)$Host))
V(igraph3)$color <- counPal(10)[cut(as.numeric(V(igraph3)$xx),breaks = 10)]
#vcount((igraph3)$Host)
#V(igraph3)$color <- shape
CounColor <- unique(cbind(V(igraph3)$Country, V(igraph3)$color))
HostShape <- unique(cbind(unique(V(igraph3)$Host), unique(V(igraph3)$shape)))
E(igraph3)$xx <- as.numeric(unlist(E(igraph3))) # make the categories of x into numeric values for color ramp
E(igraph3)$color <- rbPal(10)[cut(as.numeric(E(igraph3)$xx),breaks = 10)]
```
Plotting the igraph vector using the host information... However, I can combine the country and host information in a single plot. 


```{r echo = TRUE}
plot(igraph3,  vertex.shape="star", vertex.norays=rep(2:13), edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
legend(x=-1.7, y=0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2, cex=.8, bty="n", ncol=1)
```




## Clustering 

#Phylo
```{r echo = TRUE}
E(igraph)$weight1 <-  betweenness(igraph, V(igraph), directed = FALSE, weights = NULL, nobigint = TRUE)
E(igraph)$weight2 <- edge_betweenness(igraph, e = E(igraph), directed = FALSE, weights = NULL)
ceb <- cluster_edge_betweenness(igraph, weights = E(igraph)$weight2, directed = FALSE,
  edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,
  modularity = TRUE, membership = TRUE)
dendPlot(ceb, mode="phylo", cex= 0.2)
plot(ceb, igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_mds )
legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)

```



```{r echo = FALSE}
#Otro plot
#class(ceb)
#length(ceb)   # number of communities
#membership(ceb)
#modularity(ceb) # how modular the graph partitioning is
#crossing(ceb, igraph)   # boolean vector: TRUE for edges across communities

```

```{r echo = TRUE}
clp <- cluster_label_prop(igraph, weights = E(igraph)$weight)
plot(clp, igraph, edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk )
#otro 
```

#Hierarchichal clustering
```{r echo = TRUE}
E(igraph)$weight1 <-  betweenness(igraph, V(igraph), directed = FALSE, weights = NULL, nobigint = TRUE)
E(igraph)$weight2 <- edge_betweenness(igraph, e = E(igraph), directed = FALSE, weights = NULL)
ceb <- cluster_edge_betweenness(igraph, weights = E(igraph)$weight2, directed = FALSE,
  edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,
  modularity = TRUE, membership = TRUE)
dendPlot(ceb, mode="hclust", cex= 0.2)
plot(ceb, igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk )
legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)
```


Plotting with color and shape attributes
```{r echo = TRUE}
#igraph3<-graph_from_data_frame(x, directed = FALSE, vertices=DsMVmd)
#V(igraph3)$size=5
#V(igraph3)$xx <- as.numeric(as.factor(V(igraph3)$Country)) # make the categories of x into numeric values for color ramp
#V(igraph3)[V(igraph3)$Host]$Host  <- "square"
#V(igraph3)$color <- counPal(10)[cut(as.numeric(V(igraph3)$xx),breaks = 10)]
#CounColor <- unique(cbind(V(igraph3)$Host, V(igraph3)$color))
#E(igraph3)$xx <- as.numeric(unlist(E(igraph3))) # make the categories of x into numeric values for color ramp
#E(igraph3)$color <- rbPal(10)[cut(as.numeric(E(igraph3)$xx),breaks = 10)]
```
Plotting the igraph vector using the country and host information... 
```{r echo = TRUE}
#plot(igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
#legend(x=-0.9, y=-0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2, cex=.8, bty="n", ncol=1)
```

## SDT plot
```{r echo = TRUE}
DsMVmat2 <- sdt2plot(DsMVmat, diag = 1)
dim(DsMVmat2)
DsMVmat2[1:4,1:3]
```

#Plotting density 
prepraring data for ggplot
```{r echo = TRUE}
DsMVden <- as.data.frame(DsMVmat2[1,])
names(DsMVden)[1]<- "similarity"
ggplot(DsMVden, aes(x = similarity))+ 
  geom_density(fill="tomato", alpha = 0.5)+
  ggtitle("Pairwise Similairity Distribution")+
  theme(plot.title = element_text(hjust=0.5, color="black", size=10, face="bold"),
        panel.background = element_blank())
```

#Ploting heatmap
```{r echo=TRUE}
rSDT <- function(x, dist.dna = TRUE, model = "", size =""){
  if(dist.dna == TRUE) { 
    x <- dist.dna(as.DNAbin(x, "model", model, as.matrix = TRUE)) 
    x= round(as.matrix(-(x-(1))*100), digits=2)
    #x[ x == 0] <- 100
    diag(x) <-100
    #(x) #esto se puede borrar o mostrar un avance con print
    x[lower.tri(x)]<- NA
    #  y <- get_upper_tri(x)
    x <- melt(x, na.rm= TRUE)
  }
  else {
    #x = DsMVmat1
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
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = size , hjust = 1),
          axis.text.y = element_text(angle = 0, vjust = 1, size = size , hjust = 1),
          panel.grid.major = element_blank(), #look for a way to turn off or on the grids 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    coord_fixed()
  } 
```

The function run as follows.
```{r echo=TRUE}
rSDT(DsMVmat2, dist.dna = FALSE, size = 4)
```



```{r echo=TRUE}


```



