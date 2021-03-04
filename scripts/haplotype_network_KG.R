#Haplotypes
source("http://www.bioconductor.org/biocLite.R")
biocLite("msa")
biocLite("Biostrings")
library(msa)
library(Biostrings)
library(haplotypes)
library(seqinr)
library(igraph)
library(Cairo)

#requires Biostrings
seq <- read.fasta("Documents/Documents - ehecatl/git_db/id_net/DsMV_2017_final.fasta", as.string = FALSE, seqtype = "DNA", seqonly=FALSE, strip.desc= TRUE)
seqstr <- read.fasta("Documents/Documents - ehecatl/git_db/id_net/DsMV_2017_final.fasta", as.string = FALSE, seqtype = "DNA", seqonly=TRUE, strip.desc= FALSE)

myseq <- getSequence(seq)
myseq <- matrix(unlist(myseq), ncol = length(s2c(seqstr[[1]])) ,byrow = TRUE)
rownames(myseq)<- names(seq)
myseq[1:12,1:10]
x<-as.dna(myseq) 
#x <- x[1:40, as.matrix=FALSE]
h<-haplotype(x) #inferring haplotypes
h
as.dna(h)


# giving DNA sequences of haplotypes
d<-distance(x) #inferring haplotypes dist obj
d #to use for igraph
class(d)
#View(as.matrix(d))
dh<-haplotype(d) 
dh

# pie charts Network
xm <- as.matrix(d)
class(xm)
mygraph <- graph.adjacency(xm)
#plot(mygraph, edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black', edge.curved=T, edge.width=0.2, layout=layout_with_mds)

gvalues <- as.list(1:nrow(xm))
#gvalues <- 
for(i in 1:nrow(xm)){
       gvalues[[i]] <- xm[,i]
       }
gvalues
cairo_ps("haplo_graph.eps", width = 7.2, height = 7.2, fallback_resolution = 1200, antialias = "subpixel")
#png("haplo_graph.png")
plot(mygraph, vertex.shape="pie", vertex.pie=gvalues, vertex.pie.color=list(rainbow(5)), edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black', edge.curved=T, edge.width=0.2, layout=layout_with_mds)
#dev.off()




## Coercing a matrix to a 'Dna' object.
# all valid characters
x<-matrix(c("?","A","C","g","t","-","0","1","2","3","4","5"),4,6)
rownames(x)<-c("seq1","seq2","seq3","seq4") 
x
dna.obj<-as.dna(x) 
dna.obj
#### EXAMPLE
data("dna.obj")
x<-dna.obj[1:6,as.matrix=FALSE]

##Inferring haplotypes using 'Dna' object.
# coding gaps using simple indel coding method
h<-haplotype(x)
h
# giving DNA sequences of haplotypes
as.dna(h) 
## Inferring haplotypes using dist object. 
d<-distance(x)
d # here is a dist object
as.matrix(d) #here is a matrix object

dh<-haplotype(d) 
dh
summary(dh)
# pie charts in network indicating gene spread




# example from igraph

g <- make_ring(10)
values <- lapply(1:10, function(x) sample(1:10,3))
if (interactive()) {
  plot(g, vertex.shape="pie", vertex.pie=values,
       vertex.pie.color=list(heat.colors(5)),
       vertex.size=seq(10,30,length=10), vertex.label=NA)
}

# gene network example

gmat <- matrix(0,ncol=10,nrow=10)
gmat[3,1] <- 1
gmat[3,2] <- 1
gmat[6,4] <- 1
gmat[6,3] <- 1
gmat[6,5] <- 1
gmat[6,7] <- 1
gmat[7,6] <- 1
gmat[7,8] <- 1
gmat[7,9] <- 1
gmat[8,10] <- 1
gmati <- graph.adjacency(gmat)
plot(gmati)

gvalues <- as.list(1:10)
gvalues[[1]] <- c(0,1,1,0,0)
gvalues[[2]] <- c(0,0,1,1,0)
gvalues[[3]] <- c(0,1,1,1,0)
gvalues[[4]] <- c(1,1,1,0,0)
gvalues[[5]] <- c(1,0,1,0,0)
gvalues[[6]] <- c(1,1,1,1,0)
gvalues[[7]] <- c(0,1,1,1,1)
gvalues[[8]] <- c(0,0,0,1,1)
gvalues[[9]] <- c(0,1,1,1,0)
gvalues[[10]] <- c(0,0,0,0,1)
gvalues
plot(gmati, vertex.shape="pie", vertex.pie=gvalues, vertex.pie.color=list(rainbow(5)), vertex.size=c(10,10,10,10,10,20,20,10,10,10), vertex.label=NA)


