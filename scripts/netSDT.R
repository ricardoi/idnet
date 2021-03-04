library(igraph)
library(viridis)
library(tidyr)
library(dplyr)
library(reshape2)
library(Cairo)

setwd("/Users/RicardoI/++Manuscripts/++2017/rSDT/")
DsMV.distances <- read.csv("DsMV_pairwisematrix_codon.csv", stringsAsFactors = F, header = TRUE)
blah <- read.csv("DsMV_pairwisematrix_codon.csv", stringsAsFactors = F, header = TRUE)
blah<-blah[,2:1]
dim(DsMV.distances)
class(DsMV.distances)
mymat1 <- as.matrix(DsMV.distances[, 3:ncol(DsMV.distances)])
rownames(mymat1)<- DsMV.distances$X.1
dim(mymat1)
class(mymat1)
diag(mymat1)<-0
mymat <- (mymat1 > 88)*1
mymat[1:5,1:5]
mymat2 <- graph_from_adjacency_matrix(mymat)
#mymat2 <- graph_from_incidence_matrix(distances1) #look at this graph later
rbPal <- colorRampPalette(c("blue", "yellow", "red")) 
V(mymat2)$color <- ifelse(V(mymat2)$xx > 62, "red", "blue")
V(mymat2)$size=5
E(mymat2)$xx<-as.numeric(unlist(distances1)) # make the categories of x into numeric values for color ramp
E(mymat2)$color <- rbPal(10)[cut(as.numeric(E(mymat2)$xx),breaks = 10)]


cairo_ps("distances.eps", width = 7.2, height = 7.2,fallback_resolution = 1200, antialias = "subpixel")
plot(mymat2,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black' ,edge.curved=T, edge.width=0.2, layout=layout_with_kk)
dev.off()

plot(Papayas.network, vertex.size=10, edge.arrow.size=0.2, vertex.frame.color="black",edge.curved=F,edge.width=1, vertex.label.cex=0.5,
     vertex.label.color="blue",vertex.label.family="Times New Roman",vertex.label.font=3, layout=layout.fruchterman.reingold, main="Co-occurrence of viruses in Papaya symptoms by sample")


#################
#MELT
DsMV.dis = DsMV.distances
Country.id <- DsMV.dis[,2:1]
blah<-blah[,2:1]

#diag(DsMV.dis[,3:108])<-0
#mymat3 <- melt(DsMV.dis)
#names(mymat3)
#head(mymat3)
#mymat3<-cbind(From=as.character(mymat3$ID), To=as.character(mymat3$variable), Country=as.character(mymat3$Country), value=as.character(mymat3$value))
#mymat3<-as.data.frame(mymat3)
#write.csv(mymat3, "mat3.csv")
mymat3<- read.csv("mat3.csv", stringsAsFactors = F, header = TRUE)
match(unique(mymat3$To), unique(blah$ID))
cbind(unique(mymat3$To), unique(blah$ID))

#mymat3 <- melt(mymat1)
mymat4 <- mymat3[!rowSums(mymat3[-c(1:3)] == 0) >= 1,] #removing 0's
mymat5 <- mymat3[mymat3[-c(1:3)] > 86,] # treshold
#mymat5 <- mymat4[(mymat4$value > 99)*1]
mymat <- graph_from_edgelist(as.matrix(mymat5[,2:3]), directed = F)
#mymat <- simplify(mymat, remove.multiple = TRUE, remove.loops = TRUE)#removing something
#mymat.u <- as.undirected(mymat.s)
# colour and subsets
rbPal <- colorRampPalette(c("green", "yellow", "red")) 
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
class(mymat5)
igraph1<-graph_from_data_frame(mymat5, directed = FALSE, vertices=blah)
#subseting for colors 
#V(mymat)$color <- ifelse(V(mymat)$xx > 88, "red", "blue")
#V(igraph1)$xx <- as.numeric(V(igraph1)$Country)
#V(mymat)$xx <- as.numeric(unlist(Country.id$Country[V(mymat)]))
#V(igraph1)$color <- counPal(7)[V(igraph1)$xx]
#V(mymat)$color <- counPal(10)[Country.id$Country[V(mymat)]]
V(igraph1)$size=5
V(igraph1)$xx <- as.numeric(as.factor(V(igraph1)$Country)) # make the categories of x into numeric values for color ramp
V(igraph1)$color <- counPal(10)[cut(as.numeric(V(igraph1)$xx),breaks = 10)]
#unique(cbind(V(igraph1)$Country, V(igraph1)$color))
E(igraph1)$xx <- as.numeric(unlist(E(igraph1))) # make the categories of x into numeric values for color ramp
E(igraph1)$color <- rbPal(10)[cut(as.numeric(E(igraph1)$xx),breaks = 10)]
#plot igraph

cairo_ps("DsMVdistances86.eps", width = 7.2, height = 7.2,fallback_resolution = 1200, antialias = "subpixel")
plot(igraph1,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
legend(x=-1.5, y=-1.1, rownames(table(Country.id$Country[V(igraph1)])), pch=21,  col="black", pt.bg=c("#FF0000",  "#FF7100", "#AAAA54", "#3838C6", "#A9A9FF", "#CD8888", "#A52A2A"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


##koMV
KoMV <- read.csv("Ko_Z_JHMV_final00_ord_matrix.csv", row.names=1, header = TRUE)
KoMV.mat <-as.matrix(KoMV)
diag(KoMV.mat)<-0
KoMV.m <- melt(KoMV.mat)
KoMV.m0 <- KoMV.m[!rowSums(KoMV.m[-c(1:2)] == 0) >= 1,] #removing 0's
KoMV.n <- KoMV.m[KoMV.m[-c(1:2)] > 97,] # treshold
#mymat5 <- mymat4[(mymat4$value > 99)*1]
KoMV.net <- graph_from_edgelist(as.matrix(KoMV.n[,1:2]), directed = F)
rbPal <- colorRampPalette(c("green", "yellow", "red")) 
V(KoMV.net)$color <- ifelse(V(KoMV.net)$xx > 93, "red", "blue")
V(KoMV.net)$size=5
E(KoMV.net)$xx <- as.numeric(unlist(E(KoMV.net))) # make the categories of x into numeric values for color ramp
E(KoMV.net)$color <- rbPal(10)[cut(as.numeric(E(KoMV.net)$xx),breaks = 10)]
V(KoMV.net)
#categories <- unique(mymat4$value)
cairo_ps("KoMVdistances97.eps", width = 7.2, height = 7.2,fallback_resolution = 1200, antialias = "subpixel")
plot(KoMV.net,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black' ,edge.curved=T, edge.width=0.2, layout=layout_with_kk)
dev.off()
#### OLD CODE


#as.numeric(unlist(distances1))
mymat2<-graph_from_adjacency_matrix(mymat)
rbPal <- colorRampPalette(c("blue", "yellow", "red")) 
E(mymat2)$xx<-as.numeric(unlist(distances1)) # make the categories of x into numeric values for color ramp

# V(mymat2)$size=5

#rbPal <- viridis_pal(option="C", direction = 1)  #This adds a color ramp of colors
rbPal <- colorRampPalette(c("blue", "yellow", "red")) 
#as.numeric(unlist(a))
# V(mymat2)$xx<-as.numeric(unlist(distances1))  # make the categories of x into numeric values for color ramp
#V(mymat2)$xx<-as.numeric(V(mymat2)$xx)
#V(mymat2)$color <- ifelse(V(mymat2)$xx>90, "red", "blue")

#V(mymat2)$color <- rbPal(20)[cut( as.numeric(V(mymat2)$xx),breaks = 20)]
#E(mymat2)$color <- rbPal(20)[cut( as.numeric(E(mymat2)$value),breaks = 20)]

E(mymat2)$xx<-as.numeric(unlist(distances1)) # make the categories of x into numeric values for color ramp
#V(mymat2)$xx<-as.numeric(V(mymat2)$xx)
V(mymat2)$color <- ifelse(V(mymat2)$xx>88, "red", "blue")
V(mymat2)$size=4
E(mymat2)$color <- rbPal(10)[cut(as.numeric(E(mymat2)$xx),breaks = 10)]

plot(mymat2, edge.arrow.size=.1, vertex.label.cex=.3, vertex.label.color='black' ,edge.curved=T, layout=layout_with_kk)

cairo_ps("distances.eps", width = 7.2, height = 7.2,fallback_resolution = 1200, antialias = "subpixel")

plot(mymat2, edge.arrow.size=.1, vertex.label.cex=.3, vertex.label.color='black' ,edge.curved=T, layout=layout_with_kk)
dev.off()


#### new data
DsMV17 <- read.csv("DsMV_pairwisematrix_codon_2017.csv", stringsAsFactors = F, header = TRUE)
DsMV.dis = DsMV17 
diag(DsMV.dis[,3:129])<-0
CountryDsMV.id<-DsMV.dis[,2:1]
Country <- "Country"
ID <- "ID"
colnames(CountryDsMV.id)<-c(ID, Country)
head(CountryD.id)
Names <- DsMV.dis$Names
colnames(DsMV.dis) <- c(Country, ID, Names)
head(DsMV.dis)
tail(DsMV.dis)
DsMV.m <- melt(DsMV.dis)
head(DsMV.m)
DsMV.m<-cbind(From=as.character(DsMV.m$ID), To=as.character(DsMV.m$variable), Country=as.character(DsMV.m$Country), value=as.character(DsMV.m$value))
DsMV.m<-as.data.frame(DsMV.m)
write.csv(DsMV.m, "DsMV.m2017")
# TO LOAD
DsMV.m2017 <- read.csv("DsMV.m2017.csv", stringsAsFactors = F, header = TRUE)
DsMV.dat <- DsMV.m2017
head(DsMV.dat)
class(DsMV.dat)
#match(unique(DsMV.dat$To), unique(CountryDsMV.id$ID))
#cbind(unique(DsMV.dat$To), unique(CountryDsMV.id$ID))
#mymat3 <- melt(mymat1)
DsMV.dat <- DsMV.dat[!rowSums(DsMV.dat[-c(1:3)] == 0) >= 1,] #removing 0's
DsMV.dat <- DsMV.dat[DsMV.dat[-c(1:3)] > 85,] # treshold
#mymat5 <- mymat4[(mymat4$value > 99)*1]
#mymat <- graph_from_edgelist(as.matrix(mymat5[,2:3]), directed = F)
#mymat <- simplify(mymat, remove.multiple = TRUE, remove.loops = TRUE)#removing something
#mymat.u <- as.undirected(mymat.s)
# colour and subsets
rbPal <- colorRampPalette(c("green", "yellow", "red")) 
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
class(DsMV.dat)
igraph1<-graph_from_data_frame(DsMV.dat, directed = FALSE, vertices=CountryDsMV.id)
#subseting for colors 
#V(mymat)$color <- ifelse(V(mymat)$xx > 88, "red", "blue")
#V(igraph1)$xx <- as.numeric(V(igraph1)$Country)
#V(mymat)$xx <- as.numeric(unlist(Country.id$Country[V(mymat)]))
#V(igraph1)$color <- counPal(7)[V(igraph1)$xx]
#V(mymat)$color <- counPal(10)[Country.id$Country[V(mymat)]]
V(igraph1)$size=5
V(igraph1)$xx <- as.numeric(as.factor(V(igraph1)$Country)) # make the categories of x into numeric values for color ramp
V(igraph1)$color <- counPal(10)[cut(as.numeric(V(igraph1)$xx),breaks = 10)]
CounColor <- unique(cbind(V(igraph1)$Country, V(igraph1)$color))
E(igraph1)$xx <- as.numeric(unlist(E(igraph1))) # make the categories of x into numeric values for color ramp
E(igraph1)$color <- rbPal(10)[cut(as.numeric(E(igraph1)$xx),breaks = 10)]
#plot igraph

cairo_ps("DsMV2017_85.eps", width = 7.2, height = 7.2,fallback_resolution = 1200, antialias = "subpixel")
plot(igraph1,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
legend(x=-0.9, y=-0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()
