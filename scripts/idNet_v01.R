# R Alcala                      Nov 2019
# 
#
#
# 
# R version 3.5.1 (2018-07-02)
#
setwd('~/Documents/Documents-ehecatl/git_db/id_net/')
#---------------------- loading matrix ---------------------- 

psm <-  read.csv("Begomovirus_ex_id_net-2019.csv", header = T, row.names = 1, stringsAsFactors = F)
colnames(psm) <- rownames(psm)
m=psm
m[upper.tri(m)] <- t(m)[upper.tri(m)] 
psm=as.matrix(m) # matrix
psm.ss <- (psm >= 94)*1
psm.trs <- psm * psm.ss
# psm.l <- psm.trs * seq.lenght
#-------------
psg <- graph_from_adjacency_matrix(psm.trs, weighted = T, diag = F)

#--- network metrics
# weighted strength
psg.str <- strength(psg, vids = V(psg), mode = c("all"), loops = F)
psg.wgt <- psg.str/max(psg.str)


#--- palettes
# nodes
# nbPal <- colorRampPalette(c("yellow", "red", "purple")) # volcano type
# links
rbPal <- colorRampPalette(c("blue", "yellow", "red")) 
#--- attributes
# nodes
V(psg)$color <- "gold" # or nbPal(10)[cut(as.numeric(V(psg)$xx),breaks = 10)]
V(psg)$size <- 13*psg.wgt
# links
E(psg)$xx <- as.numeric(unlist(psm.trs)) # make the categories of x into numeric values for color ramp
E(psg)$color <- rbPal(10)[cut(as.numeric(E(psg)$xx),breaks = 10)]
E(psg)$width <- (E(psg)$xx)[as.numeric(E(psg)$xx) ifelse esto.es.igual.o.mayor.que.0 >= 0 ]
E(psg)$width 
# cairo_ps("distances.eps", width = 7.2, height = 7.2,fallback_resolution = 1200, antialias = "subpixel")
plot(psg,  edge.arrow.size= F, vertex.label.cex= (V(psg)$size)/4, vertex.label.color= 'black',
          edge.curved= F, edge.width= E(psg)$width, layout=layout_with_kk)
# dev.off()


