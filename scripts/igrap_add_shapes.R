# all vertex shapes, minus "raster", that might not be available
shapes <- setdiff(shapes(), "")
shapes
g <- make_ring(length(shapes))
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
add_shape("romboid", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=2))
add_shape("triangle", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=3))
add_shape("star4", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=4))
add_shape("star5", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))
add_shape("star6", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=6))
add_shape("star8", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=8))
add_shape("star10", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=10))
add_shape("star12", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=12))

plot(g, vertex.shape="star4", vertex.color=rainbow(vcount(g)),
     vertex.size=seq(10,20,length=vcount(g)))
plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
     vertex.size=seq(10,20,length=vcount(g)),
     vertex.norays=rep(c(3:8, 10), length=vcount(g)))

#### SECOND TRIAL: TRIANGLE
#################################################################
# generic star vertex shape, with a parameter for number of rays
mytriangle <- function(coords, v=NULL, params) {
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

#### SECOND TRIAL: TRIANGLE
#################################################################
# generic star vertex shape, with a parameter for number of rays
myromboid <- function(coords, v=NULL, params) {
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