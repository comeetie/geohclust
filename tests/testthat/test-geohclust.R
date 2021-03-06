test_that("simple ward", {
  # Fix the seed of the random number generator in order
  # to have reproducible results.
  # cf https://link.springer.com/article/10.1007/s00357-014-9161-z for the repex
  set.seed(19037561)
  # Create the input matrix to be used.
  X <- matrix(runif(20*4),nrow=20,ncol=4)
  N=nrow(X)
  nb=lapply(1:N,\(i){setdiff(1:N,i)})
  geohc=geohclust::geohclust_graph(nb,X,method = "ward")
  hc=hclust(0.5*dist(X)^2,method="ward.D")
  testthat::expect_equal(hc$merge,geohc$merge)
  testthat::expect_equal(cumsum(hc$height),geohc$height,tolerance = 10^-6)
})


test_that("simple centroid", {
  # Fix the seed of the random number generator in order
  # to have reproducible results.
  # cf https://link.springer.com/article/10.1007/s00357-014-9161-z for the repex
  set.seed(19037561)
  # Create the input matrix to be used.
  X <- matrix(runif(20*4),nrow=20,ncol=4)
  N=nrow(X)
  nb=lapply(1:N,\(i){setdiff(1:N,i)})
  geohc=geohclust::geohclust_graph(nb,X,method = "centroid")
  hc=hclust(dist(X)^2,method="centroid")
  testthat::expect_equal(hc$merge,geohc$merge)
  testthat::expect_equal(cumsum(hc$height),geohc$height,tolerance = 10^-6)
})

test_that("simple median", {
  # Fix the seed of the random number generator in order
  # to have reproducible results.
  # cf https://link.springer.com/article/10.1007/s00357-014-9161-z for the repex
  set.seed(19037561)
  # Create the input matrix to be used.
  X <- matrix(runif(20*4),nrow=20,ncol=4)
  N=nrow(X)
  nb=lapply(1:N,\(i){setdiff(1:N,i)})
  geohc=geohclust::geohclust_graph(nb,X,method = "median")
  hc=hclust(dist(X)^2,method="median")
  testthat::expect_equal(hc$merge,geohc$merge)
  testthat::expect_equal(cumsum(hc$height),geohc$height,tolerance = 10^-6)
})

test_that("simple ward zscore", {
  # Fix the seed of the random number generator in order
  # to have reproducible results.
  # cf https://link.springer.com/article/10.1007/s00357-014-9161-z for the repex
  set.seed(19037561)
  # Create the input matrix to be used.
  X <- matrix(runif(20*4),nrow=20,ncol=4)
  Xc <- apply(X,2,\(col){(col-mean(col))/sd(col)})
  N=nrow(X)
  nb=lapply(1:N,\(i){setdiff(1:N,i)})
  geohc=geohclust::geohclust_graph(nb,X,method = "ward",scaling = "zscore")
  hc=hclust(0.5*dist(Xc)^2,method="ward.D")
  testthat::expect_equal(hc$merge,geohc$merge)
  testthat::expect_equal(cumsum(hc$height),geohc$height,tolerance = 10^-6)
})


test_that("ward with constraints", {
  gr=sf::st_make_grid(sf::st_polygon(list(matrix(c(0,100,100,0,0,0,0,100,100,0),ncol=2))),cellsize = 10)
  nb=sf::st_intersects(gr)
  class(nb)="list"
  set.seed(1234)
  X = rbind(cbind(rnorm(50)+5,rnorm(50)+5),cbind(rnorm(50)-5,rnorm(50)-5))
  df.sf = sf::st_sf(geometry=gr,data.frame(X))
  geohc=geohclust::geohclust_graph(nb,X,method = "ward")
  testthat::expect_equal(cutree(geohc,2),rep(1:2,each=50))
  X[1,] = c(-5,-5) 
  geohc=geohclust::geohclust_graph(nb,X,method = "ward")
  testthat::expect_equal(cutree(geohc,3),c(1,rep(2:3,each=49),3))
})


test_that("ward polygons", {
  gr=sf::st_make_grid(sf::st_polygon(list(matrix(c(0,100,100,0,0,0,0,100,100,0),ncol=2))),cellsize = 10)
  nb=sf::st_intersects(gr)
  class(nb)="list"
  set.seed(1234)
  X = rbind(cbind(rnorm(50)+5,rnorm(50)+5),cbind(rnorm(50)-5,rnorm(50)-5))
  df.sf = sf::st_sf(geometry=gr,data.frame(X))
  geohc=geohclust::geohclust_poly(df.sf,method = "ward")
  cl=rep(1:2,each=50)
  names(cl)=1:100
  testthat::expect_equal(cutree(geohc,2),cl)
  geoagg = geocutree(geohc,2)
  testthat::expect_equal(nrow(geoagg),2)
  testthat::expect_equal(sf::st_equals(geoagg$geometry[1],sf::st_union(df.sf[cutree(geohc,2)==1,]))[[1]],1)
  geoaggX=as.matrix(geoagg[,-1] |> sf::st_drop_geometry())
  cm=colMeans(X[cutree(geohc,2)==1,])
  names(cm)=colnames(geoaggX)
  testthat::expect_equal(geoaggX[1,],cm)
  df.sf[1,1:2] = c(-5,-5) 
  geohc=geohclust::geohclust_poly(df.sf,method = "ward")
  cl=c(1,rep(2:3,each=49),3)
  names(cl)=1:100
  testthat::expect_equal(cutree(geohc,3),cl)
  geoagg = geocutree(geohc,3)
  testthat::expect_equal(nrow(geoagg),3)
  testthat::expect_equal(sf::st_equals(geoagg$geometry[1],sf::st_union(df.sf[cutree(geohc,3)==1,]))[[1]],1)
  geoaggX=as.matrix(geoagg[,-1] |> sf::st_drop_geometry())
  cm=colMeans(X[cutree(geohc,3)==2,])
  names(cm)=colnames(geoaggX)
  testthat::expect_equal(geoaggX[2,],cm)
})



test_that("ward polygons queen/rook", {
  gr=sf::st_make_grid(sf::st_polygon(list(matrix(c(0,100,100,0,0,0,0,100,100,0),ncol=2))),cellsize = 10)
  set.seed(1234)
  X = rbind(cbind(rnorm(50)+5,rnorm(50)+5),cbind(rnorm(50)-5,rnorm(50)-5))
  X[50,1:2] = c(-5,-5) 
  X[60,1:2] = c(5,5)
  df.sf = sf::st_sf(geometry=gr,data.frame(X))


  geohc=geohclust::geohclust_poly(df.sf,method = "ward",adjacency = "queen")
  cl=rep(1:2,each=50)
  cl[c(50,60)]=c(2,1)
  names(cl)=1:100
  testthat::expect_equal(cutree(geohc,2),cl)
  geoagg = geocutree(geohc,2)
  testthat::expect_equal(nrow(geoagg),2)
  testthat::expect_equal(sf::st_equals(geoagg$geometry[1],sf::st_union(df.sf[cutree(geohc,2)==1,]))[[1]],1)
  geoaggX=as.matrix(geoagg[,-1] |> sf::st_drop_geometry())
  cm=colMeans(X[cutree(geohc,2)==1,])
  names(cm)=colnames(geoaggX)
  testthat::expect_equal(geoaggX[1,],cm)
  
  
  geohc=geohclust::geohclust_poly(df.sf,method = "ward",adjacency = "rook")
  cl=rep(1:2,each=50)
  cl[60]=1
  names(cl)=1:100
  testthat::expect_equal(cutree(geohc,2),cl)
  geoagg = geocutree(geohc,2)
  testthat::expect_equal(nrow(geoagg),2)
  testthat::expect_equal(sf::st_equals(geoagg$geometry[1],sf::st_union(df.sf[cutree(geohc,2)==1,]))[[1]],1)
  geoaggX=as.matrix(geoagg[,-1] |> sf::st_drop_geometry())
  cm=colMeans(X[cutree(geohc,2)==1,])
  names(cm)=colnames(geoaggX)
  testthat::expect_equal(geoaggX[1,],cm)
})



