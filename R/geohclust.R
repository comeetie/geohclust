#' @title Hierarchical clustering with contiguity constraints between polygons
#'
#' @description This function take an \code{\link{sf:sf}} data.frame and performs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link{sf:sf}} data.frame with polygons like features
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @return an \code{\link{hclust::hclust}} like object with three additional slots
#' \describe{
#'   \item{leafs_geometry}{geometries of the dendrogram leafs as an sfc list}
#'   \item{geotree}{geometries of the dendrogram no-leafs node as an sfc list}
#'   \item{data}{The numeric data (eventually scaled used for the clustering)}
#' }
#' @export
geohclust_poly=function(df,method="ward",scaling="raw"){
  
  if(!is(df,"sf")){
    stop("The dataset must be an sf data.frame.",call. = FALSE)
  }
  
  if(!all(sapply(sf::st_geometry(df),function(u){sf::st_is(u,"MULTIPOLYGON")}) | 
          sapply(sf::st_geometry(df),function(u){sf::st_is(u,"POLYGON")}))){
    stop("The dataset must contains only POLYGONS or MULTIPOLYGONS.",call. = FALSE)
  }
  df_nogeo=sf::st_drop_geometry(df)

  
  # build graph
  nb=sf::st_intersects(df,df)
  class(nb)="list"
  hc_res=geohclust_graph(nb,df_nogeo,method,scaling)
  hc_res$call=sys.call()
  # add geographical data
  hc_res$leafs_geometry = st_geometry(df)
  hc_res$geotree = build_geotree(hc_res$merge,df)
  class(hc_res)=c(class(hc_res),"geohclust")
  hc_res
}




#' @title Hierarchical clustering with contiguity constraints
#'
#' @description This function take an data.frame and performs hierarchical clustering with contiguity 
#' constraints using a graph describing the contiguity (provided )
#' @param adjacencies_list graph describing the contiguity between the rows of df as a list of adjacencies 
#' @param df a data.frame with numeric columns
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore or raw (i.e. no scaling, the default)
#' @return an \code{\link{hclust::hclust}} like object
#' @export
geohclust_graph = function(adjacencies_list,df,method="ward",scaling="raw"){
  
  
  if(!(method %in% c("ward","average","median"))){
    stop("The method argument must be ward or average")
  }
  
  if(!(scaling %in% c("zscore","raw"))){
    stop("The scaling argument must be zscore or raw")
  }
  # select only numeric features
  num_feats = unlist(lapply(df,is.numeric))
  if(sum(num_feats)!=ncol(df)){
    warning("Some features were not numeric and have been removed from the clustering.",call. = FALSE)
    df=df[,num_feats]
  }
  # check for missing values
  if(sum(is.na(df))>0){
    stop("Some regions have missing values and missing values are not allowed.",call. = FALSE)
  }
  # scales
  if(scaling=="zscore"){
    df_scaled = t(t(df)-colMeans(df))
    df_scaled = t(t(df_scaled)/sqrt(colMeans(df_scaled^2)))
  }else{
    df_scaled = df
  }
  
  nb_c = lapply(adjacencies_list,\(nei){nei-1})
  # run the algorithm
  res=hclustcc_cpp(nb_c,as.matrix(df_scaled),method)
  
  # complete merge tree if constraints does not allow to reach K=1
  nbmissmerge = sum(rowSums(res$merge==0)==2)
  if(nbmissmerge>0){
    warning("Some regions were isolated. The hierarchy was automatically completed to reach one cluster.",call. = FALSE)
    missing_merges = setdiff(c(-1:-nrow(df_scaled),1:(nrow(df_scaled)-1-nbmissmerge)),c(res$merge[,1],res$merge[,2]))
    for(i in (nrow(df_scaled)-nbmissmerge):(nrow(df_scaled)-1)){
      cmerge = missing_merges[1:2]
      res$merge[i,]=cmerge
      res$height[i]=res$height[nrow(df_scaled)-nbmissmerge-1]*1.2
      missing_merges=c(missing_merges[3:length(missing_merges)],i)
    }
  }
  
  
  # format the results in hclust form
  hc_res = list(merge=res$merge, 
                height=cumsum(res$height),
                order=order_tree(res$merge,nrow(res$merge)),
                labels=(rownames(df)),
                call=sys.call(),
                method=method,
                dist.method="euclidean",
                data=df_scaled)
  class(hc_res)  <- "hclust"

  hc_res
}


#' @title Cut a Geograpĥic Tree into Groups of Data and return an sf data.frame 
#'
#' @description Cuts a tree, e.g., as resulting from geohclust_poly, into several groups either by specifying the desired number(s) of groups or the cut height(s).
#' @param tree a tree as produced by geohclust. cutree() only expects a list with components merge, height, and labels, of appropriate content each.
#' @param k an integer scalar or vector with the desired number of groups
#' @param h numeric scalar or vector with heights where the tree should be cut.
#' At least one of k or h must be specified, k overrides h if both are given.
#' @return an \code{\link{sf::sf}} like object
#' @export
geocutree=function(tree,k = NULL, h= NULL){
  if(!is(tree,"geohclust")){
    stop("geocutree only accepts geohclust objects.")
  }
  if(is.null(k) && is.null(h)){
    stop("At least one of k or h must be specified, k overrides h if both are given.")
  }
  if(!is.null(h) & is.null(k)){
    k = which(tree$height>h)[1]
  }
  cl = cutree(tree,k=k)
  Xg = aggregate(tree$data,list(cl),mean)
  N=nrow(tree$merge)+1
  istart = which(!duplicated(cl))
  clust_geo = list()
  ck=1
  for (i in istart){
    f = -i
    cnode = -i
    while(length(f)!=0){
      f = which(tree$merge[1:(N-k),1]==f | tree$merge[1:(N-k),2]==f)
      if(length(f)>0){
        cnode = f
      }
    }
    if(cnode>0){
      clust_geo[[ck]]=tree$geotree[[cnode]]
    }else{
      clust_geo[[ck]]=tree$leafs_geometry[[-cnode]]
    }
    ck=ck+1
  }
  st_sf(cl=1:k,Xg[,-1],geometry=sf::st_as_sfc(clust_geo,crs=st_crs(tree$leafs_geometry)))
}




order_tree=function(merge,i){
  if(merge[i,1]<0){
    left = -merge[i,1];
  }else{
    left = order_tree(merge,merge[i,1])
  }
  if(merge[i,2]<0){
    right = -merge[i,2];
  }else{
    right = order_tree(merge,merge[i,2])
  }
  c(left,right)
}


build_geotree=function(merge,df){
  geoms= sf::st_geometry(df)
  geotree = list()
  for (i in 1:nrow(merge)){
    if(merge[i,1]<0){
      left = geoms[[-merge[i,1]]];
    }else{
      left = geotree[[merge[i,1]]]
    }
    if(merge[i,2]<0){
      right = geoms[[-merge[i,2]]];
    }else{
      right = geotree[[merge[i,2]]]
    }
    geotree[[i]] = sf::st_union(left,right)
  }
  geotree
}


