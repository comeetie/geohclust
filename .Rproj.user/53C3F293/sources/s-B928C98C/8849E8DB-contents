#' @title Hierarchical clustering with contiguity constraints between polygons
#'
#' @description This function take an \code{\link{sf:sf}} data.frame and perfomrs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link{sf:sf}} data.frame with polygons like features
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @return an \code{\link{hclust::hclust}} like object
#' @export
geohclust_poly=function(df,method="ward",scaling="raw"){
  
  if(!is(df,"sf")){
    stop("The dataset must be an sf data.frame.",call. = FALSE)
  }
  
  if(!all(sf::st_is(sf::st_geometry(df),"MULTIPOLYGON")||sf::st_is(sf::st_geometry(df),"POLYGON"))){
    stop("The dataset must contains only POLYGONS or MLTIPOLYGONS.",call. = FALSE)
  }
  df_nogeo=sf::st_drop_geometry(df)

  
  # build graph
  nb=sf::st_intersects(df,df)
  class(nb)="list"
  geohclust_graph(nb,df_nogeo,method,scaling)
}




#' @title Hierarchical clustering with contiguity constraints
#'
#' @description This function take an data.frame and performs hierarchical clustering with contiguity 
#' constraints using a graph describing the contiguity (provided )
#' @param adjacencies_list graph describing the coniguities between the rows of df as a list of adjacency 
#' @param df a data.frame with numeric columns
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
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
    warning("Some regions were isolated the hierarchy was automatically completed to reach one cluster.",call. = FALSE)
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
                order=1:nrow(df),
                labels=(rownames(df)),
                call=sys.call(),
                method=method,
                dist.method="euclidean")
  class(hc_res)  <- "hclust"
  hc_res=as.hclust(reorder(as.dendrogram(hc_res),1:nrow(df)))
  hc_res$method=method
  hc_res$dist.method="euclidean"
  hc_res$call = sys.call()
  hc_res
}
order_tree = function(tree){
  cnodes = c(length(tree))
  visited = c()
  L=0.5
  x= rep(0,length(tree))
  while(length(cnodes)>0){
    newnodes = c()
    for (f in cnodes){
      visited=c(visited,f)
      children=which(tree==f)
      x[children[1]]=x[f]-L
      x[children[2]]=x[f]+L
      newnodes=c(newnodes,children)
    }
    cnodes=setdiff(newnodes,visited)
    L=L/2
  }
  order(x[1:((length(x)+1)/2)]) 
}

