#' Main function of epvclust.
#' Compute hierchical clustering and calculate AU-values and BP-values for each cluster.
#' 
#' @param filename path he of the file cointaining the dataset (it must be in CSV format)
#' @param nboot number of bootstrap replicate to do for each bootstrap session
#' @param method.hclust clustering linkage method
#' @param transpose if 0 the objects are on columns
#' 
#' @return hierachical clustering with AU-values and BP-values
epvdriver <- function(filename, nboot=10, r = c(50, 75, 100, 125, 150), method.hclust="average", transpose = 0){

  #set variables
  method.dist = "correlation";
  use.cor = "pairwise.complete.obs";
  
  #get method for clustering
  switch(method.hclust,
    average = {
      mt <- "a";  
    },
    centroid = {
      mt <- "c";  
    },
    single = {
      mt <- "s";  
    },
    {
      mt <- "a";
    }
  )

  #call epvclust
  eboot <- epvclust(filename = filename, nboot = nboot, r = r, mt = mt, transpose = as.integer(transpose));
  
  #get variables from eboot
  freqM<- eboot$f;
  rv <- eboot$r;
  r <- rv/100;
  size <- r*eboot$ncols;
  lr <- eboot$lr;
  nboot <- eboot$nboot;
  
  #get pattern (a pattern is how the dendrogram is composed)
  pattern <- eboot$pattern;
  edges.cnt <- table(factor(pattern))-table(factor(pattern));
  
  #build data.hclust
  data <- read.table(filename, header=FALSE, sep = ",");
  distance <- as.dist(cor(data, method="pearson"))
  distance <- edist(data, method=method.dist, use.cor=use.cor);
  data.hclust <- hclust(distance, method=method.hclust);

  
  
  #build mboot
  #mboot is a list of lenght eboot$lr
  mboot <- vector("list", lr);
  for (i in 1:lr){
    edges.cnt <- freqM[i,];
    mboot[[i]] <- list( edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor, method.hclust=method.hclust, 
                        nboot=nboot, size=size[i], r=r[i], store= list());
    class(mboot) <- "boot.hclust";
  }
  #   #call merge from pvclust
  result <- epvclust.merge(data=data, object.hclust=data.hclust, mboot=mboot, pattern = pattern);
  
  return(result);
}



epvclust <-function(filename, nboot = 10, r = c(50, 75, 100, 125, 150), mt ="a", transpose = 0) {
# Call C code for doing pvclust bootstrap routine
  
  lr <- length(r);
  if (lr <= 1){ #no sense to do multiscale bootstrap
    print("Lenght of r too small");
    return;
  }
  
  #calculate l as the total number of elements
  t <- read.table(file = filename, sep = ",");
  nobjs <- ncol(t)-1;  #because objects are on the columns and nrows needs to be the number of clusters (#objects-1)
  l <- nobjs*lr;
  
  lp <- nobjs * (nobjs+1);
  
  res <- .C( "epvclust", filename, nboot = as.integer(nboot),  f = as.integer(1:l), 
             r = as.integer(r), lr = as.integer(lr), pattern = as.integer(1:lp), mt = mt, 
             transpose = as.integer(transpose), PACKAGE = "epvclust");
  
  #turn f from vector into matrix
  res$f <- matrix(res$f, nrow = lr, ncol = nobjs, byrow = TRUE);
  
  #turn pattern from vector into matrix
  p <- res$pattern;
  p <- matrix(p, nrow = nobjs, ncol = nobjs+1, byrow = TRUE);
  res$pattern <- apply(p,1,paste,collapse="")

  res$nrows <- nobjs + 1;
  res$ncols <- nrow(t);
     
  return(res)
}


#---------------------------------------------------------------------------------


edist <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  else if(!is.na(pmatch(method,"uncentered"))){
    if(sum(is.na(x)) > 0){
      x <- na.omit(x)
      warning("Rows including NAs were omitted")
    }
    x  <- as.matrix(x)
    P  <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q  <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)
  }
  else
    dist(t(x),method)
}


epvclust.merge <- function(data, object.hclust, mboot, pattern){
  
#   pattern <- hc2split(object.hclust)$pattern
#   print("pat")
#   print(pattern)
  r     <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  store <- lapply(mboot,"[[", "store")
  
  rl <- length(mboot)
  ne <- length(pattern)
  
  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl))
  row.names(edges.bp) <- pattern
  names(edges.cnt) <- paste("r", 1:rl, sep="")
  
  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) 
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
  }
  
  ms.fitted <- lapply(as.list(1:ne),
                      function(x, edges.bp, r, nboot){
                        msfit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  class(ms.fitted) <- "mslist"
#   
  p    <- lapply(ms.fitted,"[[","p")
  se   <- lapply(ms.fitted,"[[","se")
  coef <- lapply(ms.fitted,"[[","coef")
  
  au    <- unlist(lapply(p,"[[","au"))
  bp    <- unlist(lapply(p,"[[","bp"))
  se.au <- unlist(lapply(se,"[[","au"))
  se.bp <- unlist(lapply(se,"[[","bp"))
  v     <- unlist(lapply(coef,"[[","v"))
  cc    <- unlist(lapply(coef,"[[","c"))
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
  
  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp, v=v, c=cc, pchi=pchi)
  
  row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne
  
  result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
                 msfit=ms.fitted, nboot=nboot, r=r, store=store)
  
  class(result) <- "pvclust"
  return(result)
}
