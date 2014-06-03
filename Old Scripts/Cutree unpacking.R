

tree <- td.clust.aver
k <- 5

function (tree, k = NULL, h = NULL) 
{
  if (is.null(n1 <- nrow(tree$merge)) || n1 < 1) 
    stop("invalid 'tree' ('merge' component)")
  n <- n1 + 1
  if (is.null(k) && is.null(h)) 
    stop("either 'k' or 'h' must be specified")
  if (is.null(k)) {
    if (is.unsorted(tree$height)) 
      stop("the 'height' component of 'tree' is not sorted (increasingly)")
    k <- n + 1L - apply(outer(c(tree$height, Inf), h, ">"), 
                        2, which.max)
    if (getOption("verbose")) 
      message("cutree(): k(h) = ", k, domain = NA)
  }
  else {
    k <- as.integer(k)
    if (min(k) < 1 || max(k) > n) 
      stop(gettextf("elements of 'k' must be between 1 and %d", 
                    n), domain = NA)
  }
  ans <- .Call(C_cutree, tree$merge, k)
  if (length(k) == 1L) {
    ans <- setNames(as.vector(ans), tree$labels)
  }
  else {
    colnames(ans) <- if (!is.null(h)) 
      h
    else k
    rownames(ans) <- tree$labels
  }
  return(ans)
}

cutree(tree, k)

group.mem[group.mem$group.num == 1, ] #september 2010 has the 900m sites in group 1. 

tree <- dendrogram.list[["2010S"]]
comm.dissim[["2010Sbsim"]]
cutree(tree, k)
tree$labels
head(site.info)