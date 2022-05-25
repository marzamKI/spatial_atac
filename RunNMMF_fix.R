#' Run Non-negative Matrix Factorization
#'
#' Decompose an expression matrix A with non-negative elements into matrices WxH, also with
#' non-negative elements. W is the feature loading matrix (features x factors) and H is the
#' low dimensional embedding of the spots (factors x spots).
#'
#' @param object Seurat object
#' @param assay Assay Name of Assay NMF is being run on
#' @param features Features to compute the NMF for. Note that these features must be present in the
#' slot used to compute the NMF. By default, the `features` is set to `VariableFeatures(object)`
#' to include the most variable features selected in the normalization step.
#' @param nfactors Total Number of factors to compute and store (20 by default)
#' @param rescale Rescale data to make sure that values of the input matrix are non-n
#' @param reduction.name Dimensional reduction name, "NMF" by default
#' @param reduction.key Dimensional reduction key, specifies the prefix of the factor ids, e.g.
#' "factor_1", "factor_2", etc.
#' @param n.cores Number of threads to use in coom parallel detectCores
#'
#' @export
#'
RunNMF_fix <- function (
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  nfactors = 20,
  rescale = TRUE,or = FALSE,
  sort.spcor.by.var = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  var.genes <- features %||% VariableFeatures(object)
  norm.coun<- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }
  if (min(norm.counts) < 0) stop("Negative values e not allowed")
  nmf.results <- rnmf(A = norm.counts[var.genes, ], k = nfactors)
  #nmf.results$W <- swne::ProjectFeatures(norm.counts, nmf.results$H, n.cores = n.cores)
  feature.loadings <- nmf.results$W
  cell.embeddings <- t(nmf.results$H)

  # Set cores
  n.cores <- n.cores %||% {detectCores() - 1}

  # Order factors based on spatial correlati <- do.call(rbind, GetSpatNet(object = object, nNeighbours = NULL, maxdist = NULL))
    resCN <- as.matrix(data.frame(reshapedcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
    resCN[resCN > 0] <- 1
    eix(0, nrow = nrow(cell.embeddings), ncol = nrow(cell.embeddings), dimnames = list(rownames(cell.embeddings), rdings)))
    colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
n = "^X", replacement = "", x = colnames(resCN))
    empty.CN[rownames(resCN), colnames(resCN)] <- resCN
    empty.CN)
    fun <- function (x) lag.listw(listw, x, TRUE)

    # Calculate the lag matrix from (cell.embeddings, 2, fun)

    # Split sp.cor by sample
    if (sort.spcor.by.var) {
      sp.cor.split (unique(GetStaffli(object)@meta.data$sample), function(s) {
        tablag.split <- tablag[GetStaffli(object)@meta.data$sample =(cell.embeddings.split), function(i) {
          cor(tablag.split[, i], cell.embeddings.split[order.vec <- order(apply(sp.cor.split, 2, var))
    } else {
      sp.cor <- unlist(lapply(1:ncol(c{
        cor(cell.embeddings[, i], tablag[, i])
      }))
      order.vec <- order(sp.cor, decreasing = TRUE)
    }

    cell.embe rownames(x = feature.loadings) <- var.genes
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nfactors)
  rowname= cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.l CreateDimReducObject (
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
n.key
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}


rnmf <- function (
  A,
  kica",
  n.cores = 1,
  loss = "mse",
  max.iter = 500,
  ica.fast = F
) {
  if (any(A < 0))
    stop("The input matrix conegative elements !")
  if (k < 3)
    stop("k must be greater than or equal to 3 to create% c("ica", "nnsvd", "random")) {
    stop("Invalid initialization method")
  }
  A <- as.matrix(A)
  if (any(A < 0)) {
    sto"Input matrix has negative values")
  }
  if (init == "ica") {
    nmf.init <- ica_init(A, k, ica.fast = ica.fas (init == "nnsvd") {
    nmf.init <- nnsvd_init(A, k, LINPACK = T)
  }
  else {
    nmf.init <- NULL
  }
  if (is.nul) {
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, n.threads = n.cores,
                          loss = loss, max.iter = max.iter)
 nit$H < zero.eps] <- 0
    zero.idx.w <- which(nmf.init$W == 0)
    zero.idx.h <- which(nmf.init$H == 0)
    nmf.w] <- runif(length(zero.idx.w), 0,
                                    A.mean/100)
    nmf.init$H[zero.idx.h] <- runif(length(zero.idx.t,
                          n.threads = n.cores, loss = loss, max.iter = max.iter)
  }
  colnames(nms$H) <- sapply(1:ncol(nmf.res$W),
                                                       function(i) paste("factor", i, sep "))
  return(nmf.res)
}

ica_init <- function (A, k, ica.fast = F)
{
  if (ica.fast) {
    pc.res.h <- irlba(t(A), nv = 50, maxit =                   tol = 1e-04)
    return(list(W = (A - Matrix::rowMeans(A)) %*% ica.res.h$S,
                H = t(ica.res.h$S)))
  }
  
  }
}

