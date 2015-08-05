library(nlme) # change the package so that it requires nlme

#' Perform the anti-logit
#' @param M matrix of M values on which to apply anti-logit.
#' @return beta matrix
#' @export
m.ilogit = function(M){
    res = exp(M) / (1 + exp(M))
    rownames(res) = rownames(M)
    colnames(res) = colnames(M)
    res
}

#' Perform the logit transform.
#' @param b matrix of beta values on which to apply logit.
#' @return M matrix
#' @export
m.logit = function(b){
    res = log(b) - log(1 - b)
    rownames(res) = rownames(b)
    colnames(res) = colnames(b)
    res
}

penFitOne = function(y, Zmat){
  adj = y
  is.obs = !is.na(y)
  Z = Zmat[is.obs,,drop=FALSE]
  y = y[is.obs]
  id = rep(1,length(y))
  lmod = try(lme(y~1, random=list(id=pdIdent(~Z-1))), silent=TRUE)
  if(!inherits(lmod,"try-error")){
    adj[is.obs] = resid(lmod) + lmod$coef$fixed[1]
  }# else {
  #    message(lmod)
  #}
  list(modelFit=lmod, adjusted=adj)
}

penFitAll = function(Ymat, Zmat){
  nFeature = dim(Ymat)[1]

  mu = sigma = tau = rep(NA, nFeature)
  beta = matrix(NA, nFeature, dim(Zmat)[2])
  adjusted = matrix(NA, nFeature, dim(Ymat)[2])
  nbad = 0
  for(i in 1:nFeature){
    pf = penFitOne(Ymat[i,], Zmat)

    if(inherits(pf$modelFit,"try-error")){
      nbad = nbad + 1
      mu[i] = mean(Ymat[i,],na.rm=TRUE)
      sigma[i] = var(Ymat[i,],na.rm=TRUE)
      tau[i] = 0
      beta[i,] = 0
    }
    else{
      mu[i] = pf$modelFit$coef$fixed[1]
      sigma[i] = pf$modelFit$sigma^2
      tau[i] = getVarCov(pf$modelFit)[1,1]
      beta[i,] = pf$modelFit$coef$random$id
    }
    adjusted[i,] = pf$adjusted
  }
  message(paste("total with error in model fit:", nbad))
  list(mu=mu, beta=beta, tau=tau, sigma=sigma, adjusted=adjusted)
}

penFit <- function(y, Zmat) {
  adj <- y
  is.obs <- !is.na(y)
  Z <- Zmat[is.obs,,drop=FALSE]
  y <- y[is.obs]
  id <- rep_len(1,length(y))
  posdMat <- pdIdent(~Z-1)
  lmod <- try(lme(y~1, random=list(id=posdMat)), silent=TRUE)
  
  if(inherits(lmod,"try-error")){
    mu <- mean(y, na.rm=TRUE)
    sigma <- var(y, na.rm=TRUE)
    tau = 0
    beta = 0
  }
  else{
    mu <- lmod$coef$fixed[1]
    sigma <- lmod$sigma^2
    tau <- getVarCov(lmod)[1,1]
    beta <- lmod$coef$random$id
    adj[is.obs] <- resid(lmod) + mu
  }

  list(mu=mu, beta=beta, tau=tau, sigma=sigma, 
       adjusted=adj, bad=inherits(lmod, "try-error"))
}

penFitOne2 = function(y, Zmat, fit){
  adj = y
  Z = Zmat
  id = rep_len(1, length(adj))
  lmod = try(update(fit, adj~1), silent=TRUE)
  if(!inherits(lmod,"try-error")){
    adj = resid(lmod) + lmod$coef$fixed[1]
  }# else {
  #    message(lmod)
  #}
  list(modelFit=lmod, adjusted=adj)
}

penFitAll2 = function(Ymat, Zmat){
  nFeature = dim(Ymat)[1]
  
  mu = sigma = tau = rep(NA, nFeature)
  beta = matrix(NA, nFeature, dim(Zmat)[2])
  adjusted = matrix(NA, nFeature, dim(Ymat)[2])
  nbad = 0
  
  pf = penFitOne(Ymat[1,], Zmat)
  fit = pf$modelFit
  
  mu[1] = fit$coef$fixed[1]
  sigma[1] = fit$sigma^2
  tau[1] = getVarCov(fit)[1,1]
  beta[1,] = fit$coef$random$id
  adjusted[1,] <- pf$adjusted
  
  for(i in 2:nFeature){
    pf = penFitOne2(Ymat[i,], Zmat, fit)
    
    if(inherits(pf$modelFit,"try-error")){
      nbad = nbad + 1
      mu[i] = mean(Ymat[i,],na.rm=TRUE)
      sigma[i] = var(Ymat[i,],na.rm=TRUE)
      tau[i] = 0
      beta[i,] = 0
    }
    else{
      mu[i] = pf$modelFit$coef$fixed[1]
      sigma[i] = pf$modelFit$sigma^2
      tau[i] = getVarCov(pf$modelFit)[1,1]
      beta[i,] = pf$modelFit$coef$random$id
    }
    adjusted[i,] = pf$adjusted
  }
  message(paste("total with error in model fit:", nbad))
  list(mu=mu, beta=beta, tau=tau, sigma=sigma, adjusted=adjusted)
}


#' Adjust Beta for cell mixture. Returns beta of the average cell-type.
#'
#' Can also return the estimates of each cell when est.only=TRUE
#'
#' @param B matrix of beta values (n_probes * n_samples) rownames must
#'        be the illumina cg id or the chrom:position
#' @param top_n number of probes from cell.coefs to use
#' @param mc.cores number of cores to use for parallelization
#' @param cell.coefs path to tab-delimited file containind DMR-probes for the
#'        cell types of interest. system.file("extdata", "houseman-dmrs.txt",
#'        package="celltypes450")
#' @param est.only if TRUE, only estimate the coefficients for the cell-types
#'        without actually adjusting.
#' @return matrix of adjusted beta values (or data.frame of cell coefs if)
#'         \code{est.only} is FALSE
#' @export
adjust.beta = function(B, top_n=500, mc.cores=2, 
                       cell.coefs=NULL,
                       est.only=FALSE, version=NULL){
    stopifnot(all((B > 0) & (B < 1), na.rm=TRUE))

    if(is.null(cell.coefs)){
      if(nrow(B) <= 200000){
        if(length(grep("^cg", rownames(B), value=TRUE, perl=TRUE)) >=
           length(grep("^chr", rownames(B), value=TRUE, perl=TRUE))){
          cell.coefs = system.file("extdata", "houseman-dmrs.txt", package="celltypes450")
        } else {
          cell.coefs = system.file("extdata", "houseman-dmrs-locs.txt", package="celltypes450")
        }
        warning("We suggest you use at least 2x10^5 probes")
      } else {
        # rownames are cg id.
        if(length(grep("^cg", rownames(B), value=TRUE, perl=TRUE)) > 200000){
            cell.coefs = system.file("extdata", "houseman-dmrs.txt", package="celltypes450")
        } else { # row names are e.g. chr17:12345
            cell.coefs = system.file("extdata", "houseman-dmrs-locs.txt", package="celltypes450")
        }
      }
    }
    # after adjusting, set values < 0 or > 1 to the smallest observed value
    epsilon.min = min(B, 1 - B)
    epsilon.max = 1 - epsilon.min
    
    dmr.coefs = read.delim(cell.coefs, row.names=1)
    # take shared probes. may be differences if beta is from 450k
    dmr.coefs = dmr.coefs[rownames(dmr.coefs) %in% rownames(B),]
    dmr.coefs = as.matrix(dmr.coefs[1:min(top_n, nrow(dmr.coefs)),])
    stopifnot(nrow(dmr.coefs) > 0)

    # reporter matrix
    Lwbc = diag(ncol(dmr.coefs))[-(1:2),]
    Lwbc[1:2,2] = 1 ; Lwbc[,1] = 1
    rownames(Lwbc) = colnames(dmr.coefs)[-(1:2)]
    colnames(Lwbc) = colnames(dmr.coefs)
    dmrs = rownames(dmr.coefs)

    message(paste0("using ", length(dmrs), " dmrs"))

    omega.mix = projectWBC(B[dmrs,, drop=FALSE],
        dmr.coefs[dmrs,],
        contrastWBC=Lwbc,
        nonnegative=TRUE,
        lessThanOne=FALSE)

    omega.mix = (1 / rowSums(omega.mix)) * omega.mix
    if(est.only){ return(omega.mix) }

    message("adjusting beta (this will take a while)...")
    if(!is.null(version) && version == "devel") {
      tmpAdj <- parallel::mclapply(1:nrow(B), function(i){ 
        penFit(B[i,, drop=FALSE], omega.mix) 
        }, mc.cores=mc.cores)
      nbad <- sum(vapply(tmpAdj, function(myList) myList$bad, logical(1)))
      message(paste("total with error in model fit:", nbad))
      
      adjBeta <- vapply(tmpAdj, function(myList) myList$adjusted, rep(0, ncol(B)))
      adjBeta <- t(adjBeta)
    } else {
      if(!is.null(version) && version == "devel2"){
        tmpList = lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(B), by=mc.cores) })
        tmpAdj = parallel::mclapply(tmpList, function(ix){ penFitAll2(B[ix,], omega.mix) }, mc.cores=mc.cores)
        
        adjBeta = matrix(NA, nrow(B), ncol(B))
        for (i in 1:length(tmpList)){
          adjBeta[tmpList[[i]],] = tmpAdj[[i]]$adjusted
        }
      } else {
        tmpList = lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(B), by=mc.cores) })
        tmpAdj = parallel::mclapply(tmpList, function(ix){ penFitAll(B[ix,], omega.mix) }, mc.cores=mc.cores)
        
        adjBeta = matrix(NA, nrow(B), ncol(B))
        for (i in 1:length(tmpList)){
          adjBeta[tmpList[[i]],] = tmpAdj[[i]]$adjusted
        }
      }
    }
    
    nCpGs = nrow(B) * ncol(B)
    message(paste("% values <= 0:", sum(adjBeta <= 0, na.rm=TRUE) / nCpGs * 100))
    message(paste("% values >= 1:", sum(adjBeta >= 1, na.rm=TRUE) / nCpGs * 100))
    message("these will be changed to 0+epsilon, 1-epsilon respectively")

    adjBeta[(adjBeta <= 0)] = epsilon.min
    adjBeta[(adjBeta >= 1)] = epsilon.max
    rownames(adjBeta) = rownames(B)
    colnames(adjBeta) = colnames(B)

    return(adjBeta)
}
