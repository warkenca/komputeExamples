# Function to simulate from MVN(0, Sigma) for a specified number of variables/samples
simulation <- function(n.samples, pheno.cor, mask.prop, method="svd", info.cutoff=0){

  npheno <- nrow(pheno.cor)
  sim.z <- mvrnorm(n=n.samples, mu=rep(0, npheno), Sigma=pheno.cor)
  sim.z <- t(sim.z)

  ## mask a % of measured z-scores
  nimp <- nrow(sim.z)*ncol(sim.z)*mask.prop
  all.i <- 1:(nrow(sim.z)*ncol(sim.z))

  mask.i <- sort(sample(all.i, nimp))
  org.z = as.matrix(sim.z)[mask.i]
  zvec <- as.vector(as.matrix(sim.z))
  zvec[mask.i] <- NA
  zmat.imp <- matrix(zvec, nrow=npheno)
  rownames(zmat.imp) <- rownames(sim.z)


  # Matrix completion method
  if(method=="svd"){
    svd.res <- svd.impute(zmat.imp, r)
    org.z <- sim.z[mask.i]
    imp.z <- svd.res[mask.i]

  } else if(method=="kompute"){ # KOMPUTE method
    kompute.res <- kompute(zmat.imp, pheno.cor, 0.01)

    imp.z <- as.matrix(kompute.res$zmat)[mask.i]
    imp.info <- as.matrix(kompute.res$infomat)[mask.i]

    imp <- data.frame(org.z=org.z, imp.z=imp.z, info=imp.info)
    imp <- imp[complete.cases(imp),]
    imp <- subset(imp, info>=0 & info <= 1)

    # Info cutoff
    imp.sub <- subset(imp, info>info.cutoff)
    org.z <- imp.sub$org.z
    imp.z <- imp.sub$imp.z
    info <- imp.sub$info
  }

  # Calculate correlation
  cor.val <- round(cor(imp.z, org.z), digits=3)
  type <- ifelse(method=="svd", "Matrix Completion", "KOMPUTE")

  # Create plot
  if(info.cutoff==0){
    plot <- ggplot() +
      geom_point(aes(x=imp.z, y=org.z), alpha=0.1) +
      labs(title=paste0(type, ", ", mask.prop*100, "% Removed, Cor=", cor.val),
           x="Imputed Z-scores", y = "Measured Z-scores") +
      theme_minimal()

  } else{

    plot <- ggplot() +
      geom_point(aes(x=imp.z, y=org.z, col=info), alpha=0.1) +
      labs(title=paste0(type, ", ", mask.prop*100, "% Removed, Cor=", cor.val),
           x="Imputed Z-scores", y = "Measured Z-scores") +
      theme_minimal()
  }

  return(list(plot=plot, cor=cor.val))
}
