mv_lclme <- function(formulas, data, classes, families, hc, RM_method, predicted, 
                     control = list(), ...){
  
  ########
  # Data #
  ########
  
  components <- lapply(unname(formulas), extractFrames, data = data)
  
  components <- unlist(components, recursive = FALSE)
  n_outcomes <- length(formulas)
  names(components) <- paste0(names(components),
                              rep(seq_len(n_outcomes),
                                  each = length(components) / n_outcomes))
  colmns_HC <- components[grep("colmns_HC", names(components), fixed = TRUE)]
  colmns_nHC <- components[grep("colmns_nHC", names(components), fixed = TRUE)]
  seq_outcomes <- seq_len(n_outcomes)
  
  nams_vars <- c("n", "offset", "Z", "X", "Xc", "means_X", "ncx", "ncz", "y")
    
    # 
    # "N", "id", "X", "Z", "Z_", "Zinv", "Zv", "Ztinv", "X", "Xhc", "ncx", "y", "n", "offset",
    #              "ZrowsStart", "ZrowsEnd","Xc", "Xs", "Zc", "Zs", "XhcC", "XhcS",
    #              "means_X", "SDs_X", "mean_sd_X", "means_Z", "SDs_Z", "mean_sd_Z",
    #              "means_Xhc", "SDs_Xhc", "mean_sd_Xhc", "colmns_nHC")
  vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
  #if (any(ind_td <- sapply(colmns_nHC, length))) {
  #  vars <- c(vars, paste0("X", which(ind_td > 0)))
  #}
  
  Data_data <- c(list(n = components$n1), components[vars])
  
  Data_data$n_RE <- sum(unlist(components[grep("ncz", names(components), fixed = TRUE)]))
  RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                    sq = components[grep("ncz", names(components), fixed = TRUE)],
                    incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                    SIMPLIFY = FALSE)
  names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
  Data1 <- c(Data_data, RE_inds)
  
  ##########
  # Priors #
  ##########
  
  K <- length(families)
  
  datamu <- "string"
  m <- 1
  for (i in 1:classes) {
    datamu[m] <- paste0("mu0", i, " = mu", i, sep = "", collapse = "")
    m <- m + 1
  }
  datamu <- paste0(datamu, sep = "", collapse = ", ")
  
  ncZ <- Data1$n_RE
  # datainvD <- "string"
  # for (i in 1:classes) {
  #    datainvD[i] <- paste0("priorR.D", i, " = diag(priorD, ncZ), priorK.D", i, " = (ncZ)", sep = "", collapse = "")
  # }
  # datainvD <- paste0(datainvD, sep = "", collapse = ", ")
  # 
  datainvD <- paste0("priorR.D = diag(priorD, ncZ), priorK.D = (ncZ + 1)")
  
  databet <- "string"
  m <- 1
  for (k in 1:K) {
    for (i in 1:classes) {
      databet[m] <- paste0("priorMean.betas", k, i, " = betas", k, i, sep = "", collapse = "")
      m <- m + 1
    }
  }
  databet <- paste0(databet, sep = "", collapse = ", ")
  
  datatau <- "string"
  if (any(families == "gaussian")){
    K_gaus <- which(families == "gaussian")
    m <- 1
    for (k in K_gaus) {
      datatau[m] <- paste0("priorA.tau", k, " = priorA.tau", k, ", priorB.tau", k, " = priorB.tau", k, sep = "", collapse = "")
      m <- m + 1
    }
    datatau <- paste0(datatau, sep = "", collapse = ", ")
  }
  
  # CGANCE
  # databetVar <- "string"
  # m <- 1
  # for (k in 1:K) {
  #   for (i in 1:classes) {
  #     databetVar[m] <- paste0("priorTau.betas", k, i, " = diag(1/var.betas", k, i, ")", sep = "", collapse = "")
  #     m <- m + 1
  #   }
  # }
  # databetVar <- paste0(databetVar, sep = "", collapse = ", ")
  
  databetVar <- "string"
  m <- 1
  for (k in 1:K) {
    for (i in 1:classes) {
      databetVar[m] <- paste0("tau_betas", k, i, " = 1/var.betas", k, i, sep = "", collapse = "")
      m <- m + 1
    }
  }
  databetVar <- paste0(databetVar, sep = "", collapse = ", ")
  
  
  
  
  nZ <- Data1[grep("ncz", names(Data1), fixed = TRUE)]
  mu <- rep(0, Data1$n_RE)
  mu <- paste0("c(", paste0(mu, collapse = ","), ")")
  mu_s <- lapply(seq_len(classes), function(x) paste0("mu", x, " <- ", mu))
  mu_s <- paste0(mu_s, sep = "; ", collapse = "")
  eval(parse(text = mu_s))
  
  priorD <- "NA"
  args <- c(1:classes)
  priorD_s <- mapply(x = args, function(x) paste0("priorD", x, " <- ", priorD))
  priorD_s <- paste0(priorD_s, sep = "; ", collapse = "")
  eval(parse(text = priorD_s))
  
  nX <- Data1[grep("ncx", names(Data1), fixed = TRUE)]
  betas <- lapply(nX, function(x) rep(0, x))
  betas <- lapply(betas, function(x) paste0("c(", paste0(x, collapse = ","), ")"))
  betas <- rep(betas, classes)
  args <- expand.grid(x = c(1:K), y = c(1:classes))
  betas_s <- mapply(x = args$x, y = args$y, z = betas, function(x, y, z) paste0("betas", x, y, " <- ", z))
  betas_s <- paste0(betas_s, sep = "; ", collapse = "")
  eval(parse(text = betas_s))
  
  # CHANGE  
  # var.betas <- lapply(nX, function(x) rep(100, x))
  # var.betas <- lapply(var.betas, function(x) paste0("c(", paste0(x, collapse = ","), ")"))
  # var.betas <- rep(var.betas, classes)
  # args <- expand.grid(x = c(1:K), y = c(1:classes))
  # var.betas_s <- mapply(x = args$x, y = args$y, z = var.betas, function(x, y, z) paste0("diag(var.betas", x, y, " <- ", z, ")"))
  # var.betas_s <- paste0(var.betas_s, sep = "; ", collapse = "")
  # eval(parse(text = var.betas_s))
  
  args <- expand.grid(x = c(1:K), y = c(1:classes))
  var.betas_s <- mapply(x = args$x, y = args$y, function(x, y, z) paste0("var.betas", x, y, " <-  100"))
  var.betas_s <- paste0(var.betas_s, sep = "; ", collapse = "")
  eval(parse(text = var.betas_s))
  
  
  
  
  if (any(families == "gaussian")){
    K_gaus <- which(families == "gaussian")
    priorA.tau <- 0.01
    args <- expand.grid(x = c(K_gaus))
    priorA.tau_s <- lapply(args$x, function(x) paste0("priorA.tau", x, " <- ", priorA.tau))
    priorA.tau_s <- paste0(priorA.tau_s, sep = "; ", collapse = "")
    eval(parse(text = priorA.tau_s))
  }
  
  if (any(families == "gaussian")){
    K_gaus <- which(families == "gaussian")
    priorB.tau <- 0.01
    args <- expand.grid(x = c(K_gaus))
    priorB.tau_s <- mapply(x = args$x,function(x) paste0("priorB.tau", x, " <- ", priorB.tau))
    priorB.tau_s <- paste0(priorB.tau_s, sep = "; ", collapse = "")
    eval(parse(text = priorB.tau_s))
  }
  
  
  
  
  data_priors <- paste0("Data2 <- list(", datamu, ",\n", datainvD, ",\n", databetVar, 
                        if (any(families == "gaussian")) { paste0(",\n", datatau) }, ")")
  
  eval(parse(text = data_priors))
  
  Data <- c(Data1, Data2)
  
  ########
  # MCMC #
  ########
  
  con <- list(n.processors = parallel::detectCores() - 1, 
              working.directory = getwd(), clear.model = TRUE,
              seed = 1L, optimize_only = FALSE, verbose = FALSE, 
              n.iter = 28000L, n.burnin = 3000L, n.thin = 50L, 
              n.adapt = 3000L, n.chains = 2L, seed = 1L, n.cores = 1L)
  
  
  
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (!any(namc == "n.thin")) {
    con$n.thin <- if (engine == "JAGS") {
      max(1, floor((con$n.iter - con$n.burnin) * con$n.chains / 1000))
    } else {
      max(1, floor((con$n.iter - con$n.warmup) * con$n.chains / 1000))
    }
  }
  if (length(noNms <- namc[!namc %in% namC]) > 0)
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  ######################
  # Parameters to save #
  ######################
  
  args <- expand.grid(x = c(1:n_outcomes), y = c(1:classes))
  params_betas <- paste0('betas', paste0(args[,1], args[,2]))
  if (any(ind_gs <- sapply(families, function (x) x == "gaussian"))) {
    params_betas <- c(params_betas, paste0("tau", which(ind_gs)), paste0("sigma", which(ind_gs)))
  }
  params_b <- paste0('b', seq_len(classes))
  params_invD <- paste0('inv.D', seq_len(classes))
  params_pred <- paste0('y', seq_len(n_outcomes), '_pred')
  
  params <- c(params_betas, params_b, params_invD, "pr", "v", if (predicted == TRUE) {params_pred}) 
  
  ############
  # Initials #
  ############
  
  # inits <- function () {
  #   ints <- lapply(components[grep("ncx", names(components), fixed = TRUE)],
  #                  rnorm, sd = 0.1)
  #   names(ints) <- paste0('betas', seq_len(n_outcomes))
  #   ints$u <- drop(matrix(rnorm(Data$n * Data$n_RE), Data$n,
  #                         Data$n_RE))
  #   ints$inv_D <- ints$D <- if (Data$n_RE > 1) diag(Data$n_RE) else 1
  #   if (any(ind_gs)) {
  #     nms <- which(ind_gs)
  #     taus <- rep(list(1), length(nms))
  #     names(taus) <- paste0("tau", nms)
  #     ints <- c(ints, taus)
  #   }
  #   ints
  # }
  
  ####################
  # Build JAGS model #
  ####################
  
  JAGSmodel(classes = classes, families, hc, RM_method, predicted)
  
  #############
  # Fit model #
  #############
  
  n_betas <- sum(unlist(lapply(Data[grep("X", names(Data), fixed = TRUE)], function(x) ncol(x))))
  #n_sigma <- length(unlist(Data[grep("tau", names(Data), fixed = TRUE)]))/2
  #n_invD <-   sum(lower.tri(Data$priorR.D1, diag = TRUE))
  a <- n_betas #+ n_invD
  Data$prior.cl <- rep(a/2 - 0.1, classes)
  
  model_name <- "mixedmodel.txt"
  # Data$priorMeanbetas1 <- rep(0, Data$ncx1)
  # Data$priorMeanbetas2 <- rep(0, Data$ncx2)
  # Data$priorMeanbetas3 <- rep(0, Data$ncx3)
  # Data$priorMeanbetas4 <- rep(0, Data$ncx4)
  # Data$priorTaubetas1 <- diag(0.01, Data$ncx1)
  # Data$priorTaubetas2 <- diag(0.01, Data$ncx2)
  # Data$priorTaubetas3 <- diag(0.01, Data$ncx3)
  # Data$priorTaubetas4 <- diag(0.01, Data$ncx4)
  
  # # centering
  # for(o in 1:n_outcomes){
  #   code_cen <- paste0("Data$X", o, "_mean <- as.vector(apply(Data$X", o, ", 2, mean))")
  #   eval(parse(text = code_cen))
  # }
  
  if (classes != 1 ){
    if (RM_method == TRUE) {
      
      # library(rjags)
      # model.fit <- jags.model(file = "mixedmodel.txt", data = Data, n.chains = 1,
      #                         n.adapt = con$n.adapt)
      # update(model.fit, con$n.burnin)
      # res <- coda.samples(model.fit, params, n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
      # codaFit <- as.mcmc.list(res)
      # 
      # bss <- do.call(rbind,codaFit)
      # colnames(bss)
      # n.sims <- nrow(bss)
      # sims.list <- vector("list", length(params))
      # names(sims.list) <- params
      # for (p in seq_along(params)) {
      #   ii <- grep(paste("^", params[p], sep = ""), colnames(bss))
      #   sims.list[[p]] <- bss[, ii]
      # }
      # W <- sims.list$v
      
      fit <- jags(data = Data, parameters.to.save = params, #inits = inits, 
                  model.file = file.path(con$working.directory, model_name),
                  n.chains = con$n.chains, parallel = TRUE, n.cores = con$cores,
                  n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
                  n.thin = con$n.thin, seed = con$seed, verbose = con$verbose)
      
      
      W <- fit$mcmc$v
      N <- dim(W)[2]
      
      ClassSel <- function(n_per, classes) {
        EPTC <- NULL
        for (i in 1:dim(W)[1]){
          EPTC[i] <- sum(sapply(X = 1:classes, FUN = function(x) sum(W[i,] == x) <= n_per*N/100))
        }
        EPTC
      }
      
      
      Mod <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }
      
      cl_0.01 <- classes - Mod(ClassSel(1, classes))
      cl_0.05 <- classes - Mod(ClassSel(5, classes))
      cl_0.10 <- classes - Mod(ClassSel(10, classes))
      cl_0.15 <- classes - Mod(ClassSel(15, classes))
      
      
      out <- c("Non-empty classes for 0.01" = cl_0.01, "Non-empty classes for 0.05" = cl_0.05, 
               "Non-empty classes for 0.10" = cl_0.10, "Non-empty classes for 0.15" = cl_0.15)
      out
    }
    
    
    if (RM_method == FALSE) {
      
      fit <- jags(data = Data, parameters.to.save = params, #inits = inits, 
                  model.file = file.path(con$working.directory, model_name),
                  n.chains = con$n.chains, parallel = TRUE, #con$n.processors > 1, 
                  n.cores = con$cores,
                  n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
                  n.thin = con$n.thin, seed = NULL, verbose = con$verbose)
      
      
      
      #######################
      # fix label switching #
      #######################
      
      
      mcmc.pars <- fit$sims.list[-which(names(fit$sims.list) %in% c("pr", "v", "deviance", params[grep("tau", params)], 
                                                                    params[grep("sigma", params)],  params[grep("inv.D", params)],
                                                                    params[grep("y", params)]))]
      iter_per_chain <- (con$n.iter - con$n.burnin)/con$n.thin
      vec.iter <- seq(1, (iter_per_chain * con$n.chains) + iter_per_chain, by = iter_per_chain)
      for (p in 1:con$n.chains){
        mcmc.chainA <- paste0("mcmc.pars.chain", p, " <- lapply(mcmc.pars[names(mcmc.pars) %in% c(params[grep('betas', params)])], function(x) x[", vec.iter[p], ":", (vec.iter[p+1]-1), ", ])")
        eval(parse(text = mcmc.chainA))
        mcmc.chainB <- paste0("mcmc.pars.chain", p, " <- c(mcmc.pars.chain", p, ", lapply(mcmc.pars[!names(mcmc.pars) %in% c(params[grep('betas', params)])], function(x) x[", vec.iter[p], ":", (vec.iter[p+1]-1), ", , ]))")
        eval(parse(text = mcmc.chainB))
        
        eval(parse(text = paste0("mcmc.pars.chain", p, " <- as.data.frame(mcmc.pars.chain", p, ")")))
      }
      
      
      for (l in 1:con$n.chains){
        eval(parse(text = paste0("mcmc.pars <- mcmc.pars.chain", l)))
        
        class_par <- list()
        for (i in 1:classes) {
          class_par[[i]] <- mcmc.pars[grep(paste0(i, "."), names(mcmc.pars), fixed = TRUE)]
        } 
        
        sort_names <- lapply(class_par, function(x) paste0(1:dim(class_par[[classes]])[2], "_", colnames(x)))
        for (i in 1:classes){
          colnames(class_par[[i]]) <- sort_names[[i]]
        }
        
        DD <- list()
        for (i in 1:dim(class_par[[1]])[2]){
          DD[[i]] <- as.data.frame(lapply(class_par, "[", , i))
          colnames(DD[[i]]) <- unlist(lapply(class_par, function(x) colnames(x)[i]))
        }
        DD <- array(as.numeric(unlist(DD)), dim = c(dim(class_par[[1]])[1], classes, dim(class_par[[1]])[2]))
        
        p <- fit$sims.list$pr[vec.iter[l]:(vec.iter[l+1]-1),,]
        run <- stephens(p)
        
        reordered.mcmc <- permute.mcmc(DD, run$permutations)
        list_reordered.mcmc <- lapply(seq(dim(reordered.mcmc$output)[3]), function(x) reordered.mcmc$output[ , , x])
        
        for (i in seq(1, dim(class_par[[1]])[2])){
          vec_name <- unlist(lapply(class_par, function(x) colnames(x)[i]))
          vec_name <- sub(".*?_","", vec_name)
          #vec_name <- gsub("\\..*","",a)
          colnames(list_reordered.mcmc[[i]]) <- vec_name
        }
        
        mat_reordered.mcmc <- do.call(cbind, list_reordered.mcmc)
        
        mat_betas <- list()
        for (i in 1:length(params_betas)){
          mat_betas[[i]] <- mat_reordered.mcmc[, grep(params_betas[i], colnames(mat_reordered.mcmc), fixed = TRUE)]
        }
        mat_betas  <- mat_betas[unlist(lapply(mat_betas, length) != 0)] 
        namesbetas1 <- params_betas[-(grep(c("tau"), params_betas, fixed = TRUE))]
        names(mat_betas) <- namesbetas1[-(grep(c("sigma"), namesbetas1, fixed = TRUE))]
        
        mat_b <- list()
        for (i in 1:length(params_b)){
          mat_b. <- mat_reordered.mcmc[, grep(params_b[i], colnames(mat_reordered.mcmc), fixed = TRUE)]
          mat_b[[i]] <- array(mat_b., dim = c(dim(mat_b.)[1], 
                                              eval(parse(text = paste0("Data1$n", i))), Data$n_RE))
        }
        names(mat_b) <- params_b
        
        
        # mat_invD <- list()
        # for (i in 1:length(params_invD)){
        #   mat_invD. <- mat_reordered.mcmc[, grep(params_invD[i], colnames(mat_reordered.mcmc), fixed = TRUE)]
        #   mat_invD[[i]] <- array(mat_invD., dim = c(dim(mat_invD.)[1], 
        #                                             Data$n_RE, Data$n_RE))
        # }
        # names(mat_invD) <- params_invD
        
        new_sims.list <- do.call(c, list(mat_betas, mat_b, #mat_invD, 
                                         fit$sims.list[which(names(fit$sims.list) %in% c("pr", "v", "deviance", params[grep("tau", params)], 
                                                                                         params[grep("sigma", params)],
                                                                                         params[grep("y", params)], params[grep("inv.D", params)]))]))
        eval(parse(text = paste0("new_sims.list", l, " <-new_sims.list")))
        
      }
      
      new_sims.list <- new_sims.list1
      # combine chains and fix label switching problem between chains
      ## incomplete (works for 1 chain only at the moment)
      # order(c(median(new_sims.list1$betas11[,1]),median(new_sims.list1$betas12[,1]), median(new_sims.list1$betas13[,1])))
      # order(c(median(new_sims.list2$betas11[,1]),median(new_sims.list2$betas12[,1]), median(new_sims.list2$betas13[,1])))
      # 
      # order(c(median(new_sims.list1$betas11[,2]),median(new_sims.list1$betas12[,2]), median(new_sims.list1$betas13[,2])))
      # order(c(median(new_sims.list2$betas11[,2]),median(new_sims.list2$betas12[,2]), median(new_sims.list2$betas13[,2])))
      
      
      
      # #### test ####
      # head(new_sims.list$betas11)
      # head(fit$sims.list$betas11)
      # 
      # dim(new_sims.list$betas21)
      # dim(fit$sims.list$betas21)
      # 
      # tail(new_sims.list$betas22)
      # tail(fit$sims.list$betas22)
      # 
      # tail(new_sims.list$b1)
      # tail(fit$sims.list$b1)
      # 
      # head(new_sims.list$b2)
      # head(fit$sims.list$b2)
      # 
      # dim(new_sims.list$inv.D1)
      # dim(fit$sims.list$inv.D1)
      # 
      # (new_sims.list$inv.D1)[1:5, 1:2, 3]
      # (fit$sims.list$inv.D1)[1:5, 1:2, 3]
      # 
      # head(new_sims.list$tau1)
      # head(fit$sims.list$tau1)
      # 
      # fit$sims.list$betas11[112,]
      # new_sims.list$betas11[112,]
      # 
      # fit$sims.list$betas12[212,]
      # new_sims.list$betas12[212,]
      # 
      # fit$sims.list$betas21[12,]
      # new_sims.list$betas21[12,]
      # 
      # fit$sims.list$betas22[345,]
      # new_sims.list$betas22[345,]
      # 
      # fit$sims.list$b1[24,40:43,1]
      # new_sims.list$b1[24,40:43,1]
      # 
      # fit$sims.list$b1[24,20:23,2]
      # new_sims.list$b1[24,20:23,2]
      # 
      # fit$sims.list$b2[24,20:23,2]
      # new_sims.list$b2[24,20:23,2]
      # 
      # fit$sims.list$inv.D1[45, , 1]
      # new_sims.list$inv.D1[45, , 1]
      # 
      # fit$sims.list$inv.D2[500, , 3]
      # new_sims.list$inv.D2[500, , 3]
      out <- list(mcmc = new_sims.list, components = components, data = data,
                  families = families, control = con, mcmc.info = fit$mcmc.info[-which(names(fit$mcmc.info) %in% c("end.values"))],
                  DIC = fit$DIC, pD = fit$pD, Rhat = fit$Rhat)
      
      Xnams <- lapply(components[grep("^X[0-9]", names(components))], colnames)
      for (i in seq_along(Xnams)) {
        for (j in 1:classes) {
          colnames(out$mcmc[[paste0("betas", i, j)]]) <- Xnams[[i]]
        }
      }
      #pat <- paste0("^Z", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
      Znams <- lapply(components[grep("^Z[0-9]", names(components))], colnames)
      Znams <- unlist(mapply(paste0, Znams, seq_len(n_outcomes), SIMPLIFY = FALSE))
      #dimnames(out$mcmc$D) <- dimnames(out$mcmc$inv.D) <- list(NULL, Znams, Znams)
      #dimnames(out$mcmc$b) <- list(NULL, NULL, Znams)
      # calculate statistics
      summary_fun <- function (FUN, ...) {
        out <- lapply(out$mcmc, function (x) {
          if (!is.null(dim(x)) && length(dim(x)) > 1) {
            d <- if (is.matrix(x)) 2L else c(2L, 3L)
            apply(x, d, FUN, ...)
          } else {
            FUN(x, ...)
          }
        })
        out[!sapply(out, is.null)]
      }
      out$postMeans <- summary_fun(mean, na.rm = TRUE)
      out$postModes <- summary_fun(modes)
      out$EffectiveSize <- summary_fun(effectiveSize)
      out$StErr <- summary_fun(stdErr)
      out$StDev <- summary_fun(sd, na.rm = TRUE)
      out$CIs <- summary_fun(quantile, probs = c(0.025, 0.975))
      out$Pvalues <- summary_fun(computeP)
      out$call <- formulas
      class(out) <- "mvlclme"
      out
      
    }
    
  } else {
    
    fit <- jags(data = Data, parameters.to.save = params, #inits = inits, 
                model.file = file.path(con$working.directory, model_name),
                n.chains = con$n.chains, parallel = TRUE, n.cores = 7,
                n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
                n.thin = con$n.thin, seed = con$seed, verbose = con$verbose)
    
    out <- list(mcmc = fit$sims.list, components = components, data = data,
                families = families, control = con, mcmc.info = fit$mcmc.info[-which(names(fit$mcmc.info) %in% c("end.values"))],
                DIC = fit$DIC, pD = fit$pD, Rhat = fit$Rhat)#,
    
    Xnams <- lapply(components[grep("^X[0-9]", names(components))], colnames)
    for (i in seq_along(Xnams)) {
      for (j in 1:classes) {
        colnames(out$mcmc[[paste0("betas", i, j)]]) <- Xnams[[i]]
      }
    }
    #pat <- paste0("^Z", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
    Znams <- lapply(components[grep("^Z[0-9]", names(components))], colnames)
    Znams <- unlist(mapply(paste0, Znams, seq_len(n_outcomes), SIMPLIFY = FALSE))
    #dimnames(out$mcmc$D) <- dimnames(out$mcmc$inv.D) <- list(NULL, Znams, Znams)
    #dimnames(out$mcmc$b) <- list(NULL, NULL, Znams)
    # calculate statistics
    summary_fun <- function (FUN, ...) {
      out <- lapply(out$mcmc, function (x) {
        if (!is.null(dim(x)) && length(dim(x)) > 1) {
          d <- if (is.matrix(x)) 2L else c(2L, 3L)
          apply(x, d, FUN, ...)
        } else {
          FUN(x, ...)
        }
      })
      out[!sapply(out, is.null)]
    }
    out$postMeans <- summary_fun(mean, na.rm = TRUE)
    out$postModes <- summary_fun(modes)
    out$EffectiveSize <- summary_fun(effectiveSize)
    out$StErr <- summary_fun(stdErr)
    out$StDev <- summary_fun(sd, na.rm = TRUE)
    out$CIs <- summary_fun(quantile, probs = c(0.025, 0.975))
    out$Pvalues <- summary_fun(computeP)
    out$call <- formulas
    class(out) <- "mvlclme"
    out
  }
  
  
  
  #########################################################################
  ### ??? ### Why fit$mean$inv.D1 has 3 columns and 3 matrices? ### ??? ###
  #########################################################################
  
  
  
  
}

