
gwasreml <- function (baseModel, ...)
UseMethod("gwasreml")

gwasreml.default <- function(baseModel, ...)
    stop("Currently the only supported method is \"asreml\"")

gwasreml.asreml <- function (baseModel, genObj, merge.by = NULL, Trait = NULL, covariate = NULL, type = "main+int", fix.lines = TRUE, chr = names(nmar(genObj)), gen.type = "interval", n.fa = 0, breakout = -1, thresh = "b-corr", TypeI = 0.05, qtl.window = 40, effects.window = 20, trace = TRUE, ...)
{
    qtlcall <- match.call()
    baseObj <- baseModel$QTL
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
        if(!baseModel$converge)
            stop("Base model not converged: Check base model before proceeding with QTL analysis.")
    }
    asremlEnv <- lapply(baseModel$formulae, function(el) attr(el, ".Environment"))
    phenoData <- eval(baseModel$call$data)
    if (missing(genObj))
        stop("genObj is a required argument.")
    if (!inherits(genObj, "cross"))
        stop("genObj is not of class \"cross\"")
    if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
    if (is.null(glines <- genObj$pheno[, merge.by]))
        stop("Genotypic data does not contain column \"", merge.by,
             "\".")
    if (is.null(plines <- phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
             "\".")
    if (all(is.na(match(glines, plines))))
        stop("Names in genotypic \"", merge.by, "\" column do not match any names in phenotypic \"",
             merge.by, "\" column.")
    if(!is.numeric(breakout) | breakout < -1 | breakout == 0)
        stop("breakout argument must be -1 or a positive integer.")
    if(!is.null(covariate)){
        if(!(covariate %in% names(phenoData)))
            stop("Named covariate argument not in data.")
        if(type %in% c("main","int"))
            stop("Scan type can only be \"main+int\", when covariate is non-NULL.")
    }
    if(is.null(Trait)){
        if(is.null(covariate)){
            if(type != "main")
                stop("Scan type can only be \"main\" when Trait = NULL and covariate = NULL.")
        }
        if(n.fa > 0){
            warning("Number of factors cannot be greater than zero when Trait = NULL. Setting n.fa = 0.")
            n.fa <- 0
        }
       lhs <- NULL; lab <- NULL; section <- FALSE
    } else {
        if(!(Trait %in% names(phenoData)))
            stop("Named Trait argument not in data.")
        if((n.trait <- length(levels(phenoData[, Trait]))) > 2) {
            n.par.fa <- (n.fa+1)*n.trait - n.fa*(n.fa-1)/2
            n.par.us <- n.trait*(n.trait+1)/2
            if(n.par.fa > n.par.us)
                stop('n.fa set too high: reset and try again\n')
        }
        lhs <- paste("diag(", Trait, ")", sep = "")
        lab <- paste(" x ", Trait, sep = "")
        section <- TRUE
    }
    if(!is.null(chr) && any(!(chr %in% names(nmar(genObj)))))
        stop("Some chromosome names do not exist inside genObj.")
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }
    if(gen.type %in% "interval"){
        gdat <- lapply(genObj$geno, function(el) el$interval.data)
        if(any(sapply(gdat, is.null)))
            stop("No interval data found in at least one chromosome.")
    }
    else gdat <- lapply(genObj$geno, function(el) el$imputed.data)
    genoData <- do.call("cbind", gdat)
    nint <- lapply(gdat, function(el) 1:ncol(el))
    lint <- unlist(lapply(nint, length))
    mnams <- paste("Chr", rep(names(genObj$geno), times = lint), unlist(nint), sep = ".")
    dimnames(genoData) <- list(as.character(glines), mnams)
    genoData <- genoData[rownames(genoData) %in% as.character(plines),]
    rterms <- unlist(strsplit(deparse(baseModel$call$random[[2]], width.cutoff = 500), " \\+ "))
    forw <- paste(Trait,".*", merge.by, sep = "")
    reve <- paste(merge.by,".*", Trait, sep = "")
    pterm <- rterms[grep(paste(forw, reve, sep = "|"), rterms)]
    rterms <- rterms[!(rterms %in% pterm)]
    whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)
    genetic.term <- merge.by
    vm <- FALSE
    if(!all(whg) & fix.lines){
#        !all(whg <- levels(phenoData[[merge.by]]) %in% rownames(genoData)))
        phenoData$Gomit <- phenoData$Gsave <- plines
        levels(phenoData$Gsave)[!whg] <- NA
        levels(phenoData$Gomit)[whg] <- "GEN"
        fix.form <- as.formula(paste(". ~ ", pasteT(Trait, "Gomit")," + .", sep = ""))
        pterm <- gsub(merge.by, "Gsave", pterm)
        ran.base <- formula(paste("~ ", paste(c(pterm, rterms), collapse = " + "), sep = ""))
        baseModel$call$data <- quote(phenoData)
        cat("\nFixing lines and updating initial base model:\n")
        cat("============================================\n")
        baseModel <- update(baseModel, fixed. = fix.form, random. = ran.base, ...)
        merge.by <- "Gsave"
    }
    if(n.fa > -1){
        message("\nQTL",lab," Diagonal Random effects model.")
        cat("========================================\n")
        qtlModel <- baseModel
        if(ncol(genoData) > nrow(genoData)){
            cov.env <- constructCM(genoData)
            covObj <- cov.env$relm
            qterm <- pasteT(lhs, paste("vm","(",merge.by,", covObj)", sep = ""))
            ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
            attr(genObj, "env") <- cov.env
            vm <- TRUE
        } else {
            covObj <- cbind.data.frame(rownames(genoData), genoData)
            names(covObj)[1] <- merge.by
            qtlModel$call$mbf$ints$key <- rep(merge.by, 2)
            qtlModel$call$mbf$ints$cov <- "covObj"
            qterm <- pasteT(lhs, paste("mbf('ints')", sep = ""))
            ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
        }
        assign("covObj", covObj, envir = parent.frame())
        qtlModel$call$data <- quote(phenoData)
        qtlModel <- update(qtlModel, random. = ran.form, ...)
        if(n.fa > 0){
            qsp <- strsplit(qterm, ":")
#            lhs <- sapply(qsp, "[", 1)
            if(n.trait == 2){
                lhs <- gsub("diag", "us", lhs)
                qterm <- paste(lhs, sapply(qsp, "[", 2), sep = ":")
                ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
                message("\nQTL x ",Trait," Bivariate Random Effects model.")
                cat("===============================================\n")
                qtlModel <- update(qtlModel, random. = ran.form, ...)
            }
            else {
                for(i in 1:n.fa){
                    if(length(grep("diag", lhs))){
                        lhs <- gsub("diag", "fa", lhs)
                        lhs <- gsub(")", ",1)", lhs)
                    } else {
                        old <- paste(i - 1, ")", sep = "")
                        new <- paste(i, ")", sep = "")
                        lhs <- gsub(old, new, lhs)
                    }
                    qterm <- paste(lhs, sapply(qsp, "[", 2), sep = ":")
                    ran.form <- as.formula(paste(c("~", qterm, pterm, rterms), collapse = " + "))
                    message("\nQTL x ",Trait," Factor Analytic (",i,") Random Effects model.")
                    cat("===================================================\n")
                    qtlModel <- update(qtlModel, random. = ran.form, ...)
                }
            }
        }
    }
    outObj <- list()
    chr <- chr[!(chr %in% names(baseObj$geno))]
    if(!is.null(covariate))
        Trait <- covariate
    if(length(chr)){
        for(c in 1:length(chr)){
            tempModel <- qtlModel
            chrn <- chr[c]
            cnams <- paste("Chr", chrn, nint[[chrn]], sep = ".")
            whn <- mnams %in% cnams
            if(n.fa > -1){
                genoRel <- genoData[,!whn]
                if(ncol(genoRel) > nrow(genoRel)){
                    cov.env <- constructCM(genoRel)
                    covObj <- cov.env$relm
                    attr(genObj, "env") <- cov.env
                }
                else {
                    covObj <- cbind.data.frame(rownames(genoRel), genoRel)
                    names(covObj)[1] <- merge.by
                    if(is.null(tempModel$call$mbf$ints) & vm){
                        attr(genObj, "env") <- NULL
                        rterms <- unlist(strsplit(deparse(tempModel$call$random[[2]], width.cutoff = 500), " \\+ "))
                        rterms <- rterms[!(rterms %in% qterm)]
                        qterm <- pasteT(lhs, paste("mbf","(",merge.by,", covObj)", sep = ""))
                        tempModel$call$mbf$ints$key <- rep(merge.by, 2)
                        tempModel$call$mbf$ints$cov <- "covObj"
                        ran.form <- as.formula(paste(c("~ ", qterm, rterms), collapse = " + "))
                        tempModel$call$random <- ran.form
                    }
                }
            }
            assign("covObj", covObj, envir = parent.frame())
            genoSet <- genoData[,whn]
            xmark <- colnames(genoSet) <- gsub("Chr","X", colnames(genoSet))
            genoSet <- cbind.data.frame(rownames(genoSet), genoSet)
            names(genoSet)[1] <- merge.by
            tempDat <- phenoData
            tempDat$ord <- 1:nrow(tempDat)
            tempDat <- merge(tempDat, genoSet, all.x = TRUE, by = merge.by)
            tempDat <- tempDat[order(tempDat$ord),]
            tempModel$call$data <- quote(tempDat)
            waldObj <- list()
            asreml.options(Cfixed = TRUE)
            for(i in 1:length(xmark)){
                tempModel$call$fixed <- qtlModel$call$fixed
                if(type %in% "int")
                    fix.mark <- paste(Trait, xmark[i], sep = ":")
                else if(type %in% "main+int"){
                    fix.mark <- paste(Trait, xmark[i], sep = ":")
                    fix.mark <- paste(xmark[i], fix.mark, sep = "+")
                } else fix.mark <- xmark[i]
                fix.form <- as.formula(paste(". ~ . + ", fix.mark, sep = ""))
                tempModel <- update(tempModel, fixed. = fix.form, ...)
                list.coefs <- tempModel$coefficients$fixed
                zind <- cind <- grep(xmark[i], rownames(list.coefs))
                cf <- list.coefs[cind, , drop = FALSE]
                if(length(int <- grep(":", rownames(cf)))){
                    cfi <- cf[int, 1]
                    iind <- cind[int]
                    izero <- list(coef = iind[cfi != 0], type = "zero")
                    wt <- waldTest(tempModel, cc = list(izero))
                    waldObj$Interaction[[i]] <- c(length(cfi), wt$Zero[1,1])
                    cind <- cind[-int]
                }
                if(length(cind)){
                    mzero <- list(coef = cind, type = "zero")
                    wt <- waldTest(tempModel, cc = list(mzero))
                    waldObj$Main[[i]] <- c(1, wt$Zero[1,1])
                }
                if(type %in% c("main","int")){
                    sub.list <- rev(list.coefs[zind, 1])
                    names(sub.list) <- rev(rownames(list.coefs)[zind])
                    outObj$geno[[chrn]]$coef[[i]] <- sub.list
                    outObj$geno[[chrn]]$vcoef[[i]] <- rev(tempModel$vcoeff$fixed[zind])
                }
                if(breakout == i)
                    break
            }
            waldObj <- lapply(waldObj, function(el){
                el <- do.call("rbind.data.frame", el)
                names(el) <- c("df","wd")
                el })
            atype <- names(waldObj)
            nr <- sapply(waldObj, nrow)
            waldObj <- do.call("rbind.data.frame", waldObj)
            waldObj$mark <- rep(xmark[1:nr], length(atype))
            dist <- genObj$geno[[chrn]]$map
            if ((gen.type == "interval") & (length(dist) > 1))
                dist <- dist[2:length(dist)] - diff(dist)/2
            waldObj$dist <- rep(dist[1:nr], length(atype))
            waldObj$chr <- chrn
            waldObj$type <- rep(atype, times = nr)
            outObj$geno[[chrn]]$wald <- waldObj
            if(type %in% c("main","int")){
                if(!is.null(Trait)){
                    outObj$geno[[chrn]]$coef <- t(do.call("cbind", outObj$geno[[chrn]]$coef))
                    outObj$geno[[chrn]]$vcoef <- t(do.call("cbind", outObj$geno[[chrn]]$vcoef))
                } else {
                    outObj$geno[[chrn]]$coef <- matrix(outObj$geno[[chrn]]$coef, ncol = 1)
                    colnames(outObj$geno[[chrn]]$coef) <- "Expt"
                    outObj$geno[[chrn]]$vcoef <- matrix(outObj$geno[[chrn]]$vcoef, ncol = 1)
                }
            }
        }
    }
    if(!is.null(baseObj$geno))
        outObj$geno <- c(outObj$geno, baseObj$geno)
    if(thresh %in% "b-corr"){
        subObj <- subset(genObj, chr = names(outObj$geno))
        if(gen.type %in% "interval")
            cdat <- do.call("cbind", lapply(subObj$geno, function(el) el$interval.data))
        else cdat <- do.call("cbind", lapply(subObj$geno, function(el) el$imputed.data))
        eig <- eigen(cor(cdat))$values
        core <- length(eig[eig > 1]) + sum(abs(eig)[abs(eig) > 1e-10] - floor(abs(eig)[abs(eig) > 1e-10]))
        adj <- TypeI/core
    }
    outObj$geno <- lapply(outObj$geno, function(el, qtl.window, effects.window, adj){
        if(type %in% c("main","int"))
            el$mark.thresh <- qchisq(1 - adj, 1)
        el$wald$thresh <- qchisq(1 - adj, el$wald$df)
        el$peaks <- findQTLPeaks(el$wald, qtl.window, effects.window)
        el
    }, qtl.window, effects.window, adj)
    peaks <- lapply(outObj$geno, function(el) el$peaks)
    peaks <- peaks[!sapply(peaks, is.null)]
    if(length(peaks)){
        tempModel <- qtlModel
        peaks <- do.call("rbind", peaks)
        state <- rep(1, ncol(genoData))
        names(state) <- colnames(genoData)
        schr <- sapply(strsplit(names(state), "\\."), "[", 2)
        marks <- strsplit(as.character(peaks$mark), "\\.")
        mchr <- sapply(marks, "[", 2)
        mint <- as.numeric(sapply(marks, "[", 3))
        for(i in 1:nrow(peaks)){
            wnams <- names(state)[schr %in% mchr[i]]
            inums <- as.numeric(sapply(strsplit(wnams, "\\."), "[", 3))
            dists <- unique(outObj$geno[[mchr[i]]]$wald$dist)
            dists <- dists[inums]
            exc <- wnams[abs(dists - dists[mint[i]]) <= qtl.window]
            state[exc] <- 0
        }
        genoFin <- genoData[,as.logical(state)]
        if(ncol(genoFin) > nrow(genoFin)){
                cov.env <- constructCM(genoFin)
                covObj <- cov.env$relm
                attr(genObj, "env") <- cov.env
            }
            else {
                covObj <- cbind.data.frame(rownames(genoFin), genoFin)
                names(covObj)[1] <- merge.by
                if(is.null(tempModel$call$mbf$ints) & vm){
                    attr(genObj, "env") <- NULL
                    rterms <- unlist(strsplit(deparse(tempModel$call$random[[2]], width.cutoff = 500), " \\+ "))
                    rterms <- rterms[!(rterms %in% qterm)]
                    qterm <- pasteT(lhs, paste("mbf","(",merge.by,", covObj)", sep = ""))
                    tempModel$call$mbf$ints$key <- rep(merge.by, 2)
                    tempModel$call$mbf$ints$cov <- "covObj"
                    ran.form <- as.formula(paste(c("~ ", qterm, rterms), collapse = " + "))
                    tempModel$call$random <- ran.form
                }
            }
        assign("covObj", covObj, envir = parent.frame())
        genoQTL <- genoData[,gsub("X","Chr", as.character(peaks$mark)), drop = FALSE]
        colnames(genoQTL) <- gsub("Chr","X", colnames(genoQTL))
        genoQTL <- cbind.data.frame(rownames(genoQTL), genoQTL)
        names(genoQTL)[1] <- merge.by
        phenoData$ord <- 1:nrow(phenoData)
        phenoData <- merge(phenoData, genoQTL, by = merge.by, all.x = TRUE)
        phenoData <- phenoData[order(phenoData$ord),]
        message("\nFinal putative QTL",lab," model.")
        cat("=====================================\n")
        if(type %in% "main")
            final.terms <- as.character(peaks$mark)
        else if(type %in% "int"){
            terms.int <- paste(Trait, peaks$mark, sep = ":")
            final.int <- paste(terms.int, collapse = " + ")
            final.main <- paste(peaks$mark, collapse = " + ")
            fix.form <- as.formula(paste(". ~ . ", final.main, final.int, sep = " + "))
            qtlModel <- update(tempModel, fixed. = fix.form, data = phenoData, ...)
            list.coefs <- qtlModel$coefficients$fixed
            peak.test <- list()
            for(i in 1:nrow(peaks)){
                forw <- paste(Trait,".*", peaks$mark[i], sep = "")
                reve <- paste(peaks$mark[i],".*", Trait, sep = "")
                zind <- grep(paste(forw, reve, sep = "|"), rownames(list.coefs))
                ci <- list.coefs[zind, 1]
                peak.test[[i]] <- list(coef = zind[ci != 0], type = "zero")
            }
            wt <- waldTest(qtlModel, cc = peak.test)
            final.terms <- ifelse(wt$Zero[,2] > 0.05, as.character(peaks$mark), terms.int)
        } else {
            if(!is.null(covariate)){
                main.terms <- sapply(strsplit(final.terms, ":"), function(el) el[grep("X\\.", el)])
                final.terms <- unique(c(main.terms, final.terms))
            } else final.terms <- ifelse(peaks$type %in% "Main", peaks$mark, paste(Trait, peaks$mark, sep = ":"))
        }
        fix.form <- as.formula(paste(". ~ . + ", paste(final.terms, collapse = " + "), sep = ""))
        qtlModel <- update(tempModel, fixed. = fix.form, data = phenoData, ...)
    } else peaks <- NULL
    qtl.list <- list()
    qtl.list$call <- qtlcall
    qtl.list$Trait <- Trait
    qtl.list$section <- section
    qtl.list$peaks <- peaks
    qtl.list$type <- gen.type
    qtl.list$qtl.window <- qtl.window
    qtl.list$effects.window <- effects.window
    qtl.list$geno <- outObj$geno
    qtlModel$QTL <- qtl.list
    class(qtlModel) <- c("gwasreml", "asreml")
    qtlModel
}

pasteT <- function(lhs, rhs){
    if(is.null(lhs)) rhs
    else paste(lhs, rhs, sep = ":")
}

update.gwasreml <- function(object, ...){
    if (is.null(newcall <- object$QTL$call))
        stop("need an object with call component or attribute")
    tempcall <- list(...)
    if(!is.null(tempcall$chr)){
        if(!is.null(chrq <- names(object$QTL$geno)))
            tempcall$chr <- c(tempcall$chr, chrq)
    }
    newcall$baseModel <- eval(newcall$baseModel, parent.frame())
    newcall$baseModel$QTL <- object$QTL
    if (length(tempcall)) {
        what <- !is.na(match(names(tempcall), names(newcall)))
        for (z in names(tempcall)[what]) newcall[[z]] <- tempcall[[z]]
        if (any(!what)) {
            newcall <- c(as.list(newcall), tempcall[!what])
            newcall <- as.call(newcall)
        }
    }
    eval(newcall, sys.parent())
}

findQTLPeaks <- function(wald, qtl.exc, eff.exc){
    walds <- split(wald, wald$type)
    peaks <- lapply(walds, function(el, qtl.exc){
        temp <- el
        thr <- unique(el$thresh)
        pos <- pracma::findpeaks(el$wd, minpeakheight = thr, npeaks = 20, sortstr =  TRUE)[,2]
        if(!is.null(pos)){
            el <- el[pos,]
            el$position <- pos
            if(nrow(el) > 1){
                i <- 1
                repeat{
                    if(nrow(el) <= i)
                        break
                    ind <- (i + 1):length(el$dist)
                    wh <- abs(el$dist[i] - el$dist[ind]) < qtl.exc
                    if(any(wh)){
                        omit <- ind[wh]
                        el <- el[-omit,,drop = FALSE]
                    }
                    i <- i + 1
                }
            }
            if(nrow(el) > 1){
                cb <- t(combn(nrow(el), 2))
                logi <- logical()
                for(i in 1:nrow(cb)){
                    marks <- as.character(el$mark[cb[i,]])
                    wind <- (1:nrow(temp))[temp$mark %in% marks]
                    wwd <- temp$wd[wind[1]:wind[2]]
                    logi[i] <- all(wwd > thr)
                }
                if(any(logi)){
                    uni <- unique(cb[logi,2])
                    el <- el[-uni,]
                }
            }
            el
        } else NULL
    }, qtl.exc)
    peaks <- peaks[!sapply(peaks, is.null)]
    if(!length(peaks))
        return(NULL)
    if(length(peaks) > 1){
        pm <- peaks[["Main"]]$dist
        pi <- peaks[["Interaction"]]$dist
        op <- outer(pm, pi, function(X, Y) abs(X - Y))
        op <- apply(op, 1, function(el) all(el > eff.exc))
        peaks[["Main"]] <- peaks[["Main"]][op,]
    }
    do.call("rbind", peaks)
}

plotProfile <- function(object, genObj, chr = names(object$QTL$geno), by.trait = FALSE, LOD = TRUE, thresh = TRUE, use.dist = TRUE, chr.lines = TRUE, annotate.peaks = TRUE){
    if (missing(genObj))
        stop("genObj is a required argument")
    if (!inherits(object, "gwasreml"))
        stop("object is not of class \"gwasreml\"")
    if (!inherits(genObj, "cross"))
        stop("genObj does not inherit class \"cross\"")
    if(any(!(chr %in% names(object$QTL$geno))))
        stop("Some chromosome names do not exist inside object.")
    genObj <- subset(genObj, chr = chr)
    object$QTL$geno <- object$QTL$geno[chr]
    gtype <- object$QTL$type
    nchr <- length(nmar(genObj))
    if(!by.trait){
        waldl <- lapply(object$QTL$geno, function(el) el$wald)
        wald <- do.call("rbind.data.frame", waldl)
        wald <- wald[order(wald$type, wald$chr),]
        grps <- length(unique(wald$type))
    }
    else {
        if(is.null(object$QTL$geno[[1]]$coef))
            stop("Complete marker effects can only be estimated if type = \"main\" or \"int\".")
        zrat <- lapply(object$QTL$geno, function(el){
            sigma2 <- ifelse(object$QTL$section, 1, object$sigma2)
            zr <- (el$coef^2)/el$vcoef*sigma2
            colnames(zr) <- sapply(strsplit(colnames(el$coef), ":"), "[", 1)
            zr
        })
        chrm <- rep(names(object$QTL$geno), times = sapply(zrat, nrow))
        distm <- unlist(lapply(object$QTL$geno, function(el) el$wald$dist))
        markm <- unlist(lapply(object$QTL$geno, function(el) el$wald$mark))
        zratd <- do.call("rbind.data.frame", zrat)
        wald <- cbind.data.frame(mark = markm, dist = distm, chr = chrm)
        wald <- cbind.data.frame(wald, thresh = object$QTL$geno[[1]]$mark.thresh, zratd)
        wald <- reshape::melt.data.frame(wald, id.vars = 1:4, variable_name = "type")
        names(wald)[6] <- "wd"
        object$QTL$peaks <- wald[wald$mark %in% object$QTL$peaks$mark,]
        grps <- ncol(zratd)
    }
    if(use.dist){
        cint <- c()
        cint[1] <- 0
        if(nchr > 1){
            chrl <- chrlen(genObj)
            for(i in 2:nchr){
                cint[i] <- cint[i - 1] + chrl[i - 1] + 3
                whc <- wald$chr == names(chrl)[i]
                wald$dist[whc] <- cint[i] + wald$dist[whc]
            }
        }
    }
    else {
        if(gtype == "interval")
            dist <- 1:sum(sapply(genObj$geno, function(el) ncol(el$interval.data)))
        else dist <- 1:totmar(genObj)
        wald$dist <- rep(dist, grps)
    }
    y.lab <- "Wald"
    if(LOD){
        wald$wd <- 0.5 * log(exp(wald$wd), base = 10)
        y.lab <- "LOD"
    }
    udist <- unique(wald$dist)
    gp <- ggplot(wald, aes_string(x = "dist", y = "wd", colour = "chr")) +
            facet_wrap( ~ type, ncol = 1) + geom_line(size = 1.3) +
            scale_x_continuous(breaks = udist, labels = rep("", length(udist))) +
            ylab(y.lab) + xlab("") + theme_scatter()
    wdc <- wald[,c("chr","dist")]
    wdc <- wdc[!duplicated(wdc),]
    if(chr.lines){
        if(length(chr) > 1){
            maxd <- tapply(wdc$dist, list(wdc$chr), function(el) el[length(el)])
            maxd <- maxd[1:(length(maxd) - 1)] + 0.5
            gp <- gp + geom_vline(xintercept = maxd, colour = "grey80", size = 1)
        }
    }
    if(thresh){
        if(LOD)
            wald$thresh <- 0.5 * log(exp(wald$thresh), base = 10)
        threshd <- wald[,c("type", "thresh")]
        threshd <- threshd[!duplicated(threshd),]
        gp <- gp + geom_hline(data = threshd, aes(yintercept = thresh), colour = "grey50", size = 1, linetype = 2)
    }
    if(annotate.peaks & !is.null(object$QTL$peaks)){
        peaks <- object$QTL$peaks
        peaks <- peaks[peaks$chr %in% chr,, drop = FALSE]
        peaks$dist <- wald$dist[match(peaks$mark, wald$mark)]
        yrc <- ggplot_build(gp)$layout$panel_scales_y[[1]]$range$range
        peaks$wd <- ifelse(rep(LOD, nrow(peaks)), 0.5 * log(exp(peaks$wd), base = 10), peaks$wd)
        peaks$wd <- peaks$wd  + (diff(yrc)/20)
        peaks$mark <- gsub("X.", "", peaks$mark)
        gp <- gp + geom_text(data = peaks, aes_string(label = "mark"), size = 4)
    }
    dist.label <- tapply(wdc$dist, list(wdc$chr), function(el) (el[1] + el[length(el)])/2)
    ann.type <- unique(wald$type)[length(unique(wald$type))]
    gb <- ggplot_build(gp)$layout
#    mr <- diff(gb$panel_params[[1]]$y$minor_breaks[1:2])
    mr <- diff(gb$panel_scales_y[[1]]$range$range)/8
    yl <- - mr*(length(unique(wald$type))/3)
    ann <- cbind.data.frame(wd = rep(yl, nchr), dist = dist.label, type = ann.type)
    ann$chr <- rownames(ann)
    yrc <- gb$panel_scales_y[[1]]$range$range
    gp <- gp + geom_text(data = ann, aes_string(label = "chr"), size = 3.5, colour = "grey40") +
        coord_cartesian(ylim = c(0, yrc[2] + diff(yrc)/12), clip="off")
    gp
}


theme_scatter <- function (base_size = 11, base_family = "") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
        legend.position = "none",
        panel.grid.minor = element_line(colour = "grey90", size = 0.4),
        panel.grid.major = element_line(colour = "grey90", size = 0.8),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = base_size),
        axis.text.y = element_text(size = base_size),
        strip.text = element_text(size = base_size))
}

constructCM <- function(genoData, scale.method = "diag") {
        tg <- t(genoData)
        relm <- crossprod(tg)
        scale <- mean(diag(relm))
        relm <- relm/scale
        attr(relm, "rowNames") <- dimnames(relm)[[2]] <- rownames(genoData)
        ch <- chol(relm)
        chol.inv <- chol2inv(ch)
        rm.env <- new.env()
        rm.env$trans <- (tg %*% chol.inv)/scale
        rm.env$relm <- relm
        rm.env$scale <- scale
        rm.env
}

getQTL <- function (object, genObj)
{
    spe <- strsplit(object$QTL$effects, "\\.")
    wchr <- sapply(spe, "[", 2)
    wint <- as.numeric(sapply(spe, "[", 3))
    qtlm <- matrix(ncol = 6, nrow = length(wchr))
    for (i in 1:length(wchr)) {
        lhmark <- genObj$geno[[wchr[i]]]$map[wint[i]]
        qtlm[i, 1:4] <- c(wchr[i], wint[i], names(lhmark), round(lhmark,
            2))
        if (object$QTL$type == "interval") {
            if (length(genObj$geno[[wchr[i]]]$map) > 1)
                rhmark <- genObj$geno[[wchr[i]]]$map[wint[i] +
                  1]
            else rhmark <- genObj$geno[[wchr[i]]]$map[wint[i]]
            qtlm[i, 5:6] <- c(names(rhmark), round(rhmark, 2))
        }
        else qtlm <- qtlm[, -c(5:6)]
    }
    qtlm
}

summary.gwasreml <- function (object, genObj, LOD = TRUE, ...)
{
    if (missing(genObj))
        stop("genObj is a required argument")
    if (!inherits(genObj, "cross"))
        stop("genObj is not of class \"cross\"")
    if (is.null(object$QTL$peaks)) {
        cat("There are no significant putative QTL's above the threshold.\n")
        return(character(0))
    }
    trait <- object$QTL$Trait
    coefs <- object$coefficients$fixed
    inds <- grep("X\\.", rownames(coefs))
    nams <- rownames(coefs)[inds]
    coefs <- coefs[inds,]
    vcoef <- object$vcoeff$fixed[inds]
    sigma2 <- ifelse(object$QTL$section, 1, object$sigma2)
    zrat <- coefs/sqrt(vcoef*sigma2)
    enams <- strsplit(nams, ":")
    object$QTL$effects <- sapply(enams, function(el) el[grep("X\\.", el)])
    traits <- sapply(enams, function(el){
        if(length(el) > 1) el[-grep("X\\.", el)]
        else "MAIN"
    })
    prefix <- paste(trait, "_", sep = "")
    traits <- gsub(prefix, "", traits)
    qtlm <- as.data.frame(getQTL(object, genObj))
    if(object$QTL$type == "interval")
        names(qtlm) <- c("Chromosome", "Interval", "Left Marker", "dist(cM)", "Right Marker", "dist(cM)")
    else names(qtlm) <- c("Chromosome", "Interval", "Marker", "dist(cM)")
    qtlm <- cbind.data.frame(Env = traits, qtlm)
    qtlm$Size <- round(coefs, 4)
    qtlm$Pvalue <- round(2 * (1 - pnorm(abs(zrat))), 4)
    if(LOD) qtlm$LOD <- round(0.5*log(exp(zrat^2), base = 10), 4)
    nints <- as.numeric(as.character(qtlm$Interval))
    qtlm <- qtlm[order(qtlm$Chromosome, nints, qtlm$Env),]
    qtlm
}

waldTest.asreml <- function(object, cc, keep.fac = TRUE)
{
    if(oldClass(object) != "asreml")
        stop("Requires an object of class asreml\n")
    if(is.null(object$Cfixed)) {
        warning("Requires C matrix from model object. Refitting test model with argument \"Cfixed = T\"\n")
        asreml.options(Cfixed = TRUE)
        object <- update(object)
    }
    vrb <- as.matrix(object$Cfixed)
    tau <- c(object$coefficients$fixed)
    names(tau) <- rownames(object$coefficients$fixed)
    nc <- length(tau)
    sigma2 <- object$sigma2
    vrb <- vrb/sigma2
    ccnams <- names(tau)
    zdf <- cdf <- NULL
    cc <- lapply(cc, function(el, ccnams){
        if(!all(names(el) %in% c("coef","type","comp","group")))
            stop("Inappropriately named argument for comparison object.")
        if(is.numeric(el$coef)) {
            if(max(el$coef) > length(ccnams))
                stop("coefficient subscript out of bounds")
            names(el$coef) <- ccnams[el$coef]
        }
        else {
            if(any(is.na(pmatch(el$coef, ccnams))))
                  stop("Names in contrast do not match the names of coefficients of object")
            temp <- pmatch(el$coef, ccnams)
            names(temp) <- el$coef
            el$coef <- temp
        }
        el
    }, ccnams)
   ## split contrasts and other available tests
    ctype <- unlist(lapply(cc, function(el) el$type))
    if(!all(ctype %in% c("con","zero")))
        stop("Contrast types must be either \"con\" for treatment comparisons or \"zero\" for testing zero equality")
    cons <- cc[ctype %in% "con"]
    zero <- cc[ctype %in% "zero"]
    cse <- ctau <- zwtest <- cwtest <- zpval <- c()
    if(length(cons)) {
       CRows <- lapply(cons, function(el, nc){
           if(length(el) < 3){
               con <- contr.helmert(length(el$coef))[, (length(el$coef) - 1)]
               names(con) <- cnam <- names(el$coef)
               cat("Warning: default contrast being taken for", cnam, "is", con, "\n")
               row <- rep(0, nc)
               row[el$coef] <- con
               row
           }
           else {
               if(is.matrix(el$comp)) {
                   if(length(el$coef) != ncol(el$comp))
                       stop("Length of contrast does not match the number of specified coefficients")
                   cons <- split(el$comp, 1:nrow(el$comp))
                   rows <- lapply(cons, function(ell, first = el$coef, nc){
                       row <- rep(0, nc)
                       row[first] <- ell
                       row
                   }, first = el$coef, nc)
                   rows <- unlist(rows, use.names = F)
                   matrix(rows, nrow = nrow(el$comp), byrow = T)
               }
               else {
                   if(length(el$coef) != length(el$comp))
                       stop("Length of contrast does not match the number of specified coefficients")
                   row <- rep(0, nc)
                   row[el$coef] <- el$comp
                   row
               }
           }
       }, nc)
       Cmat <- do.call("rbind", CRows)
       if(!keep.fac)
           ccnams <- substring(ccnams, regexpr("\\_", ccnams) + 1, nchar(ccnams))
       cnam <- lapply(split(Cmat, 1:nrow(Cmat)), function(el, ccnams){
           namr <- ccnams[ifelse(el < 0, T, F)]
           naml <- ccnams[ifelse(el > 0, T, F)]
           c(paste(naml, collapse = ":"), paste(namr, collapse = ":"))
       }, ccnams)
       Cnam <- do.call("rbind", cnam)
       gnams <- lapply(cons, function(el){
           if(!is.null(el$group)){
               if(!any(names(el$group) %in% c("left","right")))
                   stop("group names must be \"left\" and \"right\".")
               if(is.null(el$group$left)){
                   if(is.matrix(el$comp))
                       el$group$left <- rep(NA, nrow(el$comp))
                   else el$group$left <- NA
               } else {
                   if(is.matrix(el$comp)){
                       if(length(el$group$left) == 1)
                           el$group$left <- rep(el$group$left, nrow(el$comp))
                       if(length(el$group$left) != nrow(el$comp))
                          stop("No. of group names do not match the number of comparisons in object")
                   }
               }
                if(is.null(el$group$right)){
                   if(is.matrix(el$comp))
                       el$group$right <- rep(NA, nrow(el$comp))
                   else el$group$right <- NA
               } else {
                   if(is.matrix(el$comp)) {
                       if(length(el$group$right) == 1)
                           el$group$right <- rep(el$group$right, nrow(el$comp))
                       if(length(el$group$right) != nrow(el$comp))
                          stop("No. of group names do not match the number of comparisons in object")
                   }
               }
           } else {
               if(is.matrix(el$comp))
                   el$group$left <- el$group$right <- rep(NA, nrow(el$comp))
               else el$group$left <- el$group$right <- NA
           }
           cbind(el$group$left, el$group$right)
       })
       Gnam <- do.call("rbind", gnams)
       Cnam[!is.na(Gnam[,1]), 1] <- Gnam[!is.na(Gnam[,1]),1]
       Cnam[!is.na(Gnam[,2]), 2] <- Gnam[!is.na(Gnam[,2]),2]
       for(i in 1:nrow(Cmat)) {
           varmat <- sum(Cmat[i,  ]*crossprod(vrb, t(Cmat)[, i]))
           cse[i] <- sqrt(varmat * sigma2)
           ctau[i] <- sum(Cmat[i,  ]*tau)
           cwtest[i] <- (ctau[i]/cse[i])^2
       }
       cdf <- data.frame(wald = round(cwtest, 6), pval = round(1 - pchisq(cwtest, 1), 6),
                         coef = round(ctau, 6), se = round(cse, 6))
       attr(cdf, "names") <- c("Wald Statistic", "P-Value", "Cont. Coef.", "Std. Error")
       attr(cdf, "row.names") <- paste(Cnam[,1], Cnam[,2],  sep = " vs ")
       oldClass(cdf) <- "data.frame"
   }
      if(length(zero)) {
        ZRows <- lapply(zero, function(el, nc){
            rows <- rep(rep(0, nc), length(el$coef))
            dum <- seq(0, (length(el$coef) - 1) * nc, by = nc)
            rows[el$coef + dum] <- 1
            matrix(rows, nrow = length(el$coef), byrow = T)
        }, nc)
        znam <- unlist(lapply(zero, function(el, ccnams) {
            if(is.null(el$group))
                paste(ccnams[el$coef], collapse = ":")
            else el$group
            }, ccnams))
        if(any(table(znam) > 1))
            stop("Duplicate names in group structures for zero equality tests.")
        for(i in 1:length(ZRows)) {
            varmat <- ZRows[[i]] %*% crossprod(vrb, t(ZRows[[i]]))
            Ctau <- ZRows[[i]] %*% tau
            zwtest[i] <- sum(Ctau*crossprod(solve(varmat), Ctau))/sigma2
            zpval[i] <- 1 - pchisq(zwtest[i], nrow(ZRows[[i]]))
        }
        zdf <- data.frame(wald = round(zwtest, 6), pval = round(zpval, 6))
        attr(zdf, "names") <- c("Wald Statistic", "P-Value")
        attr(zdf, "row.names") <- znam
        oldClass(zdf) <- "data.frame"
      }
    res <- list(Contrasts = cdf, Zero = zdf)
    invisible(res)
}

waldTest <- function(object, ...)
    UseMethod("waldTest")
