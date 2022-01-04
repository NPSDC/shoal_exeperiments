suppressPackageStartupMessages(library(tximeta))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(fishpond))
suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(SummarizedExperiment))

edgerPrep <- function(cts, condition = NULL, cov = NULL) {
    y <- DGEList(cts)
    y <- y[filterByExpr(y),]
    y <- calcNormFactors(y)
    if (is.null(cov)){
        design <- model.matrix(~condition)
    } else {
        design <- model.matrix(~condition + cov)
    }
    list(y, design)
}

limmavoom <- function(cts, condition = NULL, cov = NULL) {
    out <- edgerPrep(cts, condition = condition, cov = cov)
    y <- out[[1]]; design <- out[[2]]
    v <- voom(y, design)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    topTable(fit, coef=2, n=nrow(y), sort="none")
}

deseq2 <- function(cts, condition = NULL, cov = NULL) {
    cts <- round(cts)
    suppressMessages({
        if (is.null(cov)) {
            samps <- data.frame(condition=condition)
            dds <- DESeqDataSetFromMatrix(cts, samps, ~condition)
        } else {
            samps <- data.frame(condition=condition, cov=cov)
            dds <- DESeqDataSetFromMatrix(cts, samps, ~condition + cov)
        }
        dds <- DESeq(dds, minReplicates=Inf, quiet=TRUE)
    })
    res <- results(dds, name="condition_2_vs_1")
    res$padj[is.na(res$padj)] <- 1
    res
}

myicobra <- function(cd, lvl="Gene", n.sub, thrs = c(.01,.05,.1)) {
    cp <- calculate_performance(cd,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve", "roc"),
                                thrs=thrs
                                )
    cobraplot <- prepare_data_for_plot(cp,
                                       #colorscheme=cols,
                                       facetted = TRUE,
                                       incloverall = TRUE)
    yrng <- c(.4, 1)
    xrng <- c(0,max(.2, max(fdrtpr(cp)$FDR)))
    p <- plot_fdrtprcurve(cobraplot, plottype="points",
                     pointsize=2.5,
                     xaxisrange=xrng,
                     yaxisrange=yrng,
                     title=paste0(lvl,"-level, n=",n.sub," vs ",n.sub))
    return(list(cp, p))
}

getQuantFiles <- function(dir, type="shoal") {
    samples <- paste(outer(c(1:6), c(1:2), function(x,y) paste(x,y,sep="_")))
    if(!type %in% c("shoal", "salmon"))
        stop(paste(type, "is, expected shoal or salmon"))
    if(type == "salmon")
        files <- file.path(dir, samples, "quant.sf")
    if(type == "shoal")
        files <- file.path(dir, samples, "adapt_quant.sf")
    
    colData <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
    return(colData)
}

getCounts <- function(dfSalFiles, dfShoalFiles) {
    se <- tximeta::tximeta(dfSalFiles)
    txi <- tximport(dfSalFiles[["files"]], type="salmon", txOut=T,
                    countsFromAbundance="lengthScaledTPM")
    y <- labelKeep(se)
    y <- y[mcols(y)[["keep"]],] 
    
    seShoal <- se
    txiShoal <- tximport(dfShoalFiles[["files"]], type="salmon", txOut=T,
                    countsFromAbundance="lengthScaledTPM")
    lShoal <- lapply(dfShoalFiles[["files"]], function(f) read.delim(f,row.names = 1))
    cShoal <- matrix(0, dim(lShoal[[1]])[1], length(lShoal), dimnames = list(rownames(lShoal[[1]]),
            dfShoalFiles[["names"]]))
    for(i in seq_along(lShoal))
        cShoal[,i] <- lShoal[[i]][["NumReads"]]
    assays(seShoal)[["counts"]] <- cShoal
    yShoal <- labelKeep(seShoal)
    yShoal <- yShoal[mcols(yShoal)[["keep"]],] 
    return(list(txiSal = txi, ySal = y, txiShoal = txiShoal, yShoal = yShoal))
}

runlimma <- function(countList, thrs) {
    allTxps <- union(rownames(countList[["ySal"]]), rownames(countList[["yShoal"]]))
    dSal <- limmavoom(countList[["txiSal"]][["counts"]][allTxps,], countList[["ySal"]][["condition"]])
    dShoal <- limmavoom(countList[["txiShoal"]][["counts"]][allTxps,], countList[["yShoal"]][["condition"]])
    
    unTxps <- union(rownames(dSal), rownames(dShoal))
    missSal <- setdiff(unTxps, rownames(dSal))
    dSalU <- rbind(dSal, matrix(1, nrow = length(missSal), ncol = ncol(dSal),
                                dimnames = list(missSal, colnames(dSal))))
    missShoal <- setdiff(unTxps, rownames(dShoal))
    dShoalU <- rbind(dShoal, matrix(1, nrow = length(missShoal), ncol = ncol(dShoal),
                                    dimnames = list(missShoal, colnames(dShoal))))
    dShoalU <- dShoalU[unTxps,]
    dSalU <- dSalU[unTxps,]
    
    padj <- data.frame(row.names = rownames(dSalU))
    iso.any <- iso.dtu | iso.dte | iso.dge
    status <- as.numeric(rownames(padj) %in% names(iso.dtu)[iso.any])
    truth <- data.frame(status, row.names=rownames(padj))
    
    padj[["Sal"]] <- dSalU[["adj.P.Val"]]
    padj[["Shoal"]] <- dShoalU[["adj.P.Val"]]
    
    cd <- COBRAData(padj=padj, truth=truth)
    return(list(cd=cd, cob=myicobra(cd, "Transcript", 6, thrs=thrs), salLimma = dSalU, shoalLimma = dShoalU))
}

runDE <- function(dirSal, dirShoal, thrs = c(0.01, 0.05, 0.1)) {
    dfSal <- getQuantFiles(dirSal, type="salmon")
    dfShoal <- getQuantFiles(dirShoal, type="shoal")
    
    countsList <- getCounts(dfSal, dfShoal)
    cob <- runlimma(countsList, thrs)
    return(list(counts=countsList, cob=cob))
}

computeFactor <- function(cts, eps =1e-05) {
    rat <- rowSds(cts)/(rowMeans(cts))
    ifelse(is.nan(rat), 0, exp(-rat))
}