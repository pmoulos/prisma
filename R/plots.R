plotCvMetrics <- function(cvo,what=c("r2","rmse","mae","pr2","crl")) {
    gcv <- .plotCvMetricsCv(cvo,what)
    gpr <- .plotCvMetricsCv(.cvToPrsForPlot(cvo),what)
    return(list(cvBased=gcv,prsBased=gpr))
}

.plotCvMetricsCv <- function(cvo,what=c("r2","rmse","mae","pr2","crl")) {
    what <- what[1]
    labs <- .getPlotColsLabs(what)
    sel <- labs$sel
    nams <- labs$nams
    yl <- labs$yl
    
    # The numbers
    M <- unlist(lapply(cvo,function(x) {
        tmp <- lapply(x,function(y) {
            apply(y[,sel],2,mean) 
        })
        return(as.numeric(do.call("rbind",tmp)))
    }))
    SD <- unlist(lapply(cvo,function(x) {
        tmp <- lapply(x,function(y) {
            apply(y[,sel],2,sd) 
        })
        return(as.numeric(do.call("rbind",tmp)))
    }))
    
    # The labels
    nsnpLab <- rep(names(cvo),each=length(cvo[[1]])*length(sel))
    nsnpLab <- paste(nsnpLab,"SNPs")
    nsnpLab <- factor(nsnpLab,levels=paste(names(cvo),"SNPs"))
    
    metricLab <- rep(rep(nams,each=length(cvo[[1]])),length(cvo))
    metricLab <- factor(metricLab,levels=nams)
    
    index <- rep(rep(seq_along(cvo[[1]]),length(sel)),length(cvo))
    
    # The data
    ggdata <- data.frame(
        Mean=M,
        SD=SD,
        Metric=metricLab,
        NSNP=nsnpLab,
        Index=index
    )
    
    nc <- ifelse(length(cvMetrics[[1]]) > 5,2,3)
    
    g <- ggplot(ggdata,aes(x=Index,y=Mean,colour=Metric)) +
        geom_line() + 
        geom_point(size=2) +
        geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width=0.2,size=0.25) +
        facet_wrap(. ~ NSNP,ncol=nc) +
        scale_x_continuous(breaks=seq_along(cvo[[1]]),
            labels=paste0(100*as.numeric(names(cvo[[1]])),"%")) +
        xlab("\nPercentage of leave-out samples in cross-validation") +
        ylab(yl) +
        theme_bw() +
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.text.x=element_text(size=11,face="bold"),
            axis.text.y=element_text(size=12),
            strip.text.x=element_text(size=10,face="bold"),
            legend.position="bottom",
            legend.text=element_text(size=10),
            legend.title=element_text(size=12)
        )
    
    return(g)
}

.plotCvMetricsPrs <- function(cvm,what=c("r2","rmse","mae","pr2","crl")) {
    # Labels etc based on What to plot
    what <- what[1]
    labs <- .getPlotColsLabs(what)
    sel <- labs$sel
    nams <- labs$nams
    yl <- labs$yl
    
    # Numbers
    M <- unlist(lapply(cvm,function(x,s) {
        apply(x[,s],2,mean)
    },sel))
    S <- unlist(lapply(cvm,function(x,s) {
        apply(x[,s],2,sd)
    },sel))
    
    # Assemble the data
    ggdata <- data.frame(
        X=factor(rep(paste(names(cvm),"SNPs"),each=length(nams)),
            levels=paste(names(cvm),"SNPs")),
        Mean=M,
        SD=S,
        Source=factor(rep(nams,length(cvm)),levels=nams)
    )
    if (length(cvm) > 1) {
        ggdata$Type <- c(rep("Main",length(nams)),
            rep("Others",length(nams)*(length(cvm)-1)))
        g <- ggplot(ggdata,aes(x=X,y=Mean,fill=Source,linetype=Type))
    }
    else
        g <- ggplot(ggdata,aes(x=X,y=Mean,fill=Source))
        
    
    g <- g + geom_bar(stat="identity",position=position_dodge(width=0.6),
        width=0.7,colour="black") +
        geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),
            position=position_dodge(width=0.6),width=0.3,colour="black",
            linetype="solid") +
        xlab("\nNumber of SNPs in PRS candidates") +
        ylab(yl) +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.text.x=element_text(size=11,face="bold",angle=45,vjust=0.9,
                hjust=1),
            axis.text.y=element_text(size=12,face="bold")
        ) +
        guides(linetype=guide_legend(override.aes=list(fill=NA,col="black")))
    
    return(g)
}

.plotPrsEvaluation <- function(base,metrics,by="prs_r2",pval="p_ttest",
    stat="mean") {
    r2 <- metrics[,paste0(stat,"_",by)]
    s <- metrics[,paste0(ifelse(stat=="mean","sd","iqr"),"_",by)]
    p <- .zeroFix(metrics[,pval],.Machine$double.eps)
    fr <- metrics[,"freq"]
    N <- metrics[,"n_snp"]
    r2adj <- sqrt(r2/log(N))
    logp <- -log10(p)
    names(r2) <- as.character(metrics[,"n_snp"])
    helpNames <- paste(names(r2)," (",fr,")",sep="")

    ggdata <- data.frame(
        Index=seq_along(r2),
        PR2=r2,
        APR2=r2adj,
        SD=s,
        Significance=logp,
        NFSNP=factor(helpNames,levels=helpNames),
        N=metrics[,"n_snp"],
        Frequency=fr
    )
    
    if (stat=="mean") {
        yl1 <- expression("PRS regression adjusted mean R"^"2")
        yl2 <- expression("PRISMA regression adjusted mean R"^"2")
    }
    else {
        yl1 <- expression("PRS regression adjusted median R"^"2")
        yl2 <- expression("PRISMA regression adjusted median R"^"2")
    }
    
    # PRS R2 vs #SNPs and siginificance
    g1 <- ggplot(ggdata,aes(x=NFSNP,y=PR2)) + 
        geom_bar(aes(fill=Significance),stat="identity",width=0.7) +
        scale_fill_gradient2(low="#264653",mid="#E9C46A",high="#FD3025") +
        geom_errorbar(aes(ymin=PR2-SD,ymax=PR2+SD),width=0.2,colour="#264653") +
        geom_hline(yintercept=mean(base),colour="#727272",linetype="longdash",
            size=0.75) +
        geom_hline(yintercept=mean(base)+sd(base),linetype="dotted") +
        geom_hline(yintercept=mean(base)-sd(base),linetype="dotted") + 
        xlab("\nNumber of PRS candidate SNPs (# appearances)") +
        ylab(yl1) +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14,face="bold"),
            axis.text.x=element_text(size=8,face="bold",angle=45,vjust=0.9,
                hjust=1),
            axis.text.y=element_text(size=12,face="bold")
        )
    
    # Adjusted (our) PRS R2 vs #SNPs and siginificance
    g2 <- ggplot(ggdata,aes(x=NFSNP,y=APR2)) + 
        geom_point(aes(colour=Significance),size=3) +
        #scale_colour_gradient2(low="#E00000",mid="#FFAA28",high="#00E000") +
        scale_colour_gradient2(low="#264653",mid="#E9C46A",high="#FD3025") +
        xlab("\nNumber of PRS candidate SNPs (# appearances)") +
        ylab(yl2) +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14,face="bold"),
            axis.text.x=element_text(size=8,face="bold",angle=45,vjust=0.9,     
                hjust=1),
            axis.text.y=element_text(size=12,face="bold")
        )
    
    # PRS R2 vs #SNPs and frequency
    g3 <- ggplot(ggdata,aes(x=N,y=PR2)) + 
        geom_point(aes(colour=Frequency),size=3) +
        scale_x_log10(labels=comma) +
        xlab("\nNumber of SNPs in PRS") +
        ylab(yl1) +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14,face="bold"),
            axis.text.x=element_text(size=12,face="bold",angle=45,vjust=0.9,     
                hjust=1),
            axis.text.y=element_text(size=12,face="bold")
        )
    
    return(list(r2snp=g1,ar2snp=g2,r2fr=g3))
}

.plotPrsTrait <- function(x,y,n) {
    ggdata <- data.frame(PRS=x,Y=y)
    
    p1 <- ggplot(ggdata,aes(x=PRS,y=Y)) + 
        geom_point(aes(colour=PRS),size=2) +
        scale_colour_gradient(low="#00A000",high="#A00000") +
        xlab("\nPolygenic Risk Score") +
        ylab(n) +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.text.x=element_text(size=12,face="bold"),
            axis.text.y=element_text(size=12,face="bold")
        )
    
    p2 <- ggplot(ggdata,aes(x=PRS,fill=..x..)) + 
        geom_histogram(color="black",size=0.5) +
        scale_fill_gradient(low="#00A000",high="#A00000",name="Risk") +
        xlab("\nPolygenic Risk Score") +
        ylab("Frequency\n") +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.text.x=element_text(size=12,face="bold"),
            axis.text.y=element_text(size=12,face="bold"),
            legend.title=element_text(size=11,face="bold"),
            legend.text=element_text(size=10)
        )
    
    return(list(scatter=p1,hist=p2))
}

# Just for fun
.plotFreqDensities <- function(base,metrics,by="prs_r2") {
    qu <- lapply(metrics,function(x) return(x[,by]))
    
    if (length(qu[[1]] < length(base)))
        base <- .subsampleBase(qu[[1]],base)
    
    ggdata <- data.frame(
        PR2=c(base,unlist(qu)),
        Source=factor(c(rep("Baseline",length(bnase)),rep(names(qu),
            lengths(qu))),levels=c("Baseline",names(qu))),
        Type=factor(c(rep("Analysis",length(base)),rep("Selection",
            sum(lengths(qu)))),levels=c("Selection","Analysis"))
    )

    g <- ggplot(ggdata,aes(x=PR2)) +
        geom_density(data=ggdata[ggdata$Type=="Selection",],
            aes(colour=Source)) +
        geom_density(data=ggdata[ggdata$Type=="Analysis",],aes(colour=Source),
            colour="black",linetype="longdash",size=0.8) +
        xlim(-0.05,round(max(ggdata$PR2),digits=1)+0.05) +
        xlab("PRSice adjusted R2") +
        ylab("Density") +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.text.x=element_text(size=12,face="bold"),
            axis.text.y=element_text(size=12,face="bold")
        )
    return(g)
}

.getPlotColsLabs <- function(what) {
    switch(what,
        r2 = { 
            sel <- c("full_r2","reduced_r2","prs_r2") 
            nams <- c("Full R2","Reduced R2","PRS R2")
            yl <- expression("Mean R"^"2")
        },
        rmse = { 
            sel <- c("full_rmse","reduced_rmse")
            nams <- c("Full RMSE","Reduced RMSE")
            yl <- "Root Mean Square Error"
        },
        mae = { 
            sel <- c("full_mae","reduced_mae")
            nams <- c("Full MAE","Reduced MAE")
            yl <- "Mean Absolute Error"
        },
        crl = {
            sel <- c("full_pred_cor","reduced_pred_cor") 
            nams <- c("Full R","Reduced R")
            yl <- "Correlation"
        },
        pr2 = {
            sel <- c("full_pred_r2","reduced_pred_r2","prs_pred_r2") 
            nams <- c("Full PR2","Reduced PR2","PRS PR2")
            yl <- "Adjusted correlation"
        }
    )
    return(list(sel=sel,nams=nams,yl=yl))
}

.subsampleBase <- function(x,b) {
    nx <- length(x)
    nb <- length(b)
    if (nx < nb) {
        ii <- sample(nb,nx)
        return(b[ii])
    }
    return(b)
}

#~ #https://stackoverflow.com/questions/5294955/
#~ #   how-to-scale-down-a-range-of-numbers-with-a-known-min-and-max-value
#~ .scaleP2R <- function(x,minp,maxp,minb,maxb) {
#~     return(((maxb-minb)*(x-minp))/(maxp-minp) + minb)
#~ }

#~ .plotQuantileEvaluation <- function(qu,p,b) {
#~  qu <- lapply(metrics,function(x) return(x[,by]))
#~     r2 <- unlist(lapply(qu,mean))
#~     s <- unlist(lapply(qu,sd))
#~     p <- .zeroFix(p,.Machine$double.eps)
#~     logp <- -log10(p)
#~     if (max(logp) == 0)
#~         normP <- rep(0,length(logp))
#~     else
#~         normP <- logp/max(logp)
#~     names(r2) <- names(s) <- 
#~         paste0(100*suppressWarnings(as.numeric(gsub("Q","",names(qu)))),"%")

#~     ggdata <- data.frame(
#~         Index=seq_along(r2),
#~         PR2=r2,
#~         SD=s,
#~         Quantile=names(r2),
#~         P=.scaleP2R(normP,min(normP-s),max(normP+s),min(r2-s),max(r2+s)),
#~         Ptext=as.character(round(logp,3))
#~     )
     
#~     if (max(b) > max(ggdata$PR2+s)) {
#~         breaksR <- pretty(c(ggdata$PR2-s,b+sd(b)))
#~         labelsR <- paste0(round(100*breaksR,digits=1),"%")
#~         breaksP <- pretty(c(ggdata$P,ggdata$P))
#~         tmp <- normP*max(logp)
#~         labelsP <- round(seq(min(tmp),max(tmp),length.out=length(breaksP)),
#~             digits=1)
#~     }
#~     else {
#~         breaksR <- pretty(c(b-sd(b),ggdata$PR2-s,ggdata$PR2+s,b+sd(b)))
#~         labelsR <- paste0(round(100*breaksR,digits=1),"%")
#~         breaksP <- pretty(c(ggdata$P,ggdata$P))
#~         tmp <- normP*max(logp)
#~         labelsP <- round(seq(min(tmp),max(tmp),length.out=length(breaksP)),
#~             digits=1)
#~     }
    
    
#~     g <- ggplot(ggdata,aes(x=Quantile,y=PR2)) + 
#~         geom_bar(stat="identity",fill="#2A9D8F",width=0.6) +
#~         geom_errorbar(aes(ymin=PR2-SD,ymax=PR2+SD),width=0.1,colour="#264653") +
#~         geom_hline(yintercept=mean(b),colour="#264653") +
#~         geom_hline(yintercept=mean(b)+sd(b),linetype="dashed",
#~             colour="#264653") +
#~         geom_hline(yintercept=mean(b)-sd(b),linetype="dashed",
#~             colour="#264653") +
#~         geom_point(aes(x=Index,y=P),size=3,colour="#E76F51") +
#~         geom_line(aes(x=Index,y=P),linetype="dashed",colour="#E76F51") +
#~         scale_y_continuous(name=expression("PRSice adjusted mean R"^"2"),
#~             breaks=breaksR,labels=labelsR,limits=c(min(breaksR),max(breaksR)),
#~             sec.axis=sec_axis(~ . *1,name=expression("-log"[10]*"(p-value)"),
#~             breaks=breaksP,labels=labelsP)) +
#~         #geom_text(aes(x=Index,y=P,label=Ptext),colour="#E76F51") +
#~         xlab("\nPRS candidates frequency quantile cutoff") +
#~         theme_bw() + 
#~         theme(
#~             axis.title.x=element_text(size=14),
#~             axis.title.y=element_text(size=14,face="bold"),
#~             axis.text.x=element_text(size=12,face="bold"),
#~             axis.text.y=element_text(size=12,face="bold",color="#2A9D8F"),
#~             axis.line.y=element_line(color="#2A9D8F",size=1),
#~             axis.ticks.y=element_line(color="#2A9D8F"),
#~             axis.text.y.right=element_text(color="#E76F51"),
#~             axis.line.y.right=element_line(color="#E76F51",size=1),
#~             axis.ticks.y.right=element_line(color="#E76F51")
#~         )
#~     g
#~ }
