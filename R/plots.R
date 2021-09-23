.evalPrismaParts <- function(dn,qus,r2=c("r2_adj","r2_full"),
    test=c("empirical","wilcoxon","ttest")) {
    
    # Inputs
    r2 <- r2[1]
    test <- test[1]
    sel <- ifelse(r2=="r2_adj","r2p","r2m")
    
    # Plot/selection data
    baseline <- .getR2(dn)[,sel]
    qu <- lapply(qus,function(x,s) {
        return(x[,s])
    },sel)
    
    # Siginificance of the result
    switch(test,
        empirical = {
            p <- unlist(lapply(qu,function(x,b) {
                b <- .subsampleBase(x,b)
                return(length(which(x < b))/length(b))
            },baseline))
        },
        wilcoxon = {
            p <- unlist(lapply(qu,function(x,b) {
                if (length(b) > 2*length(x))
                    b <- .subsampleBase(x,b)
                h <- wilcox.test(x,b,alternative="greater")
                return(h$p.value)
            },baseline))
        },
        ttest = {
            p <- unlist(lapply(qu,function(x,b) {
                if (length(b) > 2*length(x))
                    b <- .subsampleBase(x,b)
                h <- t.test(x,b,alternative="greater")
                return(h$p.value)
            },baseline))
        }
    )
    
    gQuanDens <- .plotQuantileDensities(baseline,qu)
    gQuanSel <- .plotQuantileEvaluation(qu,p,baseline)
    
    return(list(density=gQuanDens,selection=gQuanSel))
}

.plotQuantileEvaluation <- function(qu,p,b) {
    r2 <- unlist(lapply(qu,mean))
    s <- unlist(lapply(qu,sd))
    p <- .zeroFix(p,.Machine$double.eps)
    logp <- -log10(p)
    if (max(logp) == 0)
        normP <- rep(0,length(logp))
    else
        normP <- logp/max(logp)
    names(r2) <- names(s) <- 
        paste0(100*suppressWarnings(as.numeric(gsub("Q","",names(qu)))),"%")

    ggdata <- data.frame(
        Index=seq_along(r2),
        PR2=r2,
        SD=s,
        Quantile=names(r2),
        P=.scaleP2R(normP,min(normP-s),max(normP+s),min(r2-s),max(r2+s)),
        Ptext=as.character(round(logp,3))
    )
     
    #breaksR <- pretty(c(ggdata$PR2-s,ggdata$PR2+s))
    if (max(b) > max(ggdata$PR2+s)) {
        breaksR <- pretty(c(ggdata$PR2-s,b+sd(b)))
        labelsR <- paste0(round(100*breaksR,digits=1),"%")
        breaksP <- pretty(c(ggdata$P,ggdata$P))
        tmp <- normP*max(logp)
        labelsP <- round(seq(min(tmp),max(tmp),length.out=length(breaksP)),
            digits=1)
    }
    else {
        breaksR <- pretty(c(b-sd(b),ggdata$PR2-s,ggdata$PR2+s,b+sd(b)))
        labelsR <- paste0(round(100*breaksR,digits=1),"%")
        breaksP <- pretty(c(ggdata$P,ggdata$P))
        tmp <- normP*max(logp)
        labelsP <- round(seq(min(tmp),max(tmp),length.out=length(breaksP)),
            digits=1)
    }
    
    
    g <- ggplot(ggdata,aes(x=Quantile,y=PR2)) + 
        geom_bar(stat="identity",fill="#2A9D8F",width=0.6) +
        geom_errorbar(aes(ymin=PR2-SD,ymax=PR2+SD),width=0.1,colour="#264653") +
        geom_hline(yintercept=mean(b),colour="#264653") +
        geom_hline(yintercept=mean(b)+sd(b),linetype="dashed",
            colour="#264653") +
        geom_hline(yintercept=mean(b)-sd(b),linetype="dashed",
            colour="#264653") +
        geom_point(aes(x=Index,y=P),size=3,colour="#E76F51") +
        geom_line(aes(x=Index,y=P),linetype="dashed",colour="#E76F51") +
        scale_y_continuous(name=expression("PRSice adjusted mean R"^"2"),
            breaks=breaksR,labels=labelsR,limits=c(min(breaksR),max(breaksR)),
            sec.axis=sec_axis(~ . *1,name=expression("-log"[10]*"(p-value)"),
            breaks=breaksP,labels=labelsP)) +
        #geom_text(aes(x=Index,y=P,label=Ptext),colour="#E76F51") +
        xlab("\nPRS candidates frequency quantile cutoff") +
        theme_bw() + 
        theme(
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14,face="bold"),
            axis.text.x=element_text(size=12,face="bold"),
            axis.text.y=element_text(size=12,face="bold",color="#2A9D8F"),
            axis.line.y=element_line(color="#2A9D8F",size=1),
            axis.ticks.y=element_line(color="#2A9D8F"),
            axis.text.y.right=element_text(color="#E76F51"),
            axis.line.y.right=element_line(color="#E76F51",size=1),
            axis.ticks.y.right=element_line(color="#E76F51")
        )
    g
    
    #ggsave(file=".png"),plot=g,width=8,height=6)
}

# Just for fun
.plotQuantileDensities <- function(base,qu) {
    ggdata <- data.frame(
        PR2=c(base,unlist(qu)),
        Source=c(rep("Baseline",length(base)),rep(names(qu),lengths(qu))),
        Type=factor(c(rep("Analysis",length(base)),rep("Selection",
            sum(lengths(qu)))),levels=c("Selection","Analysis"))
    )
    helper <- ggdata$Source
    helper <- paste0(100*suppressWarnings(as.numeric(gsub("Q","",helper))),"%")
    helper[helper=="NA%"] <- "Baseline"
    ggdata$Quantile <- helper

    g <- ggplot(ggdata,aes(x=PR2)) +
        geom_density(data=ggdata[ggdata$Type=="Selection",],
            aes(colour=Quantile)) +
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
    g
}

#https://stackoverflow.com/questions/5294955/
#   how-to-scale-down-a-range-of-numbers-with-a-known-min-and-max-value
.scaleP2R <- function(x,minp,maxp,minb,maxb) {
    return(((maxb-minb)*(x-minp))/(maxp-minp) + minb)
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
