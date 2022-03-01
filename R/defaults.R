getDefaults <- function(what) {
    allowed <- c("filters","glm","rrblup","statgen","snptest","plink",
        "externalTools","lassosum","prsice")
    
    if (!(what %in% allowed))
        stop("what must be one of ",paste(allowed,collapse=", "))
    
    switch(what,
        filters = {
            return(list(
                snpCallRate=0.98,
                sampleCallRate=0.95,
                maf=0.05,
                hwe=1e-6,
                heteroStat="median",
                heteroFac=3,
                heteroHard=NA,
                pcaOut=TRUE,
                pcaRobust="hubert",
                nPC=NA,
                LD=0.2,
                IBD=0.1,
                inbreed=0.1
            ))
        },
        glm = {
            return(list(
                #family="gaussian"
                family=NULL,
                size=1
            ))
        },
        rrblup = {
            return(list(
                pcblup="auto",
                npcs=NULL
            ))
        },
        statgen = {
            return(list(
                kinship="astle",
                reml="EMMA",
                gls="single"
            ))            
        },
        snptest = {
            return(list(
                test="frequentist",
                model="additive",
                workspace=NULL
            ))
        },
        plink = {
            return(list(
                effect="genotypic",
                seed=42,
                workspace=NULL,
                cliopts=NULL
            ))
        },
        externalTools = {
            return(list(
                snptest="v2.5.6",
                plink="v1.90",
                prsice="v2.3.3",
                impute="v2.3.2",
                gtool="v0.7.5"
            ))
        },
        lassosum = {
            return(list(
                anc="eur",
                valid="auto"
            ))
        },
        prsice = {
            return(list(
                clump_kb=250,
                clump_r2=0.1,
                clump_p=1,
                score="avg",
                perm=10000,
                seed=42
            ))
        }
    )
}

.getPrismaMainDefaults <- function() {
    return(list(
        gwe=NULL,
        phenotype=NULL,
        covariates=NULL,
        pcs=FALSE,
        npcs=0,
        trainSize=0.8,
        niter=10,
        resolution="frequency",
        step=1,
        minFreq=2,
        minSnps=5,
        dropSameQuantiles=TRUE,
        aggregation="intersection",
        effectWeight="mean",
        prsSelectMethod="maxima",
        prsSelectCrit="prs_r2",
        prsSelectStat="mean",
        prsSelectR2="adjusted",
        filters=getDefaults("filters"),
        pcaMethod="snprel",
        imputeMissing=FALSE,
        imputeMethod="single",
        gwaMethods="glm",
        gwaCombine="simes",
        glmOpts=getDefaults("glm"),
        rrblupOpts=getDefaults("rrblup"),
        statgenOpts=getDefaults("statgen"),
        snptestOpts=getDefaults("snptest"),
        plinkOpts=getDefaults("plink"),
        prsMethods=c("lassosum","prsice"),
        lassosumOpts=getDefaults("lassosum"),
        prsiceOpts=getDefaults("prsice"),
        prsWorkspace=NULL,
        cleanup="intermediate",
        logging="file",
        dnOutput="gwaslist",
        output="normal",
        continue=FALSE,
        useDenovoWorkspace=NULL,
        runId=NULL,
        evalOnSplit="original",
        evalWith="vanilla",
        rc=NULL
    ))
}

.getPrsSelectionDefaults <- function() {
    return(list(
        dnList=NULL,
        gwe=NULL,
        phenotype=NULL,
        covariates=NULL,
        pcs=FALSE,
        npcs=0,
        trainSize=0.8,
        niter=10,
        resolution="frequency",
        step=1,
        minFreq=2,
        minSnps=5,
        dropSameQuantiles=TRUE,
        aggregation="intersection",
        effectWeight="mean",
        filters=getDefaults("filters"),
        pcaMethod="snprel",
        imputeMissing=FALSE,
        imputeMethod="single",
        gwaMethods="glm",
        gwaCombine="simes",
        glmOpts=getDefaults("glm"),
        rrblupOpts=getDefaults("rrblup"),
        statgenOpts=getDefaults("statgen"),
        snptestOpts=getDefaults("snptest"),
        plinkOpts=getDefaults("plink"),
        prsMethods=c("lassosum","prsice"),
        lassosumOpts=getDefaults("lassosum"),
        prsiceOpts=getDefaults("prsice"),
        prsWorkspace=NULL,
        cleanup="intermediate",
        logging="file",
        output="gwaslist",
        continue=FALSE,
        useDenovoWorkspace=NULL,
        runId=NULL,
        evalWith="vanilla",
        rc=NULL
    ))
}

.getPrsPipelineDefaults <- function() {
    return(list(
        gwe=NULL,
        phenotype=NULL,
        covariates=NULL,
        pcs=FALSE,
        npcs=0,
        snpSelection=NULL,
        trainSize=0.8,
        niter=10,
        filters=getDefaults("filters"),
        pcaMethod="snprel",
        imputeMissing=FALSE,
        imputeMethod="single",
        gwaMethods="glm",
        gwaCombine="simes",
        glmOpts=getDefaults("glm"),
        rrblupOpts=getDefaults("rrblup"),
        statgenOpts=getDefaults("statgen"),
        snptestOpts=getDefaults("snptest"),
        plinkOpts=getDefaults("plink"),
        prsMethods=c("lassosum","prsice"),
        lassosumOpts=getDefaults("lassosum"),
        prsiceOpts=getDefaults("prsice"),
        prsWorkspace=NULL,
        cleanup="intermediate",
        logging="screen",
        output="gwaslist",
        continue=FALSE,
        useDenovoWorkspace=NULL,
        runId=NULL,
        evalWith="vanilla",
        rc=NULL
    ))
}
