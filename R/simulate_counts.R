.onLoad <- function(libname, pkgname) {
    utils::data("voom_sim", package=pkgname, envir=parent.env(environment()))
}


simulate<-function(ngenes=10000,invChisq=TRUE,fc=2,equal=TRUE,n1=4,n2=4,base.lib.size=11e6,rseed=2013,num.dif=200,min.nonzero=min(n1,n2),min.counts=10)
{
    
    set.seed(rseed)
    
    num.isoforms<-sample(size=ngenes,x = 2:10,prob = c(0.4,0.2,0.1,0.05,0.05,0.05,0.05,0.05,0.05),replace = TRUE)
    
    baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
    baselineprop <- baselineprop/sum(baselineprop)
    
    # Design
    group <- factor(c(rep(1,n1),rep(2,n2)))
    design <- model.matrix(~group)
    
    nlibs <- n1+n2
    
    # Library size
    if(equal){
        expected.lib.size <- rep(base.lib.size,n1+n2)
    } else {
        expected.lib.size <- base.lib.size*2 * sample(c(1,0.1),n1+n2,replace = TRUE)
    }
    
    # Expected counts, group basis
    i <- sample(1:ngenes,num.dif)
    i1 <- i[1:(num.dif/2)]
    i2 <- i[(num.dif/2+1):num.dif]
    
    ds<-sample(1:ngenes,num.dif)
    baselineprop1 <- baselineprop2 <- baselineprop
    baselineprop1[i1] <- baselineprop1[i1]*fc
    baselineprop2[i2] <- baselineprop2[i2]*fc
    baselineprop1=unlist(lapply(1:length(num.isoforms),function(i,baselineprop1,num.isoforms){rep(baselineprop1[i],num.isoforms[i])},baselineprop1,num.isoforms))
    baselineprop1=baselineprop1/unlist(lapply(num.isoforms,function(x){rep(x,x)}))
    baselineprop2=unlist(lapply(1:length(num.isoforms),function(i,baselineprop2,num.isoforms){rep(baselineprop2[i],num.isoforms[i])},baselineprop2,num.isoforms))
    baselineprop2=baselineprop2/unlist(lapply(num.isoforms,function(x){rep(x,x)}))
    ds.isoforms=rep(1,sum(num.isoforms))
    for (ds.gene.idx in ds)
    {
        idxs.in.prop=sum(num.isoforms[(1:length(num.isoforms))<ds.gene.idx])+1
        idxs.in.prop=idxs.in.prop:(idxs.in.prop+num.isoforms[ds.gene.idx]-1)
        ds.case=sample(1:num.isoforms[ds.gene.idx],1)
        ds.control=sample(setdiff(1:num.isoforms[ds.gene.idx],ds.case),1)
        baselineprop1[idxs.in.prop[ds.case]]=baselineprop1[idxs.in.prop[ds.case]]*2
        baselineprop2[idxs.in.prop[ds.control]]=baselineprop2[idxs.in.prop[ds.control]]*2
        ds.isoforms[c(idxs.in.prop[ds.control],idxs.in.prop[ds.case])]=c(0.5,2)
        
    }
    
    mu0.1 <- matrix(baselineprop1,sum(num.isoforms),1) %*% matrix(expected.lib.size[1:n1],1,n1)
    mu0.2 <- matrix(baselineprop2,sum(num.isoforms),1) %*% matrix(expected.lib.size[(n1+1):(n1+n2)],1,n2)
    mu0 <- cbind(mu0.1,mu0.2)
    status <- rep(0,ngenes)
    status[i1] <- -1
    status[i2] <- 1
    
    # Biological variation
    BCV0 <- 0.2+1/sqrt(mu0)
    if(invChisq){
        df.BCV <- 40
        BCV <- BCV0*sqrt(df.BCV/rchisq(sum(num.isoforms),df=df.BCV))
    } else {
        BCV <- BCV0*exp( rnorm(sum(num.isoforms),mean=0,sd=0.25)/2 )
    }
    
    if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
    shape <- 1/BCV^2
    scale <- mu0/shape
    
    mu <- matrix(rgamma(sum(num.isoforms)*nlibs,shape=shape,scale=scale),sum(num.isoforms),nlibs)
    
    # Technical variation
    counts <- matrix(rpois(sum(num.isoforms)*nlibs,lambda=mu),sum(num.isoforms),nlibs)
    rownames(counts)<-paste0('Gene-',unlist(lapply(1:length(num.isoforms),function(i,num.isoforms){rep(i,num.isoforms[i])},num.isoforms)))
    
    # Filter
    keep <- rowSums(counts[,1:n1]>0)>=min.nonzero & rowSums(counts[,(n1+1):(n1+n2)]>0)>=min.nonzero & rowSums(counts)>=min.counts
    nkeep <- sum(keep)
    counts2 <- counts[keep,]
    ds.isoforms<-ds.isoforms[keep]
    keep<-rownames(counts2) %in% rownames(counts2)[duplicated(rownames(counts2))]
    counts2<-counts2[keep,]
    ds.isoforms<-ds.isoforms[keep]
    # Add gene name and isoform name columns
    
    counts2<-cbind(rownames(counts2),1:nrow(counts2),counts2)
    
    next.gene<-counts2[1,1]
    iso.itr<-1
    
    for (i in (1:nrow(counts2)))
    {
        
        if (counts2[i,1]!=next.gene)
        {
            iso.itr<-1
            next.gene<-counts2[i,1]
        }
        
        counts2[i,2]<-paste0('iso-',iso.itr)
        iso.itr<-iso.itr+1
        
    }
    
    write.table(counts2,'counts.txt',sep='\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(paste0('Gene-',which(status!=0)),'de_genes.txt',quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(cbind(counts2[which(ds.isoforms!=1),1:2],ds.isoforms[which(ds.isoforms!=1)]),'ds_isoforms.txt',sep='\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
}

