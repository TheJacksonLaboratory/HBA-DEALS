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

  ds<-(1:ngenes %in% sample(1:ngenes,num.dif))
  baselineprop1 <- baselineprop2 <- baselineprop
  baselineprop1[i1] <- baselineprop1[i1]*fc
  baselineprop2[i2] <- baselineprop2[i2]*fc
  mu0.1 <- matrix(baselineprop1,ngenes,1) %*% matrix(expected.lib.size[1:n1],1,n1)
  mu0.2 <- matrix(baselineprop2,ngenes,1) %*% matrix(expected.lib.size[(n1+1):(n1+n2)],1,n2)
  mu0 <- cbind(mu0.1,mu0.2)
  status <- rep(0,ngenes)
  status[i1] <- -1
  status[i2] <- 1

  # Biological variation
  BCV0 <- 0.2+1/sqrt(mu0)
  if(invChisq){
      df.BCV <- 40
      BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
    } else {
      BCV <- BCV0*exp( rnorm(ngenes,mean=0,sd=0.25)/2 )
  }

  if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
  shape <- 1/BCV^2
  scale <- mu0/shape

  mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)

  split.row<-function(gene.vec,is.control,num.iso,ds) #Apply it to mu, record gene ids
  {
    frac<-abs(rnorm(num.iso))
    frac<-frac/sum(frac)

    if (ds==FALSE){

      new.rows<-matrix(rep(gene.vec,num.iso),byrow = TRUE,nrow=num.iso)*frac
      ds.isoforms<-rep(1,num.iso)

    }else{
      num.diff.up<-sample(1:(num.iso-1),1)
      ds.up<-sample(c(rep(TRUE,num.diff.up),rep(FALSE,num.iso-num.diff.up)),num.iso)
      num.diff.down<-sample(1:sum(!ds.up),1)
      ds.down<-rep(FALSE,num.iso)
      if (sum(!ds.up)==1){ds.down[which(!ds.up)]=TRUE}else{ds.down[sample(which(!ds.up),num.diff.down)]<-TRUE}
      frac.case<-frac
      frac.case[ds.up]<-frac.case[ds.up]+sum(frac.case[ds.down]/2)/num.diff.up
      frac.case[ds.down]<-frac.case[ds.down]-frac.case[ds.down]/2
      new.rows<-cbind(matrix(rep(gene.vec[is.control],num.iso),byrow = TRUE,nrow=num.iso)*frac,matrix(rep(gene.vec[!is.control],num.iso),byrow = TRUE,nrow=num.iso)*frac.case)
      ds.isoforms<-frac.case/frac
      ds.isoforms[!(ds.up | ds.down)]<-1
    }

    return(list(new.rows,ds.isoforms))
  }


  new.mu<-NULL
  prev.rows<-0
  isoform.status<-c()
  for (i in (1:ngenes))
  {

    split.res<-split.row(mu[i,],c(rep(TRUE,n1),rep(FALSE,n2)),num.isoforms[i],ds[i])

    new.mu<-rbind(new.mu,split.res[[1]])

    isoform.status<-c(isoform.status,split.res[[2]])

    rownames(new.mu)[(prev.rows+1):nrow(new.mu)]<-paste0('Gene-',i)

    prev.rows<-nrow(new.mu)

  }

  # Technical variation
  counts <- matrix(rpois(sum(num.isoforms)*nlibs,lambda=new.mu),sum(num.isoforms),nlibs)
  rownames(counts)<-rownames(new.mu)

  # Filter
  keep <- rowSums(counts[,1:n1]>0)>=min.nonzero & rowSums(counts[,(n1+1):(n1+n2)]>0)>=min.nonzero & rowSums(counts)>=min.counts
  nkeep <- sum(keep)
  counts2 <- counts[keep,]
  isoform.status<-isoform.status[keep]
  keep<-rownames(counts2) %in% rownames(counts2)[duplicated(rownames(counts2))]
  counts2<-counts2[keep,]
  isoform.status<-isoform.status[keep]
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
  write.table(cbind(counts2[which(isoform.status!=1),1:2],isoform.status[which(isoform.status!=1)]),'ds_isoforms.txt',sep='\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
}
