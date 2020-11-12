\name{hbadeals}
\alias{hbadeals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical Bayesian Analysis of Differential Expression and ALternative SPlicing (HBA-DEALS)
}
\description{
The function determines differential expression and transcript proportions between two experimental conditions.  It takes as input a table that contains the gene names in its first column, the transcript names in its second column,
and the corresponding transcript counts under every sample in the rest of its columns.  In addition, it takes as input an ordered numeric vector of 1's and 2's, where 1 indicates that the corresponding sample belongs to the first experimental condition and 2 to the second. The function returns a table that provides 1-probability values for differential expression of each gene, and either differential splicing for the gene or differential proportion for each transcript.
}
\usage{
hbadeals(countsData,labels,n.cores,isoform.level,mcmc.iter=3000,mcmc.warmup=4000,hierarchy='auto',lib.size=NULL,mtc=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countsData}{
   A table of gene names, transcript names and transcript counts in each sample.  At least two transcripts must correspond to
   each gene.
}
\item{labels}{
  An ordered vector of 1's and 2's.  Its length is ncol(countsdata)-2.  Each entry indicates whether the corresponding sample/column
  of countsData belongs to the first experimental condition or the second.
}
\item{n.cores}{
  The number of cores to use in the calculation.  It is recommended to dedicate as many cores as possible.
}
\item{isoform.level}{
  Should 1-probability of differential proportion for each transcript be returned.  FALSE by deafult, which returns 1-probability
  of differential splicing for each gene.
}
\item{mcmc.iter}{
  The number of iterations of the MCMC algorithm after warmup.
}
\item{mcmc.warmup}{
  The number of warmup iterations of the MCMC algorithm.
}
\item{hierarchy}{
  Determines whether a hierarchical model will be used (hierarchy='yes'), a flat one (hierarchy='no') or will the decision
  will be made automatically (hierarchy='auto')
}
\item{lib.size}{
  A numeric vector containing total library sizes for each sample. If not provided, the default is  columnwise count totals.
}
\item{mtc}{
  A logical argument (default FALSE) that indicates whether the output probabolities should be corrected for multiple testing.
}

}
\value{
The table returned by this function contains 4 columns.  The first column is the gene name, the second is the transcript name,
the third is the fold change, and the fourth is 1-probability of differential expression or proportion.  Entries that
refer to expression have 'Expression' in their second column.  If isoform.level is FALSE, entries that refer to differential splicing of the gene will have 'Splicing' in their second column entry.  The fold change for expression is given as log2 fold change, and for splicing as fold change.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Guy Karlebach, Peter N Robinson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
simulate(rseed=1,fc=3,equal = TRUE)
n.samples=4 #per condition
countsData=read.table('counts.txt',sep='\t',header=FALSE)
labels=c(rep(1,n.samples),rep(2,n.samples))
res=hbadeals(countsData,labels)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expression }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~splicing }% __ONLY ONE__ keyword per line
