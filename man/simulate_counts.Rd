\name{simulate}
\alias{simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for simulating counts of transcripts.
}
\description{
This function prints a table of genes, their corresponding transcripts and simulated counts in multiple sample
and two different experimental conditions.  It also outputs tables of the differentially expressed genes and the transcripts whose proportions differ between the two conditions.
}
\usage{
simulate(ngenes=10000,invChisq=TRUE,fc=2,equal=TRUE,n1=4,n2=4,base.lib.size=11e6,rseed=2013,num.dif=200,min.nonzero=min(n1,n2),min.counts=10)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ngenes}{The number of genes to simulate.}
  \item{invChisq}{Whether to use the inverse chi^2 distribution for dispersion.  The alternative is log normal.}
  \item{fc}{Fold change of differentially expressed genes.}
  \item{equal}{Whether to simulate equal library sizes.  The alternative is fold 10 difference in library size.}
  \item{n1}{Number of sample in the first experimental condition.}
  \item{n2}{Number of sample in the second experimental condition.}
  \item{base.lib.size}{Base library size to simulate.}
  \item{rseed}{Seed for the simulation.  The same data will be generated in multiple calls with the same rseed.}
  \item{num.dif}{Num of differentially expressed genes and genes that have differential transcript proportions. }
  \item{min.nonzero}{Minimal number of non-zero entries in each condition that are necessary to keep a row in the generated table.}
  \item{min.counts}{Minimal number counts are necessary to keep a row in the generated table.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function prints a counts table 'counts.txt', differential genes table 'de_genes.txt' and differential transcripts table
'ds_isoforms.txt' to the working directory.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Guy Karlebach and Peter N Robinson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 simulate(rseed=1,fc=3,equal=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~simulation }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~expressions }% __ONLY ONE__ keyword per line
