\name{DCN-package}
\alias{DCN-package}
\alias{DCN}
\docType{package}
\title{
	Identifying differential subnetworks
}
\description{
	Identification of differential subnetworks from a differential co-expression network.
}
\author{
\packageAuthor{DCN}

Maintainer: \packageMaintainer{DCN}
}
\references{
	Liu et al., Identification of lncRNA associated differential subnetworks in esophageal squamous 
	cell carcinoma by differential co-expression analysis, 2020.
}

\keyword{ package }
\examples{
	library("igraph")
	data(g)
	res <- getDiffSubnetwork(graph = g, seed = "AL121899.1", 
                          d = 2, r = 0.1, alpha=0.7 ,th = 0.01)
	plot(res$subg)
}
