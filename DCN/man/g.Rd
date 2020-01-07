\name{g}
\alias{g}
\docType{data}
\title{
	The differential co-expression network
}
\description{
	The differential co-expression network. It is an igraph object with 
	at least two vertex attributes ("name" and "DE") and one edge attribute ("DC").
}
\usage{data("g")}
\format{
	An \code{igraph} R object.
}
\details{
\code{g} is a differential co-expression network. There are 2074 nodes and 3917 edges in \code{g}. 
Each node in the graph represents a PCG/lncRNA. Each node has three attributes: "name", "DE", and "Type". 
"name" represents gene name of the node."DE" represents the differential expression score of the node. 
"Type" represents the type (PCG/lncRNA) of the gene. Each edge has one attribute "DC" which represents 
the differential co-expression score of the edge.
}

\examples{
	data(g)
}
\keyword{datasets}
