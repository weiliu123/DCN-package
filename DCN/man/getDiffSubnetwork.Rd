\name{getDiffSubnetwork}
\alias{getDiffSubnetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get differential co-expression subnetwork
}
\description{
A differential subnetwork identification method for differential co-expression network.
}
\usage{
getDiffSubnetwork(graph, seed, d = 2, r = 0.05, alpha = 0.7, th = 0.5^6)
}

\arguments{
  \item{graph}{
	An igraph object with at least two vertex attributes ("name" and "DE")
	and one edge attribute ("DC").
}
  \item{seed}{
	The seed node. Should be the name of a vertex in \code{graph}.
}
  \item{d}{
	A positive integer that controls the searching space.
}
  \item{r}{
	A positive real number that controls the increasing rate of the subnetwork score.
}
  \item{alpha}{
	A real number between 0 and 1 that controls the relative weight of
	differential expression score and differential co-expression score.
}
  \item{th}{
	A threshold (> 0) of differential co-expression score. Edges with differential co-expression
	score < th will not be added into the differential subnetwork.
}
}
\details{
	This function implements a differential subnetwork searching algorithm which considers 
	both the differential expression score and the differential co-expression score. The search 
	starts from a seed and iteratively adds neighboring nodes. At each iteration, the search 
	considers addition of a node from the neighbors of nodes in the current subnetwork and 
	the corresponding edges connect this node and the current subnetwork. The addition which 
	yields maximum score is adopted. The search will stop if no node satisfies the following 
	two conditions: (i) the number of edges in the shortest path between this node and the seed 
	is less than or equal to d; and (ii) the addition of this node increase the score of the 
	subnetwork over an improvement rate r.
	The \code{graph} should have at least two vertex attributes ("name" and "DE") and one edge
	attribute ("DC"). "name" of a vertex represents the gene name of this vertex; "DE" of a vertex
	represents the t statistic in a paired, two-tailed t-test comparing the expression values
	of this vertex (gene) between tumor samples and normal samples; "DC" of an edge represents
	the weight of differential co-expression (differential co-expression score) between the 
	corresponding two vertexs.
}
\value{
	\item{subg }{The identified differential co-expression subnetwork.}
	\item{D }{The subnetwork score of \code{subg}.}
	\item{DE}{The differential expression score of \code{subg}.}
	\item{DC}{The differential co-expression score of \code{subg}.}
}
\references{
	Liu et al., Identification of lncRNA associated differential subnetworks in esophageal squamous 
	cell carcinoma by differential co-expression analysis, 2020.
}
\author{
	Wei Liu <freelw@qq.com>
}


\examples{
	library("igraph")
	data(g)
	res <- getDiffSubnetwork(graph = g, seed = "AL121899.1", 
                          d = 2, r = 0.1, alpha=0.7 ,th = 0.01)
	plot(res$subg)
}

