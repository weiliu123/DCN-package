\name{getDC}
\alias{getDC}
\title{
	Get differential co-expression score.
}
\description{
	Get differential co-expression score of a graph.
}
\usage{
getDC(g)
}

\arguments{
  \item{g}{
	An igraph object with two vertex attributes ("name" and "DE") and one edge attribute ("DC").
}
}


\value{
	The differential co-expression score of \code{g}.
}
