\name{getDE}
\alias{getDE}

\title{
	Get differential expression score.
}
\description{
	Get differential expression score of a graph.
}
\usage{
getDE(g)
}

\arguments{
  \item{g}{
	An igraph object with two vertex attributes ("name" and "DE") and one edge attribute ("DC").
}
}

\value{
	The differential expression score of \code{g}.
}
