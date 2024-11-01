% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichmentAnalysis.R
\name{enrich_dotplot}
\alias{enrich_dotplot}
\title{Create a Dot Plot for Enrichment Results}
\usage{
enrich_dotplot(e.res, topN = 4)
}
\arguments{
\item{e.res}{A list of enrichment results (clusterProfiler::enricher objects).}

\item{topN}{Integer. The number of top terms to include for each group. Default is 4.}
}
\value{
A ggplot object representing the dot plot of enrichment results.
}
\description{
This function generates a dot plot visualization for enrichment analysis results.
}
\details{
The function creates a dot plot where each dot represents an enriched term.
The size of the dot indicates the count of genes, and the color represents the
adjusted p-value (-log10 transformed).
}
\examples{
# Assuming 'enrichment_results' is your list of enrichment results
# plot <- enrich_dotplot(enrichment_results)
# print(plot)

}
