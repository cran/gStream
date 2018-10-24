\name{gstream}
\alias{gstream}
\title{Sequential Change-Point Detection based on k-Nearest Neighbors}
\description{This function finds change-points in the sequence when the underlying distribution changes. It reports four graph-based test statistics and the analytical approximations for thresholds used in their corresponding stopping rules.
}
\usage{
gstream(distM, L, N0, k, statistics = c("all", "o", "w", "g", "m"),
n0 = 0.3*L, n1 = 0.7*L, ARL = 10000, alpha = 0.05, skew.corr = TRUE)
}
\arguments{
  \item{distM}{A distance matrix constructed based on some distance measure.}
  \item{L}{The number of observations the k-NN graph will be constructed from.}
  \item{N0}{The number of historical observations.}
  \item{k}{A fixed integer used to construct k-NN graph.}
  \item{statistics}{The scan statistic to be computed. A character indicating the type of of scan statistic desired. The default is \code{"all"}.

   \code{"all"}: specifies to compute \bold{a}ll of the scan statistics: original, weighted, generalized, and max-type;

   \code{"o", "ori"} or \code{"original"}: specifies the \bold{o}riginal edge-count scan statistic;

   \code{"w"} or \code{"weighted"}: specifies the \bold{w}eighted edge-count scan statistic;

   \code{"g"} or \code{"generalized"}: specifies the \bold{g}eneralized edge-count scan statistic; and

   \code{"m"} or \code{"max"}: specifies the \bold{m}ax-type edge-count scan statistic.
   }
  \item{n0}{The starting index to be considered as a candidate for the change-point. We recommend you set this to be 0.2*L}
  \item{n1}{The ending index to be considered as a candidate for the change-point. For example, n1 = L-n0.}
  \item{ARL}{The average run length: the expectation of the stopping rule when there is no change-point.}
  \item{alpha}{The probability of an early stop.}
  \item{skew.corr}{If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
}

\value{
  Returns a list with items \code{scanZ}, \code{b} and \code{tauhat} for each type of statistic specified. See below for more details.
  \item{scanZ}{A vector of the test statistic (maximum of the scan statistics) for each time n = N0+1,..N.

  \code{ori}: A vector of the original scan statistics (standardized counts) if statistic specified is "all" or "o".

  \code{weighted}: A vector of the weighted scan statistics (standardized counts) if statistic specified is "all" or "w".

  \code{generalized}: A vector of the generalized scan statistics (standardized counts) if statistic specified is "all" or "g".

  \code{max.type}: A vector of the max-type scan statistics (standardized counts) if statistic specified is "all" or "m".}

  \item{b}{Thresholds used in the stopping rules for each test statistic. These thresholds are based on analytical approximations of the average run length.}
  \item{tauhat}{Estimate of the locations of change-points based on the thresholds.}
}

\seealso{
  \code{\link{gStream}}
}
\examples{
# This example contains two distance matrices (distM1 and distM2).
# Information on how distM1 and distM2 are generated can be found in gStream.

data(Example)

# Example:
 # distM1 is a distance matrix constructed from a dataset with n=40 observation.
  # The first 20 observations are treated as historical observations.
  # It has been determined that there are no change-points among the
  # first 20 observations (see package gSeg for offline change-point detection).
  # There is change in mean when tau = 20 (This means a change happens 20 observations
  # after we start the tests. We start the test at N0+1 = 21.)
  N0 = 20
  L = 20 # the k-nn graph is constructed on only the L most recent observations.
  k = 1

  r1= gstream(distM1, L, N0, k, statistics="all", n0=0.3*L, n1=L-0.3*L,
  ARL=200,alpha=0.05, skew.corr=FALSE)
  # output results based on all four statistics; the scan statistics can be found in r1$scanZ
  r1$tauhat # reports the locations where a change-point is detected
  r1$b # reports the analytical approximations of the thresholds used in the stopping rules
  
  # Set ARL = 10,000
  # r1= gstream(distM1, L, N0, k, statistics="all", n0=0.3*L, n1=L-0.3*L,
  # ARL=10000,alpha=0.05, skew.corr=TRUE) # uncomment to run this function
}
