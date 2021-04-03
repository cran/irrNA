#' @title Kendall's coefficient of concordance W -- generalized for randomly incomplete datasets
#' @description This function computes Kendall's coefficient of concordance W that is an index of 
#' interrater reliability for ordinal ratings. This function also works on incomplete datasets without 
#' any imputation of missing values or (row- or cloumn-wise) omissions of data.
#' 
#' @param X n*m matrix or dataframe; n objects (rows), k  raters (columns)
#' @return 
#'   \item{amrho}{mean Spearman's \eqn{\rho}}
#'   \item{amk}{mean number of (pairwise) ratings per object}
#'   \item{W}{Kendall's coefficient of concordance among raters}
#'   \item{chisqu}{value of the  \eqn{\chi}-squared test statistic}
#'   \item{df}{degrees of freedom}
#'   \item{p}{one-tailed type I error probability (statistical significance)}
#' @details This function is able to calculate W, also on randomly incomplete (i.e. unbalanced) 
#' data sets. Therefor it uses the mean Spearman's \eqn{\rho} of all pairwise comparisons, see Kendall
#' (1962):
#' \deqn{W = [1 + mean \rho_S * (k-1)] / k}
#' where k is the mean number of (pairwise) ratings per object and \eqn{mean \rho_S} is calculated 
#' weighted, according to Taylor (1987), since the pairwise \eqn{\rho_S} are possibly based on a 
#' different number of ratings, what must be reflected in weights.\cr
#' Thus, an imputation of missing values or (row- or cloumn-wise) omissions of data are obsolete. In 
#' case of complete datasets, it yields the same results as usual implementations of Kendall's W, 
#' except for tied ranks. In case of tied ranks, the (pairwise) correction of \eqn{\rho_S} is used, 
#' which (already with complete datasets) results in slightly different values than the tie correction 
#' explicitly specified for W.\cr
#' More details are given in Brueckl (2011).
#' @references Brueckl, M. (2011). Statistische Verfahren zur Ermittlung der Urteileruebereinstimmung. 
#' in: Altersbedingte Veraenderungen der Stimme und Sprechweise von Frauen, Berlin: Logos, 88--103.
#' @references Kendall, M.G. (1962). Rank correlation methods (3rd ed.). London: Griffin.
#' @references Taylor, J.M.G. (1987). Kendall's and Spearman's correlation coefficients in the 
#' presence of a blocking variable. Biometrics, 43, 409--416.
#'
#' @author Markus Brueckl
#' @seealso \code{\link[irrNA]{iccNA}, \link[irr]{kendall}}
#' @export kendallNA
#' @examples
#' # Example 1:
#' data(ConsistNA)
#' # ConsistNA exhibits missing values and a perfect concordance
#' # between raters:
#' ConsistNA
#' # Common W-algorithms fail, since each row as well as each 
#' # column of ConsistNA exhibits unfilled cells and these missing 
#' # data are omitted column-wise or row-wise:
#' library(irr)
#' # try here: kendall(ConsistNA)
#' # But the generalization of Kendall's W implemeted in irrNA 
#' # is able to assess the perfect concordance, assuming that 
#' # the data were at least ordinally scaled and not tied, e.g. 
#' # that each rater really ranked the objects that he rated 
#' # without giving equal ranks to two or more objects.
#' kendallNA(ConsistNA)
#' #
#' # Example 2:
#' data(IndepNA)
#' # IndepNA exhibits missing values and zero variance between 
#' # the raters (just as well as between the objects):
#' IndepNA
#' # Common W-algorithms fail:
#' # try here: kendall(IndepNA)
#' # kendallNA includes all (rater-pairwise) available data in 
#' # its calculation (e.g. only Objects 1--4 when Rater1 and 
#' # Rater2 are correlated):
#' kendallNA(IndepNA)
#' #
#' # Example 3:
#' data(IndepW)
#' # IndepW exhibits missing values and a mean Spearman's rho,
#' # that equals zero:
#' IndepW
#' # Again, common W-algorithms fail:
#' # try here: kendall(IndepW)
#' # kendallNA includes all (rater-pairwise) available 
#' # data:
#' kendallNA(IndepW)
"kendallNA" <- function(X){
  
  N <- nrow(X) # N is number of rows
  m <- ncol(X) # m is number of columns

# create a working Matrix: delete columns with none or only one score (raters providing no 
# information)
  nonNAs <- function(x) {
    as.vector(apply(x, 2, function(x) length(which(!is.na(x)))))
  }  
  X <- X[nonNAs(X) > 1]
  m <- ncol(X) # rewrite m with number of columns
  
# calculate Spearmans rhos of matrix X
# disabling warnings for "standard deviation is zero"
  oldw <- getOption("warn")
  options(warn = -1)
    rho <- stats::cor(X, method="spearman", use="pairwise.complete")
  options(warn = oldw)

# create a matrix "rho2" with same no. of columns and rows and if element is NA set it to 0
  j <- 1:m #replace with seq_along(rho)
  i <- 1:m
  rho2 <- ifelse (is.na(rho[j,i]) == TRUE, 0, rho[j,i])

# upper triangle of matrix is set to NA
  rho3 <- ifelse(upper.tri(rho2, diag = TRUE) == TRUE, NA, rho2) 

# create matrix of weights 
  f <- m-2
  X_teil <- X[-(1:f)]
  n <- t(!is.na(X)) %*% (!is.na(X))
  n <- ifelse(upper.tri(n, diag = TRUE) == TRUE, 0, n)
  w <- n-1 # cp. Taylor (1987)

# replace "-1" (upper half) values by NA
  w <- ifelse(upper.tri(w, diag = TRUE) == TRUE, NA, w)
  
# correct eventual -1 values in (the lower half of) w by 0
  w[w<0] <- 0
  
# weighted average of correlation coefficients 
  wsum <- sum(colSums(w, na.rm = TRUE))
  kq <- mean(nonNAs(t(X)), na.rm = TRUE)
  amrho <- sum(rho3 * w, na.rm = TRUE) / wsum
  
# Kendalls coefficient of concordance 
  W <- (1+amrho*(kq-1))/kq # cp. Kendall (1962, p.95)
  chisqu <- kq*(N-1)*W
  df <- N-1
  p <- 1- stats::pchisq(chisqu, df, lower.tail = TRUE, log.p = FALSE)

  answ <- list("amrho" = amrho, "amk" = kq, "Kendall's W" = W, "chisqu" = chisqu,"df" = df, "p" = p)
return(answ)

}
