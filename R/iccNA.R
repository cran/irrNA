#' @title Intraclass correlation coefficients (ICCs) -- generalized for randomly incomplete datasets
#' 
#' @description This function computes intraclass correlation coefficients (ICCs) as indices of 
#' interrater reliability or agreement based on cardinally scaled data. This function also works on 
#' (unbalanced) incomplete datasets without any imputation of missing values (\emph{NA}s) or (row- or 
#' cloumn-wise) omissions of data! p-values and confidence intervals are provided. In case of 
#' extreme input data (e.g. zero variances) output \emph{NaN}s are avoided by approximation.
#' 
#' @param ratings n*m matrix or data frame; n objects (rows), m raters (columns)
#' @param rho0 numeric value; correlation in population (\eqn{\rho}) according to the null 
#' hypothesis (0 is default)
#' @param conf numeric value; confidence level (95\% is default)
#' @param detail logical; if TRUE, returns additional information (sums of squares, degrees of 
#' freedom, the means per object, data corrected for the raters' biases)
#' @param oneG logical; if TRUE, the ipsation (correction for the raters' effects) is done the 
#' simple way, using the difference of each raters mean to the one grand mean (\eqn{G}) of all 
#' values to estimate the raters' biases. If FALSE the weighted sub-means (\eqn{G_j} of those
#' objects that an individual rater \eqn{j} rated are used instead (cp. Brueckl, 2011, 
#' Equation 4.30).
#' @param Cs numeric value; denominator (10000 is default) of the effect-size-criterion to stop 
#' iteration of the correction for the raters' biases; the enumerator denotes a small effect 
#' (\eqn{\eta}-squared = 1\%)
#' @return
#'    \item{ICCs}{data frame containing the intraclass correlation coefficients, the corresponding 
#'    p-values, and confidence intervals}
#'    \item{n}{number of rated objects}
#'    \item{k}{maximum number of raters per object}
#'    \item{amk}{mean number of ratings per object}
#'    \item{k_0}{approximate harmonic mean (cp. Ebel, 1951) of the number of ratings per object}
#'    \item{n_iter}{number of iterations for correcting for the raters' biases}
#'    \item{corr_ratings}{ratings, corrected for the individual raters' biases}
#'    \item{amO}{means of ratings for each object, based on (1) the original data and on (2) the 
#'    data that 
#'    are corrected for the raters' biases}
#'    \item{oneway}{statistics for the oneway ICCs}
#'    \item{twoway}{statistics for the twoway ICCs}
#' @details This function is able to compute ICCs on randomly incomplete (i.e. unbalanced) data sets. 
#' Thus, both an imputation of missing values (\emph{NA}s) and row-wise or column-wise omissions of 
#' data are obsolete. Working on complete datasets, it yields the same results as the common 
#' functions, e.g. \link[irrNA]{icc_corr}.\cr
#' The method of Ebel (1951) is used to calculate the oneway ICCs. The solution for the twoway 
#' ICCs is derived from the oneway solution (cp. Brueckl, 2011, p. 96 ff.): The raters' individual 
#' effects (biases) are estimated, reducing this problem again to the oneway problem (cp. Greer & 
#' Dunlap, 1997).\cr
#' This estimation can be done using the difference of a certain (\eqn{j}) rater's mean to the 
#' grand mean (\eqn{G}) or to the sub-mean (\eqn{G_j}) representing only those objects 
#' that were rated by this rater. The first method is fail-safe. The second method is thought to 
#' provide the more precise estimates (of the raters' biases), the more the mean of the true values 
#' of the objects that each rater rated differ from the grand mean, e.g. if there are raters that 
#' only rate objects with low true values (and therefore also other raters that only rate objects 
#' with high true values).\cr
#' If the second method is chosen and if the ratings are unbalanced, which happens most of the time 
#' if not all raters rated all objects, the raters' biases cannot be determined exactly -- but as 
#' approximately as desired. This approximation needs an iteration, thus a stop criterion 
#' (\code{Cs}): 
#' The iteration is stopped, when the difference in the raters' effect size (\eqn{\eta}-squared) 
#' between subsequent iterations would be equal to or smaller than the \code{Cs}th part of a small 
#' effect (i.e. \eqn{\eta}-squared = 1\%).\cr
#' \cr
#' Just as in \link[irrNA]{icc_corr} and \link[irr]{icc}, the designation established by McGraw 
#' & Wong (1996) -- \emph{A} for \emph{absolute agreement} and \emph{C} for \emph{consistency} 
#' -- is used to differ between the (twoway) ICCs that rely on different cases and thus must be 
#' interpreted differently.\cr
#' \cr
#' The generalization of the procedure entails a generalization of the three cases that 
#' differentiate the ICCs (cp. Shrout & Fleiss, 1979):\cr
#'  - Case 1 (oneway case, treated by ICC(1) and ICC(k)): \cr
#'  Each object -- of a sample that was randomly drawn from the population of objects; also 
#'  holds true for case 2 and case 3 -- is rated by (a different number of) different raters 
#'  that were randomly drawn from the population of raters.\cr
#'  - Case 2 (twoway case, treated by ICC(A,1) and ICC(A,k)): \cr
#'  Each object is rated by a random subset of the group of raters that is drawn randomly from 
#'  the population of raters.\cr
#'  - Case 3 (twoway case, treated by ICC(C,1) and ICC(C,k)): \cr
#'  Each object is rated by a random subset of the group of all relevant (i.e. fixed) raters.\cr
#' \cr
#'  Output NaNs, that usually occur (see e.g. \link[irr]{icc} or \link[irrNA]{icc_corr}) in case 
#'  of extreme input data (e.g. in case of zero variance(s), within or between objects) are 
#'  avoided by approximation from little less extreme input data. Warning messages are given in 
#'  these cases.\cr
#'  \cr
#' @references Brueckl, M. (2011). Statistische Verfahren zur Ermittlung der 
#' Urteileruebereinstimmung. in: Altersbedingte Veraenderungen der Stimme und Sprechweise von 
#' Frauen, Berlin: Logos, 88--103.
#' @references Ebel, R.L. (1951). Estimation of the reliability of ratings. Psychometrika, 16(4), 
#' 407--424.
#' @references Greer, T., & Dunlap, W.P. (1997). Analysis of variance with ipsative measures. 
#' Psychological Methods, 2, 200--207.
#' @references McGraw, K.O., & Wong, S.P. (1996). Forming inferences about some intraclass 
#' correlation coefficients. Psychological Methods, 1, 30--46.
#' @references Shrout, P.E., & Fleiss, J.L. (1979). Intraclass correlations: uses in assessing rater
#' reliability. Psychological Bulletin, 86(2), 420--428.
#' @author Markus Brueckl
#' @seealso \code{\link[irrNA]{kendallNA}, \link[irrNA]{icc_corr}, \link[irr]{icc}}
#' @importFrom stats complete.cases pf qf var
#' @export iccNA
#' @examples
#' # Example 1:
#' data(ConsistNA)
#' # ConsistNA exhibits missing values, a perfect consistency, and 
#' # a moderate agreement between raters:
#' ConsistNA
#' # Common ICC-algorithms fail, since each row as well as each 
#' # column of ConsistNA exhibits unfilled cells and these missing 
#' # data are omitted column-wise or row-wise:
#' library(irr)
#' icc(ConsistNA, r0=0.3)
#' # Ebel's (1951) method for computing ICC(1) and ICC(k) that is 
#' # implemented in iccNA can cope with such data without an 
#' # omission or an imputation of missing values, but still can 
#' # not depict the raters' interdependency...
#' iccNA(ConsistNA, rho0=0.3)
#' # ...but generalizations of Ebel's method for the twoway ICCs 
#' # are able to assess moderate agreement (ICC(A,1) and ICC(A,k)) 
#' # and perfect consistency (ICC(C,1) and ICC(C,k)), assuming that 
#' # the data were acquired under case 2 or case 3, see Details in 
#' # the Help file.
#' #
#' # Example 2:
#' data(IndepNA)
#' # IndepNA exhibits missing values and zero variance between 
#' # the raters just as well as between the objects:
#' IndepNA
#' # Again, common ICC-algorithms fail:
#' icc(IndepNA)
#' # But iccNA is able to include all available data in its 
#' # calculation and thereby to show the perfect independence of 
#' # the ratings:
#' iccNA(IndepNA)
#' #
#' # Example 3:
#' # The example provided by Ebel (1951, Tables 2 and 3):
#' # data(Ebel51)
#' Ebel51
#' # iCCNA achieves to include all available ratings and to assess 
#' # twoway ICCs, assuming that the data were acquired under 
#' # case 2 or case 3:
#' iccNA(Ebel51, detail=TRUE)
"iccNA" <- function (ratings, rho0 = 0, conf = 0.95, detail=FALSE , oneG=TRUE, Cs = 10000){
  if (1 <= rho0 || rho0 < 0){
    write ("Argument error: rho0 must be a value of the interval [0,1)!", stderr())
    stop()
  }
  if (1 <= conf || conf <= 0){
    write ("Argument error: conf must be a value between 0 and 1!", stderr())
    stop()
  }
  if (Cs <= 0){
    write ("Argument error: Cs must be greater than 0!", stderr())
    stop()
  }
  
  X <- ratings
  if (is.matrix(X)==TRUE){
    X <- as.data.frame(X)
  }
  
  # output the names of columns/rows without score
  for(x in seq(ncol(X))){
    if (((colSums(!is.na(X)) > 0)[x])==FALSE){
      nom <- names((colSums(!is.na(X)) > 0)[x])
      write (c("Warning: Rater \'",nom,"\' provides no rating!"), stderr(), sep="",ncolumns=3)
      write (c("  -> Rater \'",nom,"\' will be excluded from ICC calculation."), stderr(), sep="",ncolumns=3)
    }
  }
  for(x in seq(nrow(X))){
    if (((rowSums(!is.na(X)) > 0)[x])==FALSE){
      nom <- names((rowSums(!is.na(X)) > 0)[x])
      write (c("Warning: Object \'",nom,"\' receives no rating!"), stderr(), sep="",ncolumns=3)
      write (c("  -> Object \'",nom,"\' will be excluded from ICC calculation."), stderr(), sep="",ncolumns=3)
    }
  }
  # delete columns and rows without any score (raters that did not rate and objects that were not 
  # rated)
  X <- X[,(colSums(!is.na(X)) > 0)]
  X <- X[(rowSums(!is.na(X)) > 0),]
  
  # output the names of columns without variance
  for(x in seq(ncol(X))){
    if ((var(X[,x], na.rm=TRUE))==0 || is.na(var(X[,x], na.rm=TRUE))){
      nom <- names((colSums(!is.na(X)) > 0)[x])
      if (is.na(var(X[,x], na.rm=TRUE))) {
        write (c("Warning: Rater \"",nom,"\" only rates one object!"), stderr(), sep="",ncolumns=3)
      }
      else {
        write (c("Warning: There is no variance in the ratings of rater \"",nom,"\"!"), stderr(), sep="",ncolumns=3) 
      }
    }
  }

  cnam <- colnames(X)
  rnam <- rownames(X)
  
  # transform the (per definition one-tailed) confidence level into the two-tailed value needed in PDFs
  conf <- (1 - conf) / 2 + conf
  
  # calculation of oneway ICCs, cp. Ebel (1951)
  n <- nrow(X)
  k <- ncol(X)
  
  nj <- colSums(!is.na(X))
  ki <- rowSums(!is.na(X))

  olk <- mean(ki)
  k_0 <- 1/(n-1)*(sum(ki)-sum(ki^2)/sum(ki))# cp. Ebel (1951) or Snedecor & Cochran (1989), p. 245
  
  G <- sum(rowSums(X, na.rm = TRUE))
  olG <- G / (olk*n)
  olO <- rowMeans(X, na.rm = TRUE)
  olO <- as.matrix(olO)#t(t(olO)) #todo: fix: double transpose

  olR <- colMeans(X, na.rm = TRUE)
  
  SS_bO <- sum(ki*(olO^2 - olG^2))
  df_bO <- n-1
  MS_bO <- SS_bO / df_bO

  SS_iO <- sum(colSums(X^2, na.rm = TRUE)) - sum(ki*(olO^2))
  df_iO <- n*(olk-1)
  MS_iO <- SS_iO / df_iO
  
  ICC11 <- (MS_bO - MS_iO) / (MS_bO + (k_0-1)*MS_iO)
  ICC1k <- (MS_bO - MS_iO) / MS_bO
  
  F11 <- (MS_bO / MS_iO) * (1-rho0) / (1 + (k_0-1)*rho0)
  F1k <- (MS_bO / MS_iO) * (1-rho0)
  
  dfZ <- n-1
  dfN1 <- n*(olk-1)
  
  p11 <- 1 - pf(F11,dfZ,dfN1)
  p1k <- 1 - pf(F1k,dfZ,dfN1)
  
  # confidence intervals of the oneway ICCs
  Femp <- MS_bO / MS_iO
  Ftablo <- qf(conf,dfZ,dfN1)
  Ftabup <- qf(conf,dfN1,dfZ)

  FL1 <- Femp / Ftablo
  FU1 <- Femp * Ftabup

  ICCL11 <- (FL1 - 1)/(FL1 + k_0 -1)
  ICCU11 <- (FU1 - 1)/(FU1 + k_0 -1)

  ICCL1k <- 1 - 1/FL1
  ICCU1k <- 1 - 1/FU1
  

  # calculation of twoway ICCs
  # 1st ipsation (correction of the rating data for the raters' individual bias), cp. Bortz & 
  # Doering (2006, p. 274 ff.) or Greer & Dunlap (1997)
  SS_tot <- sum(as.matrix(colSums((X-olG)^2, na.rm = TRUE)))
  X_rowMeans <- rowMeans(X, na.rm = TRUE)
  
  Result <- vector(mode="list", k)
  for(i in seq(k)){
    Result[[i]] <- k
    Result[[i]] <- which(!is.na(X[i]))
  }
  
  olGj <- vector()

  for(i in seq(k)){
    f <- 1/sum(!is.na(X[i]))# Flo schrieb hier: 1/colSums(!is.na(X[i]))
    e <- sum(X_rowMeans[Result[[i]]])
    olGj[i] <- f * e
    olGj <- as.vector(olGj)
  }

  # determines, which grand mean(s) is/are to be used in the ipsation
  if (oneG == TRUE){
    olD <- olG - olR}
  else{
    olD <- olGj - olR}
  
  Y <- vector()
  for(i in seq(k)){
    Y[i] <- olD[i] + X[i]
  }

  Y <- do.call(rbind, Y)
  Y <- t(Y)
  
  olOy <- rowMeans(Y, na.rm = TRUE)
  olOy <- as.matrix(olOy)#t(t(olOy))
  
  SS_bR_diff <- sum(nj*olD^2)

  eta2_diff <- SS_bR_diff / SS_tot

  # iteration of the ipsation (cp. Brueckl, 2011, p. 96 ff.) until the difference 
  # between subsequent ipsations would be irrelevant, i.e. until the effect of this 
  # correction is less than 1/Cs of a "small" effect (eta2=0.01, cp. Cohen, 1988)
  icount <- 0
  Gbqy <- 0
#  while ((eta2_diff > 0.01/Cs) && (sum(nj*(olR^2-Gbqy^2)) >= 0)){
  while (eta2_diff > 0.01/Cs){  
    Y_teil <- Y[complete.cases(Y), ]
    Y <- as.data.frame(Y)
    Y_rowMeans <- rowMeans(Y, na.rm = TRUE)
  
    Result <- vector(mode="list", k)
    for(i in seq(k)){
      Result[[i]] <- k
      Result[[i]] <- which(!is.na(Y[i]))
    }  
    
    Gbqy <- vector()
    for(i in seq(k)){    
      f <- 1/sum(!is.na(Y[i]))# Flo schrieb hier: 1/colSums(!is.na(Y[i]))
      e <- sum(Y_rowMeans[Result[[i]]])
      Gbqy[i] <- f * e
      Gbqy <- as.vector(Gbqy)
    }
  
    Bqy <- colMeans(Y, na.rm = TRUE)
    # determines, which grand mean(s) is/are to be used in the ipsation
    if (oneG == TRUE){
      olD <- olG - Bqy[seq(k)]}
    else{
      olD <- Gbqy - Bqy[seq(k)]}

    Z <- vector()
    for(i in seq(k)){
      Z[i] <- olD[i] + Y[i]
    }
    Z <- do.call(rbind, Z)
    Z <- t(Z)
    
    olOzj <- rowMeans(Z, na.rm = TRUE)
    olOzj <- as.matrix(olOzj)#t(t(olOzj))
    SS_bR_diff <- sum(nj*olD^2)

    eta2_diff <- SS_bR_diff / SS_tot
    
    if (eta2_diff > 0.01/Cs) {
      olGj <- Gbqy
      olOy <- olOzj
      Y <- Z
    }

  icount <- icount+1
  }

  # calculation of the consistency ICCs
  SS_bOy <- sum(ki*(olOy^2 - olG^2))
#  SS_bOy <- sum(ki*(olOy-olG)^2)
  df_bOy <- n-1
  MS_bOy <- SS_bOy / df_bOy
  # for reasons of computational precision SS_bOy and thus MS_bOy may wrongly equal to 
  # a value that is (slightly) below zero (cp. near_oVbR_all)
  if (MS_bOy < 0) {
    MS_bOy = 0
  } 
  SS_res <- sum(colSums(Y^2, na.rm = TRUE)) - sum(ki*(olOy^2))
  df_res <- (n-1)*(olk-1)
  MS_res <- SS_res / df_res

  ICC31 <- (MS_bOy - MS_res) / (MS_bOy + (k_0-1)*MS_res)
  ICC3k <- (MS_bOy - MS_res) / MS_bOy

  # p-values of the consistency ICCs
  F31 <- (MS_bOy / MS_res) * (1-rho0) / (1 + (k_0-1)*rho0)
  F3k <- (MS_bOy / MS_res) * (1-rho0)

  dfN3 <- (n-1)*(olk-1)

  p31 <- 1 - pf(F31,dfZ,dfN3)
  p3k <- 1 - pf(F3k,dfZ,dfN3)
  
  # confidence intervals of the consistency ICCs
  FL3 <- (MS_bOy / MS_res) / qf(conf,dfZ,dfN3)
  FU3 <- (MS_bOy / MS_res) * qf(conf,dfN3,dfZ)

  ICCL31 <- (FL3 - 1)/(FL3 + k_0 -1)
  ICCU31 <- (FU3 - 1)/(FU3 + k_0 -1)

  ICCL3k <- 1 - 1/FL3
  ICCU3k <- 1 - 1/FU3

  
  # calculation of the absolute agreement ICCs
  # determines, which great mean(s) is/are to be used in the ipsation
  if (oneG == TRUE){
    SS_bR <- sum(nj*(olR^2 - olG^2))}
  else{
    SS_bR <- sum(nj*(olR^2 - olGj^2))}

  df_bR <- k-1
  MS_bR <- SS_bR / df_bR

  ICC21 <- (MS_bOy - MS_res) / (MS_bOy + (k_0-1)*MS_res + k_0/n*(MS_bR - MS_res))
  # avoid double negations (on both sides of the ratio) and therewith (high, >1) positive ICCs, 
  # that would occur on very unreliable data, cp. kack, ICC(A,k)
  nenn <- MS_bOy + (MS_bR - MS_res)/n
  if (nenn < 0) { # in this case the ICC2k enumerator is always < 0 (MS_res > MS_bOy)!!!
     nenn = 0     # therewith zero is the most extreme meaningful value
  }
  ICC2k <- (MS_bOy - MS_res) / nenn

  # p-values for the absolute agreement ICCs according to McGraw & Wong (1996):
  a <- k_0*rho0 / (n*(1-rho0))
  b <- 1 + (k_0*rho0*(n-1)) / (n*(1-rho0))
  c <- rho0 / (n*(1-rho0))
  d <- 1 + (rho0*(n-1))/(n*(1-rho0))

  F21 <- MS_bOy / (a*MS_bR + b*MS_res)
  F2k <- MS_bOy / (c*MS_bR + d*MS_res)
  
  dfN21 <- (a*MS_bR + b*MS_res)^2 / ((a*MS_bR)^2/df_bR + (b*MS_res)^2/df_res)
  dfN2k <- (c*MS_bR + d*MS_res)^2 / ((c*MS_bR)^2/df_bR + (d*MS_res)^2/df_res)

  p21 <- 1 - pf(F21,dfZ,dfN21)
  p2k <- 1 - pf(F2k,dfZ,dfN2k)

  # confidence intervals of the absolute agreement ICCs, cp. McGraw & Wong (1996) 
  aKI1 <- k_0*ICC21 / (n*(1-ICC21))
  bKI1 <- 1 + (k_0*ICC21*(n-1)) / (n*(1-ICC21))
  dfN2KI1 <- (aKI1*MS_bR + bKI1*MS_res)^2 / ((aKI1*MS_bR)^2 /df_bR + (bKI1*MS_res)^2/df_res)
  # correction for the case of dfN2KI1 approximating 0; ...qf can't handle!
  # AND degrees of freedom (df) smaller than 1 do not make much sense anyway, so, 
  # assuming 0.1 to be the smallest df value that qf can handle properly
  if ((dfN2KI1 < 0.1) || (is.nan(dfN2KI1)==1)){
    write ("Warning: df for estimating the CI of ICC(A,1) reaches value far below one!", stderr())
    write ("  -> The confidence limits of ICC(A,1) are derived by approximation.", stderr())
    dfN2KI1 <- 0.1
  }
  FstarL1 <- qf(conf,dfZ,dfN2KI1)
  FstarU1 <- qf(conf,dfN2KI1,dfZ)

  ICCL21 <- n*(MS_bOy - FstarL1*MS_res) / (FstarL1*(k_0*MS_bR + (k_0*n - k_0 - n)*MS_res) + 
                                             n*MS_bOy)
  # further correction for the case of dfN2KI1 approximating 0, then FstarL1 may be Inf
  # dispensable??? because of the above
  if (FstarL1==Inf){
    write ("Warning: The F-Ratio for determining ICC(A,1) is infinite!", stderr())
    write ("  -> The lower confidence limit of ICC(A,1) is derived by approximation.", stderr())
    ICCL21 <- ICC21
  }
  ICCU21 <- n*(FstarU1*MS_bOy - MS_res) / (k_0*MS_bR + (k_0*n - k_0 - n)*MS_res + 
                                             n*FstarU1*MS_bOy)

  aKIk <- k_0*ICC2k / (n*(1-ICC2k))
  bKIk <- 1 + (k_0*ICC2k*(n-1)) / (n*(1-ICC2k))
  dfN2KIk <- (aKIk*MS_bR + bKIk*MS_res)^2 / ((aKIk*MS_bR)^2 /df_bR + (bKIk*MS_res)^2/df_res)

  FstarLk <- qf(conf,dfZ,dfN2KIk)
  FstarUk <- qf(conf,dfN2KIk,dfZ)

  ICCL2k <- n*(MS_bOy - FstarLk*MS_res) / (FstarLk*(MS_bR - MS_res) + n*MS_bOy)
  ICCU2k <- n*(FstarUk*MS_bOy - MS_res) / (MS_bR - MS_res + n*FstarUk*MS_bOy)
  
  # if (MS_iO == 0 && MS_bO >= 0) then 
  # Femp = Inf, thus ICCL11 = ICCU11 = Inf/Inf, but approximating 1 when Femp takes higher values
  # and dfN2KIk = NaN (since MS_bR and MS_res both equal zero) but approximating (n-1)*(k-1)
  # cp. oVbO_all
  if (MS_iO == 0 && MS_bO >= 0) {
    write ("Warning: The variance within objects (residual and between raters) equals zero!", stderr())
    write ("  -> The confidence limits of ICC(1) and ICC(A,k) are derived by approximation.", stderr())
    ICCL11 <- 1
    ICCU11 <- 1
    ICCL2k <- 1
    ICCU2k <- 1
    write ("  -> The p-values for ICC(A,1) and ICC(A,k) are derived by approximation.", stderr())
    p21 <- 0
    p2k <- 0
  }
  
  # correction for the case of MS_bOy + (MS_bR - MS_res)/n <= 0, i.e. there 
  # is nearly no variance between objects nor between raters, but a lot residual variance
  # cp. kack, ICC(A,k)
  if (ICC2k == -Inf) {
    write ("Warning: The ICC(A,k) is negatively infinite!", stderr())
    write ("  -> The confidence limits of ICC(A,k) are derived by approximation.", stderr())
    ICCL2k <- -Inf
    ICCU2k <- -Inf
  }

  # correction for the case of MS_res = 0, i.e. there is no unexplained 
  # variance, all variance is either between objects or between raters
  if (MS_bOy > 0 && MS_res == 0) {
    write ("Warning: The residual variance equals zero!", stderr())
    write ("  -> The confidence limits of ICC(C,1) and ICC(C,k) are derived by approximation.", stderr())
    ICCL31 <- 1
    ICCU31 <- 1
  }
  
  # since in the case of (MS_res=0 AND rho0=0) (cp. ConsistNA) F21 and F2k are infinite
  # except the case of MS_bOy=0, too (see below, cp. oVbR_all), then F21 and F2k are 0
  if (MS_res == 0 && rho0 == 0 && MS_bOy > 0) {
    write ("Warning: The residual variance and rho0 equal zero!", stderr())
    write ("  -> The p-values for ICC(A,1) and ICC(A,k) are derived by approximation.", stderr())
    p21 <- 0
    p2k <- 0
  }

  # correction for the case of there is only variance between raters, cp. oVbR_all
  if (MS_bOy == 0 && MS_res == 0) {
    write ("Warning: The variance between objects and the residual variance equal zero!", stderr())
    write ("  -> The confidence limits of ICC(A,k) are derived by approximation.", stderr())
    ICCL2k <- 0
    ICCU2k <- 0
    write ("  -> The p-values for ICC(A,1) and ICC(A,k) are derived by approximation.", stderr())
    p21 <- 1
    p2k <- 1
    write ("  -> ICC(C,1) and ICC(C,k) are derived by approximation.", stderr())
    ICC31 <- -1/(k_0-1)
    ICC3k <- -Inf
    write ("  -> The confidence limits of ICC(C,1) and ICC(C,k) are derived by approximation.", stderr())
    ICCL31 <- ICC31
    ICCU31 <- ICC31
    ICCL3k <- ICC3k
    ICCU3k <- ICC3k
    write ("  -> The p-values of ICC(C,1) and ICC(C,k) are derived by approximation.", stderr())
    p31 <- 1
    p3k <- 1
  }

  # all results at a glance
  einseins <- cbind(ICC11, p11, ICCL11, ICCU11)
  einsk <- cbind(ICC1k, p1k, ICCL1k, ICCU1k)
  zweieins <- cbind(ICC21, p21, ICCL21, ICCU21)
  zweik <- cbind(ICC2k, p2k, ICCL2k, ICCU2k)
  dreieins <- cbind(ICC31, p31, ICCL31, ICCU31)
  dreik <- cbind(ICC3k, p3k, ICCL3k, ICCU3k)
  ICC <- rbind(einseins, einsk, zweieins, zweik, dreieins, dreik)
  rownames(ICC) <- c("ICC(1)", "ICC(k)", "ICC(A,1)", "ICC(A,k)", "ICC(C,1)", "ICC(C,k)")
  colnames(ICC) <- c("ICC", "p-value", "lower CI limit", "upper CI limit")

  amO <- cbind(olO, olOy)
  colnames(amO) <- c("original", "corrected")
  
  data_ips <- Y[,1:k]
  n_ips <- icount

  if (detail == FALSE) {
    answ <- list(ICCs = ICC, n_iter = n_ips-1, amk = olk, k_0 = k_0)
  }
  else if (detail == TRUE) {
    bO <- cbind(SS_bO, df_bO)
    iO <- cbind(SS_iO, df_iO)
    df_tot <- df_bO + df_iO
    tot <- cbind(SS_tot, df_tot)
    oneway <- rbind(bO, iO, tot)
    colnames(oneway) <- c("SS", "df")
    rownames(oneway) <- c("between objects", "within objects", "total")
    
    bOy <- cbind(SS_bOy, df_bOy)
    bR <- cbind(SS_bR, df_bR)
    res <- cbind(SS_res, df_res)
    SS_iOy <- SS_bR + SS_res
    df_iOy <- df_bR + df_res
    iOy <- cbind(SS_iOy, df_iOy)
    df_toty <- df_bOy + df_iOy
    toty <- cbind(SS_tot, df_toty)
    twoway <- rbind(bOy, iOy, bR, res, toty)
    colnames(twoway) <- c("SS", "df")
    rownames(twoway) <- c("between objects", "within objects","   between raters", "   residual", 
                          "total")
    
    colnames(data_ips) <- cnam
    rownames(data_ips) <- rnam
    
    answ <- list(ICCs = ICC, n = n, k = k, amk = olk, k_0 = k_0, n_iter = n_ips-1,
                 amO = amO, corr_data = data_ips, oneway_Stats = oneway, twoway_Stats = twoway)
  }
  return(answ)
}