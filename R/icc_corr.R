#' @title Intraclass correlation coefficients (ICCs) for oneway and twoway models -- corrected version 
#' of icc\{irr\}
#' @description Computes single score or average score ICCs as an index of interrater reliability of 
#' quantitative data. Additionally, F-test and confidence interval are computed. icc_corr\{irrNA\}
#' corrects 3 errors of Matthias Gamer's function \link[irr]{icc} (version 0.84.1).
#' 
#' @param ratings n*m matrix or dataframe, n subjects m raters.
#' @param model	a character string specifying if a "oneway" model (default) with row effects random, 
#' or a "twoway" model with column and row effects random should be applied. You can specify just 
#' the initial letter.
#' @param type a character string specifying if "consistency" (default) or "agreement" between 
#' raters should be estimated. If a '"oneway"' model is used, only "consistency" could be computed. 
#' You can specify just the initial letter.
#' @param unit a character string specifying the unit of analysis: Must be one of "single" (default) 
#' or "average". You can specify just the initial letter.
#' @param r0 specification of the null hypothesis r \eqn{\le} r0. Note that a one sided test 
#' (H1: r > r0) is performed.
#' @param conf.level confidence level of the interval.
#' 
#' @return 
#' A list with class '"icclist"' containing the following components:
#' \item{$subjects}{the number of subjects examined.}
#' \item{$raters}{the number of raters.}
#' \item{$model}{a character string describing the selected model for the analysis.}
#' \item{$type}{a character string describing the selected type of interrater reliability.}
#' \item{$unit}{a character string describing the unit of analysis.}
#' \item{$icc.name}{a character string specifying the name of ICC according to McGraw & Wong (1996).}
#' \item{$value}{the intraclass correlation coefficient.}
#' \item{$r0}{the specified null hypothesis.}
#' \item{$Fvalue}{the value of the F-statistic.}
#' \item{$df1}{the numerator degrees of freedom.}
#' \item{$df2}{the denominator degrees of freedom.}
#' \item{$p.value}{the p-value for a two-sided test.}
#' \item{$conf.level}{the confidence level for the interval.}
#' \item{$lbound}{the lower bound of the confidence interval.}
#' \item{$ubound}{the upper bound of the confidence interval.}
#'   
#' @details By this ICC-function three bugs are corrected that were found in the function 
#' \link[irr]{icc} of the irr package (version 0.84.1): \cr
#' Due to the first bug the p-values of ICC(A,1) and ICC(A,k) are computed wrongly: 
#' McGraw & Wong (1996) use the variable "v" both for the computation of the CIs and for the 
#' computation of the p-values. But "v" takes different values in these calculations. In the 
#' implementation of icc\{irr\} (version 0.84.1) this fact is missed. \cr
#' The second correction only affects the rare case of the residual mean square (of the twoway 
#' model) being zero, i.e. the case that the variance in the data may be explained completely 
#' by the two factors (Raters and Objects). In this case the F-value for determining all four
#' twoway p-values is not correctly computet by \link[irr]{icc}.\cr
#' The third correction addresses the problems arising in the rare cases of (a) no part or (b) 
#' nearly no part of variance may be explained by both factors.
#' 
#' @references McGraw, K.O., & Wong, S.P. (1996). Forming inferences about some intraclass 
#' correlation coefficients. Psychological Methods, 1, 30--46.
#' @references Shrout, P.E., & Fleiss, J.L. (1979), Intraclass correlation: 
#' uses in assessing rater reliability. Psychological Bulletin, 86, 420--428.
#'
#' @author Matthias Gamer, Markus Brueckl
#' @seealso \code{\link[irr]{icc}}, \code{\link[irrNA]{iccNA}}
#' @importFrom stats pf qf var
#' @export icc_corr
#' @examples
#' # Example 1:
#' data(EbelFILL)
#' # EbelFILL is a rather arbitrary data set:
#' EbelFILL
#' # If twoway agreement ICCs are computed (e.g. the single 
#' # measure) with icc{irr}, the 2nd df of F and thus the 
#' # p-value is erroneous (please install and load the irr 
#' # package):
#' #icc(EbelFILL, model="twoway", type="agreement")
#' # icc_corr calculates correctly: 
#' icc_corr(EbelFILL, model="twoway", type="agreement")
#' # 
#' # Example 2:
#' data(Consist)
#' # Consist exhibits a perfect consistency and 
#' # a moderate absolute agreement between raters:
#' Consist
#' # If twoway ICCs are computed with icc{irr}, the F-value is smaller
#' # than zero (!) and thus the p-value is enourmously erroneous:
#' #icc(Consist, model="twoway", type="consistency", unit="average")
#' # icc_corr calculates correctly: 
#' icc_corr(Consist, model="twoway", type="consistency", unit="average")
#' #
#' # Example 3:
#' data(Indep)
#' # Indep exhibits zero variance between the raters just as 
#' # well as between the objects:
#' Indep
#' # Errors occur, if twoway agreement ICCs are computed with icc{irr}:
#' # ICC(A,k) just as well as its CI-bounds are (falsely) positive 
#' # and greater than 1...
#' #icc(Indep, model="twoway", type="agreement", unit="average")
#' # ...but must be -Inf, just as icc_corr shows:
#' icc_corr(Indep, model="twoway", type="agreement", unit="average")
#' # ICC(A,1): 2nd df of F and thus the p-value are NaN
#' #icc(Indep, model="twoway", type="agreement")
#' # icc_corr calculates correlctly:
#' icc_corr(Indep, model="twoway", type="agreement")
"icc_corr"<-function (ratings, model = c("oneway", "twoway"), type = c("consistency", 
                                                           "agreement"), unit = c("single", "average"), r0 = 0, conf.level = 0.95) 
{
  ratings <- as.matrix(ratings)#as.matrix(na.omit(ratings))
  model <- match.arg(model)
  type <- match.arg(type)
  unit <- match.arg(unit)
  alpha <- 1 - conf.level
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  SStotal <- var(as.numeric(ratings)) * (ns * nr - 1)
  MSr <- var(apply(ratings, 1, mean)) * nr
  MSw <- sum(apply(ratings, 1, var)/ns)
  MSc <- var(apply(ratings, 2, mean)) * ns
  MSe <- (SStotal - MSr * (ns - 1) - MSc * (nr - 1))/((ns - 
                                                         1) * (nr - 1))
  # for the case of missing data
  if (is.na(SStotal)==1){
    warning ('Data is missing in your data set. Use iccNA instead!')
  }
  else {
  # since -- computed this way -- MSe may be smaller than zero due to rounding errors
  # cp. Consist
  if (MSe < 0){
    MSe = 0
  }

  if (unit == "single") {
    if (model == "oneway") {
      icc.name <- "ICC(1)"
      coeff <- (MSr - MSw)/(MSr + (nr - 1) * MSw)
      Fvalue <- MSr/MSw * ((1 - r0)/(1 + (nr - 1) * r0))
      df1 <- ns - 1
      df2 <- ns * (nr - 1)
      p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
      FL <- (MSr/MSw)/qf(1 - alpha/2, ns - 1, ns * (nr - 
                                                      1))
      FU <- (MSr/MSw) * qf(1 - alpha/2, ns * (nr - 1), 
                           ns - 1)
      lbound <- (FL - 1)/(FL + (nr - 1))
      ubound <- (FU - 1)/(FU + (nr - 1))
    }
    else if (model == "twoway") {
      if (type == "consistency") {
        icc.name <- "ICC(C,1)"
        coeff <- (MSr - MSe)/(MSr + (nr - 1) * MSe)
        Fvalue <- MSr/MSe * ((1 - r0)/(1 + (nr - 1) * 
                                         r0))
        df1 <- ns - 1
        df2 <- (ns - 1) * (nr - 1)
        p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
        FL <- (MSr/MSe)/qf(1 - alpha/2, ns - 1, (ns - 
                                                   1) * (nr - 1))
        FU <- (MSr/MSe) * qf(1 - alpha/2, (ns - 1) * 
                               (nr - 1), ns - 1)
        lbound <- (FL - 1)/(FL + (nr - 1))
        ubound <- (FU - 1)/(FU + (nr - 1))
       }
      else if (type == "agreement") {
        icc.name <- "ICC(A,1)"
        coeff <- (MSr - MSe)/(MSr + (nr - 1) * MSe + 
                                (nr/ns) * (MSc - MSe))
        a <- (nr * r0)/(ns * (1 - r0))
        b <- 1 + (nr * r0 * (ns - 1))/(ns * (1 - r0))
        Fvalue <- MSr/(a * MSc + b * MSe)
     
        # v must be calculated here (in addition to below), since a and b are set here for the
        # computation of the p.value, cp. McGraw & Wong (1996).
        # df1, df2, and p.value must be computed before a, b, and thus v take the new values
        # that are needed for the determination of the confidence limits.
        v <- (a * MSc + b * MSe)^2/((a * MSc)^2/(nr - 
                                                   1) + (b * MSe)^2/((ns - 1) * (nr - 1)))
        df1 <- ns - 1
        df2 <- v
        p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)

        a <- (nr * coeff)/(ns * (1 - coeff))
        b <- 1 + (nr * coeff * (ns - 1))/(ns * (1 - coeff))
        v <- (a * MSc + b * MSe)^2/((a * MSc)^2/(nr - 
                                                   1) + (b * MSe)^2/((ns - 1) * (nr - 1)))
        FL <- qf(1 - alpha/2, ns - 1, v)
        FU <- qf(1 - alpha/2, v, ns - 1)
        lbound <- (ns * (MSr - FL * MSe))/(FL * (nr * 
                                                   MSc + (nr * ns - nr - ns) * MSe) + ns * MSr)
        ubound <- (ns * (FU * MSr - MSe))/(nr * MSc + 
                                             (nr * ns - nr - ns) * MSe + ns * FU * MSr)
      }
    }
  }
  else if (unit == "average") {
    if (model == "oneway") {
      icc.name <- paste("ICC(", nr, ")", sep = "")
      coeff <- (MSr - MSw)/MSr
      Fvalue <- MSr/MSw * (1 - r0)
      df1 <- ns - 1
      df2 <- ns * (nr - 1)
      p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
      FL <- (MSr/MSw)/qf(1 - alpha/2, ns - 1, ns * (nr - 
                                                      1))
      FU <- (MSr/MSw) * qf(1 - alpha/2, ns * (nr - 1), 
                           ns - 1)
      lbound <- 1 - 1/FL
      ubound <- 1 - 1/FU
    }
    else if (model == "twoway") {
      if (type == "consistency") {
        icc.name <- paste("ICC(C,", nr, ")", sep = "")
        coeff <- (MSr - MSe)/MSr
        Fvalue <- MSr/MSe * (1 - r0)
        df1 <- ns - 1
        df2 <- (ns - 1) * (nr - 1)
        p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)
        FL <- (MSr/MSe)/qf(1 - alpha/2, ns - 1, (ns - 
                                                   1) * (nr - 1))
        FU <- (MSr/MSe) * qf(1 - alpha/2, (ns - 1) * 
                               (nr - 1), ns - 1)
        lbound <- 1 - 1/FL
        ubound <- 1 - 1/FU
      }
      else if (type == "agreement") {
        icc.name <- paste("ICC(A,", nr, ")", sep = "")
        # avoid double negations (on both sides of the ratio) and therewith (high, >1) positive 
        # ICCs, that would occur on very unreliable data, cp. Indep
        denom <- MSr + (MSc - MSe)/ns
        if (denom < 0) { # in this case the ICC2k enumerator always is less than 0!!!
          denom = 0      # therewith zero is the most extreme meaningful value
        }
        coeff <- (MSr - MSe)/denom
        a <- r0/(ns * (1 - r0))
        b <- 1 + (r0 * (ns - 1))/(ns * (1 - r0))
        Fvalue <- MSr/(a * MSc + b * MSe)
        # v must be calculated here (in addition to below), since a and b are set here for the
        # computation of the p.value, cp. McGraw & Wong (1996).
        # df1, df2, and p.value must be computed before a, b, and thus v take the new values
        # that are needed for the determination of the confidence limits.
        v <- (a * MSc + b * MSe)^2/((a * MSc)^2/(nr - 
                                                   1) + (b * MSe)^2/((ns - 1) * (nr - 1)))
        df1 <- ns - 1
        df2 <- v
        p.value <- pf(Fvalue, df1, df2, lower.tail = FALSE)

        a <- (nr * coeff)/(ns * (1 - coeff))
        b <- 1 + (nr * coeff * (ns - 1))/(ns * (1 - coeff))
        v <- (a * MSc + b * MSe)^2/((a * MSc)^2/(nr - 
                                                   1) + (b * MSe)^2/((ns - 1) * (nr - 1)))
        FL <- qf(1 - alpha/2, ns - 1, v)
        FU <- qf(1 - alpha/2, v, ns - 1)
        lbound <- (ns * (MSr - FL * MSe))/(FL * (MSc - 
                                                   MSe) + ns * MSr)
        ubound <- (ns * (FU * MSr - MSe))/(MSc - MSe + 
                                             ns * FU * MSr)
        # correction for the case of MSr + (MSc - MSe)/ns <= 0, i.e. there 
        # is nearly no variance between objects or between raters, but a lot residual variance,
        # cp. Indep
        if (coeff == -Inf) {
          lbound <- -Inf
          ubound <- -Inf
        }
      }
    }
  }
  rval <- structure(list(subjects = ns, raters = nr, model = model, 
                         type = type, unit = unit, icc.name = icc.name, value = coeff, 
                         r0 = r0, Fvalue = Fvalue, df1 = df1, df2 = df2, p.value = p.value, 
                         conf.level = conf.level, lbound = lbound, ubound = ubound), 
                    class = "icclist")
  return(rval)
  }
}