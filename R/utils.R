#............................................................
# SIR NE S3 Class Overloading
#...........................................................
#' @description function for determining if object is of class GillespieSIRne
#' @noRd
#' @export
is.GillespieSIRne <- function(x) {
  inherits(x, "is.GillespieSIRne")
}


#' @description overload print() function to print summary only
#' @noRd
#' @export
print.GillespieSIRne <- function(x, ...) {

  # print summary only
  cat(crayon::red("Final Epidemic Size:"), summary(x)$FinalEpidemicSize, "\n")
  cat(crayon::blue("Run Time:"), x$runTime, attr(x$runTime, "units"), "\n")

}

#' @description overload summary() function.
#' @noRd
#' @export
summary.GillespieSIRne <- function(x, ...) {
  # send summary only
  return(tidyout.GillespieSIRne(x))
}


#............................................................
# organzing SIR NE output
#...........................................................
#' @description Method assignment
#' @noRd
#' @export
tidyout <- function(x) {
  UseMethod("tidyout")
}

#' @description Function for taking output of SIR NE and lifting it over
#' @noRd
#' @export

tidyout.GillespieSIRne <- function(x) {
  # final epidemic size
  finalsize <- sum(x$I_traj[1,]) + sum(x$Event_traj == "transmission")
  # summary table
  cumeventdf <- tibble::tibble(
    Time = x$Time_traj,
    Event = x$Event_traj,
    Susc = apply(x$S_traj, 1, sum),
    Infxn = apply(x$I_traj, 1, sum),
    Recov = apply(x$R_traj, 1, sum)
  )
  # out
  ret <- list(FinalEpidemicSize = finalsize,
              CumEvents = cumeventdf)
  return(ret)

}


#' @title Bind SIR trajectories
#' @param x list of vectors matrix
#' @noMd
#' @export
tidy_traj_out <- function(x) {
  if (length(x) == 1) {
    ret <- matrix(data = unlist(x), ncol = 1)
    rownames(ret) <- NULL
    colnames(ret) <- "1"
  } else {
    ret <- t(do.call("cbind", (lapply(x, as.data.frame))))
    rownames(ret) <- NULL
    colnames(ret) <- as.character(1:ncol(ret))
  }
  return(as.data.frame(ret))
}


#' @title Calculate KL Divergence by Forcing Continuous Distribution to Matched Uniform
#' @param p probability vector p, considered the empiric, or "data", distribution
#' @param q probability vector q, considered the reference, or "model", distribution
#' @description Calculation is a distance between two probability distributions (\emph{N.B.} lacks symmetry).
#' Assumes use has inputted \emph{R >= Z_0} integer realizations of a test and a model that are drawn
#' from probability distributions P and Q, respectively. The code then calculates an empiric
#' distribution for each realization and compares their respective distance.
#' \emph{NB} Calculation was performed in log space and output metric is in nats (versus log-2 and bits)
#' @details This is not a statistically correct test because we force a continuous distribution to
#' a uniform/discrete state space and drop non-matching overlaps to avoid the issue of log(0)
#' @noMd
#' @importFrom magrittr %>%
#' @export

Cheap_KLdivergence_UnifDist <- function(p, q) {
  goodegg::assert_vector(p)
  goodegg::assert_vector(q)
  goodegg::assert_int(p)
  goodegg::assert_int(q)
  goodegg::assert_same_length(p, q)

  # identify/use empiric min and max
  lb <- min(c(p,q))
  ub <- max(c(p,q))
  probtib <- tibble::tibble(X = as.character(lb:ub))

  #......................
  # convert into prob dist
  #......................
  # get table
  ptab <- table(p)/length(p)
  ptab <- tibble::as_tibble(ptab) %>%
    dplyr::rename(X = p,
                  pFreq = n)
  qtab <- table(q)/length(q)
  qtab <- tibble::as_tibble(qtab) %>%
    dplyr::rename(X = q,
                  qFreq = n)
  # match
  probtib <- probtib %>%
    dplyr::left_join(., ptab, by = "X") %>%
    dplyr::left_join(., qtab, by = "X") %>%
    dplyr::filter(!is.na(pFreq)) %>%
    dplyr::filter(!is.na(qFreq))

  #......................
  # perform KL calculation
  #......................
  KLstat <- probtib$pFreq * (log(probtib$pFreq) - log(probtib$qFreq))
  KLstat <- sum(KLstat)
  # out
  return(KLstat)
}

