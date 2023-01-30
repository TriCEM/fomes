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
  finalsize <- sum(x$Event_traj == "transmission")
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
