#............................................................
# SIR NE S3 Class Overloading
#...........................................................
#' @title Check if GillespieSIRne S3 Class
#' @description Overload is: function for determining if object is of class GillespieSIRne
#' @param x GillespieSIRne Simulation
#' @noMd
#' @export
is.GillespieSIRne <- function(x) {
  inherits(x, "is.GillespieSIRne")
}

#' @title print GillespieSIRne S3 Class
#' @description overload print() function to print summary only
#' @inheritParams is.GillespieSIRne
#' @param ... further arguments passed to or from other methods.
#' @noMd
#' @export
print.GillespieSIRne <- function(x, ...) {

  # print summary only
  cat(crayon::red("Final Epidemic Size:"), summary(x)$FinalEpidemicSize, "\n")
  cat(crayon::blue("Run Time:"), x$runTime, attr(x$runTime, "units"), "\n")

}

#' @title Summary of GillespieSIRne S3 Class
#' @description overload summary() function.
#' @param object GillespieSIRne Simulation
#' @param ... further arguments passed to or from other methods.
#' @noMd
#' @export
summary.GillespieSIRne <- function(object, ...) {
  # send summary only
  return(tidyout.GillespieSIRne(object))
}


#............................................................
# organzing SIR NE output
#...........................................................
#' @title Tidy Out Sim Method
#' @description Method assignment
#' @inheritParams is.GillespieSIRne
#' @noMd
#' @export
tidyout <- function(x) {
  UseMethod("tidyout")
}

#' @title Tidy Out Sim
#' @description Function for taking output of SIR NE and lifting it over
#' @inheritParams is.GillespieSIRne
#' @noMd
#' @export

tidyout.GillespieSIRne <- function(x) {
  #......................
  # ind lvl
  #......................

  #......................
  # summaries
  #......................
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
#' @details Internal function; not meant for general use
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
