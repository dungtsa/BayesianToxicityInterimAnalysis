#' @name Toxicity

#' @title Toxicity Interim Analysis by Bayesian Posterior Probability

#' @description This function will start a "Toxicity Interim Analysis by Bayesian Posterior Probability" Shiny application

#' @export

ToxicityShinyApp <- function() {

  appDir <- system.file("shiny-examples", "myapp", package = "ToxicityAnalysis")

  if (appDir == "") {

    stop("Could not find example directory. Try re-installing `ToxicityAnalysis`.", call. = FALSE)

  }



  shiny::runApp(paste(appDir,'/app.R',sep=''), launch.browser =T,host = getOption( "127.0.0.1"))

}
