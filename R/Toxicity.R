#' @name ToxicityShinyApp

#' @title Toxicity Interim Analysis by Bayesian Posterior Probability

#' @description This function will start a "Toxicity Interim Analysis by Bayesian Posterior Probability" Shiny application

#' @export

ToxicityShinyApp<-function () 
{
  appDir <- system.file("shiny-examples", "myapp", package = "BayesianToxicityInterimAnalysis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `BayesianToxicityInterimAnalysis`.", 
         call. = FALSE)
  }
  shiny::runApp(paste(appDir, "/app.R", sep = ""), launch.browser = T, 
                host = getOption("127.0.0.1"))
}
