# require r packages
library(shiny)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(DT)   
library(knitr)
library(rmarkdown)

# 4 R functions for toxicity
##########################
### function toxicity ###
##########################
DTC.bayesian.toxicity.fun<-function(
  p.prior.par, p.tox, block.n, rep.n, p.target, p.cutoff, sim.n)
{
  ans.list<-list()
  boudary.tmp<-numeric(0)
  if(!is.null(block.n))
  {
    #--determine the boundary---
    index1<-(1:(rep.n/block.n))*block.n
    for(mm in index1)
    {
      n1<-(0:mm)
      p.tmp<-pbeta(p.tox,n1+p.prior.par[1], (mm-n1)+p.prior.par[2])
      if(any(p.tmp<p.cutoff)) boudary.tmp<-rbind(boudary.tmp,c(mm,min(n1[p.tmp<p.cutoff]))) else boudary.tmp<-rbind(boudary.tmp,c(mm,NA))
    }
    
    for(j in 1:sim.n)
    {
      k1<-event.n<-p.ans<-numeric(0)
      for(i in 1:(rep.n/block.n))
      {
        k1<-c(k1,rbinom(block.n,1,p.target))
        event.n<-sum(k1)
        p.posterior.less.p_tox<-pbeta(p.tox,p.prior.par[1]+event.n,p.prior.par[2]+length(k1)-event.n)
        p.ans<-c(p.ans,p.posterior.less.p_tox)
      }
      ans.list[[j]]<-list(k1=k1,p.ans.less.tox=p.ans)
    }
  } else {
    
#--determine the boundary---
    index1<-cumsum(rep.n)
    
    for(mm in index1)
    {
      n1<-(0:mm)
      p.tmp<-pbeta(p.tox,n1+p.prior.par[1], (mm-n1)+p.prior.par[2])
      if(any(p.tmp<p.cutoff)) boudary.tmp<-rbind(boudary.tmp,c(mm,min(n1[p.tmp<p.cutoff]))) else boudary.tmp<-rbind(boudary.tmp,c(mm,NA))
    }
    
    for(j in 1:sim.n)
    {
      k1<-event.n<-p.ans<-numeric(0)
      for(i in rep.n)
      {
        k1<-c(k1,rbinom(i,1,p.target))
        event.n<-sum(k1)
        p.posterior.less.p_tox<-pbeta(p.tox,p.prior.par[1]+event.n,p.prior.par[2]+length(k1)-event.n)
        p.ans<-c(p.ans,p.posterior.less.p_tox)
      }
      ans.list[[j]]<-list(k1=k1,p.ans.less.tox=p.ans)
    }
  }
  
  a1<-sapply(ans.list,function(x) any(x$p.ans<p.cutoff))
  dimnames(boudary.tmp)[[2]]<-c('number of patients','boundary')
  list(prob=mean(a1),boudary=boudary.tmp,data=ans.list)
}

############################
### function sensitivity ###
############################
bayesian.toxicity.sensitivity.fun<-function(
  p.prior.par.tmp, p.tox.tmp, block.n.tmp, rep.n.tmp, p.target, p.cutoff, sim.n.tmp)
{
  b = seq(p.tox.tmp-0.1*4, p.tox.tmp+0.1*4, 0.1)
  p.target.list <- b[!b <= 0]
  
  p.cutoff.list <- seq(p.cutoff,p.cutoff+0.05*p.target, 0.05)
  
  ans.list.all<-list()
  
  for(i in 1:length(p.target.list))
  {
    p.target.tmp<-p.target.list[i]
    p.cutoff.tmp<-p.cutoff.list[1]
    ans.list.all[[i]]<-DTC.bayesian.toxicity.fun(p.prior.par=p.prior.par.tmp,p.tox=p.tox.tmp,block.n=block.n.tmp,rep.n=rep.n.tmp, p.target=p.target.tmp,  p.cutoff=p.cutoff.tmp ,sim.n=sim.n.tmp)[3]
  }
  names(ans.list.all)<-p.target.list 
  
  tmp<-numeric(0)
  tmp1<-list()
  for(i in 1:length(p.cutoff.list))
  {
    tmp<-rbind(tmp,sapply(ans.list.all,function(y) {a1<-sapply(y$data,function(x) {x1<-x$p.ans.less.tox;ans<-(1:length(x1))[cumsum(x1<p.cutoff.list[i])==1];ifelse(length(ans)==0,0,ans) });a1<-factor(a1,level=0:ifelse(is.null(block.n.tmp),length(rep.n.tmp),rep.n.tmp));  mean(a1!='0')} ))
  }
  for(i in 1:length(p.cutoff.list))
  {
    tmp1[[i]]<-sapply(ans.list.all,function(y) {a1<-sapply(y$data,function(x) {x1<-x$p.ans.less.tox;ans<-(1:length(x1))[cumsum(x1<p.cutoff.list[i])==1];ifelse(length(ans)==0,0,ans) });a1<-factor(a1,level=0:ifelse(is.null(block.n.tmp),length(rep.n.tmp),rep.n.tmp)); table(a1)/length(a1)}) 
  }
  dimnames(tmp)[[1]]<-p.cutoff.list
  names(tmp1)<-p.cutoff.list
  
  # find more than 80% value for 30% true toxicity
  col_n <- which(p.target.list == p.tox.tmp)
  row_n <- which(p.cutoff.list == p.cutoff)
  
  value <- tmp[row_n,col_n]*100
  
 matplot_s <- as.data.frame(tmp) %>% 
    rownames_to_column("Cutoff") %>%  
    mutate(Cutoff = factor(Cutoff)) %>% 
    gather(rate,value1, - Cutoff) %>% 
    mutate(value = value1 *100) %>% 
    ggplot(aes(x = rate, y = value, col = Cutoff, group = Cutoff)) +
    geom_point(aes(shape=Cutoff), size = 3) +
    scale_shape_manual(values=c(15:25)) +
    geom_line(aes(linetype=Cutoff), size = 1) +
    labs(title = "Sensitivity Analysis of Toxicity", x = "True Toxicity Rate", y = "Probability of Stopping the Trial (%)") +
    theme_bw() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 16))
  
 table2 <- as.data.frame(tmp) %>% 
   rownames_to_column("True Toxicity Rate") %>%  
   mutate_at(vars(-`True Toxicity Rate`),funs(. * 100)) %>% 
   filter(`True Toxicity Rate` == p.cutoff)
 table2[1,1] <- paste('Cutoff =', p.cutoff, sep = " ")
 
  results <- list(matplot_s = matplot_s,
                  value = value,
                  table2 = table2)
  return(results)
  #legend(0.1,1,expression(phi))
  #list(ind=ans.list.all,cutoff=tmp,cutoff.ind=tmp1)
}

##########################
### function boundary ###
##########################

tox.boundary.fun<-function(
  p.prior.par, p.tox, block.n, rep.n, p.cutoff)
{
  boudary.tmp<-numeric(0)
  if(!is.null(block.n))
  {
    #--determine the boundary---
    index1<-(1:(rep.n/block.n))*block.n
    for(mm in index1)
    {
      n1<-(0:mm)
      p.tmp<-pbeta(p.tox,n1+p.prior.par[1], (mm-n1)+p.prior.par[2])
      if(any(p.tmp<p.cutoff)) boudary.tmp<-rbind(boudary.tmp,c(mm,min(n1[p.tmp<p.cutoff]))) else boudary.tmp<-rbind(boudary.tmp,c(mm,NA))
    }
  }else {
    
    #--determine the boundary---
    index1<-cumsum(rep.n)
    
    for(mm in index1)
    {
      n1<-(0:mm)
      p.tmp<-pbeta(p.tox,n1+p.prior.par[1], (mm-n1)+p.prior.par[2])
      if(any(p.tmp<p.cutoff)) boudary.tmp<-rbind(boudary.tmp,c(mm,min(n1[p.tmp<p.cutoff]))) else boudary.tmp<-rbind(boudary.tmp,c(mm,NA))
    }
    
  }
  boudary.tmp
}

##########################
### function tmp_tab ###
##########################


tmp_tab <- function(p.prior.par, p.tox, block.n, rep.n, p.cutoff, p.tmp){
  
  tmp1<-tox.boundary.fun(p.prior.par, p.tox, block.n, rep.n, p.cutoff)
  
  for(i in 2:length(p.tmp)) {
    tmp1<-cbind(tmp1,tox.boundary.fun(p.prior.par, p.tox, block.n, rep.n, p.cutoff=p.tmp[i])[,2])
  }
  dimnames(tmp1)[[2]]<-c('n',p.tmp)
 
  tab <- as.data.frame(tmp1)
  
  sub_tab <- as.data.frame(t(tab[,1:2])) %>% 
    rownames_to_column("No.stage") 
  colnames(sub_tab) <- c("No.stage",seq(1:dim(tab)[1]))  
  cutoff_value <- sub_tab[2,1]
  sub_tab[2,1] <- paste('Cutoff =', cutoff_value, sep = " ")
  
  tab1 <-  tab[round(dim(tab)[1]/2),1:2] 
  
  tab2 <- list(tab = tab, sub_tab = sub_tab,pats_n1 = tab1[,1], pats_n2 = tab1[,2])
  
  return(tab2)
}


############################
### Shiny UI and Server ###
############################

# Define UI ####
# for application that plots sensitivity analysis of toxicity
ui <- fluidPage( titlePanel("Toxicity Interim Analysis by Bayesian Posterior Probability"),
                 
                 # Sidebar layout with a input and output definitions
                 sidebarLayout(
                   
                   # Inputs
                   
                   sidebarPanel(
                     
                     hr(),
                     
                     # Set toxicity rate
                     sliderInput(inputId = "tox", 
                                 label = "Toxicity Rate", 
                                 min = 0, max = 1, 
                                 value = 0.17),
                     hr(),
                     
                     # Set beta level
                     
                     numericInput(inputId = "betaA", 
                                 label = "Beta A parameter", 
                                 min = 0, max = Inf, 
                                 value = 0),
                     hr(),
                     
                     numericInput(inputId = "betaB", 
                                 label = "Beta B parameter", 
                                 min = 0, max = Inf, 
                                 value = 1),
                     hr(),
                     
                     # Set block number 
                     numericInput(inputId = "block_n", 
                                 label = "Number patients per interim analysis", 
                                  min = 0, max = 100, 
                                  value = 1),
                    hr(),
                    
                    # Set repeat number
                     numericInput(inputId = "rep", 
                                label = "Total sample size", 
                                min = 10, max = 100, 
                                value = 14),    
                     
                    hr(),
                   # Select p. target list number
                   numericInput(inputId = "p_target", 
                       label = "Length of cutoffs of the posterior probability", 
                       min = 2, max = 10, 
                       value = 6),
                   hr(), 
                  # Input p. cutoff value
                   numericInput(inputId = "p_cutoff", 
                        label = "Cutoffs of the posterior probability for toxicity rate", 
                        min = 0, max = 0.5, 
                        value = 0.05),
                  
                  hr(),
                  # Input numer of simulation
                  numericInput(inputId = "simN", 
                               label = "number of simulation", 
                               value = 100),
                 
                   hr(),
                   actionButton("go", "Submit"),
                   radioButtons('format', 'Document format', c('Word', 'PDF'),
                                  inline = TRUE),
                   downloadButton('downloadReport')
                  
                   ),
                   
                   # Output:
                   mainPanel(
                       tabPanel(br(),
                                h4(uiOutput(outputId = "text1")),
                                br(),
                                h4(uiOutput("tab1")),
                                tableOutput(outputId = "table1"),
                                br(),
                                h4(uiOutput("tab2")),
                                tableOutput(outputId = "table2"),
                                br(),
                                h4(uiOutput("fig1")),
                                plotOutput(outputId = "plot"),
                                br(),
                                h4(uiOutput(outputId = "tab3")),
                                dataTableOutput(outputId = "table3")
                              )
                   )
                 )
)

# Define server to summarize and view the results
server <- function(input, output) {
  
  # Show plot
  # Reactive expression to compose a data frame containing all of the values
  get.result.toxicity <- eventReactive(input$go,{
    req(input$betaA)
    req(input$tox) # ensure availablity of value before proceeding
    
    sens <- bayesian.toxicity.sensitivity.fun(
      p.prior.par.tmp = c(input$betaA,input$betaB),
      p.tox.tmp = input$tox,
      block.n.tmp = input$block_n,
      rep.n.tmp = input$rep,
      p.target = input$p_target,
      p.cutoff = input$p_cutoff,
      sim.n.tmp=input$simN
      )
    
    tab <- tmp_tab(
      p.prior.par = c(input$betaA,input$betaB),
      p.tox = input$tox,
      block.n = input$block_n,
      rep.n = input$rep,
      p.cutoff = input$p_cutoff,
      p.tmp = seq(input$p_cutoff,input$p_cutoff + 0.05 * input$p_target, 0.05))
    
   
    ls <- list(df = tab$tab,
               table1 = tab$sub_tab,
               pats_n1 = tab$pats_n1,
               pats_n2 = tab$pats_n2,
               plot = sens$matplot_s,
               value = sens$value,
               table2 = sens$table2)
    return(ls)
    })
 
  # Print statistics plan 
   text1 <- eventReactive(input$go,{
    tox.tmp <- paste(round(input$tox*100, 2),'%', sep = '')
    txt <- paste("Statistics Plan : The study will enroll a total number of",input$rep ,"patients. 
                 We will evaulate toxicity at every",input$block_n," patient(s). We consider ",tox.tmp," as
                 the maximum allowable toxicity rate.  Using Bayesian posterior probability with a beta prior,
                 beta (", input$betaA,",",input$betaB,"), if the posterior probability of toxicity rate  < ",tox.tmp,"
                 is less than ",input$p_cutoff,", we will conclude the treatment is too toxic and should be terminated.
                 For example (Table 1), if ", get.result.toxicity()$pats_n2," patients or more have unacceptable toxicity
                 in ", get.result.toxicity()$pats_n1,"  patients (about 50%) enrolled and evaluated, we consider the
                 treatment is too toxic and the study will be halted. Sensitivity analysis (Table 2) shows that if 
                 the true toxic rate is ",tox.tmp,", the chance of early stopping the trial is ",
                 get.result.toxicity()$value,"% based on ", input$simN, "times simulation with the
                 posterior probability cutoff value at",input$p_cutoff,". Evaluation of different cutoffs of
                 the posterior probability for toxicity rate is given in the Figure 1 for the sensitivity 
                 analysis and in Table 3 for the stopping boundary analysis.", sep = " ")
  })
  
  
  text2 <- eventReactive(input$go,{
    txt2 <- paste("Table 1. The stopping boundary for cutoff at", input$p_cutoff, sep = " ") 
  })
  
  text3 <- eventReactive(input$go,{
    txt3 <- paste("Table 2. The posterior probability of stopping the trial for cutoff at", input$p_cutoff, sep = " ") 
  })
  
  text4 <- eventReactive(input$go,{
    txt4 <- "Figure 1. Sensitivity analysis at the different listed cutoffs of the posterior probability for toxicity rate"
  })
  
  text5 <- eventReactive(input$go,{
    txt5 <- "Table 3.The stopping boundary at the different listed cutoffs of the posterior probability for toxicity rate"
  })
  
  output$text1 <- renderUI({
    HTML(text1())
  })
  
  output$tab1 <- renderUI({
    HTML(text2())
  })
  
  output$tab2 <- renderUI({
    HTML(text3())
  })
  
  output$fig1 <- renderUI({
    HTML(text4())
  })
  
  output$tab3 <- renderUI({
    HTML(text5())
  })
  
  # Show table1 
  output$table1  <- renderTable({
    get.result.toxicity()$table1}, digits = 0, rownames = TRUE)
  
  
  # Show table2 
  output$table2  <- renderTable({
    get.result.toxicity()$table2}, digits = 0)
 
  
    #Show plot
  output$plot <- renderPlot({
    get.result.toxicity()$plot
  })
  
  #show table 3
  output$table3  <- renderDataTable({
    datatable(data = get.result.toxicity()$df,
              options = list(pageLength = 10, lengthMenu = c(10, 25, 40)), 
              rownames = FALSE)
  }) 
  
  
  #---output report----
  
  output$downloadReport <- downloadHandler(
    
    filename = function() {
      
      paste('ToxicityReport', sep = '.',
          switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'))
      
    },
    
    content = function(file) {
      
      out <- rmarkdown::render(input = 'ToxicityReportMarkdownFile.Rmd',
                               
                               output_format =
                                 switch(input$format,
                                        PDF = rmarkdown::pdf_document(),
                                        HTML = rmarkdown::html_document(),
                                        Word = rmarkdown::word_document()
                                        
                                 ),
                               
                  params = list(set_title = input$project_title, set_author = input$author_input)
                               
      )
      
      #file.rename(out, file)
       file.copy(out,file)
      
    }
    
  )
  
}  
  



# Run the application 
shinyApp(ui = ui, server = server)
