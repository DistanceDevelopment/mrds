
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)
source("helper.r")
shinyServer(function(input,output,session) {
  # Load observation data from file; update width with max distance
  data <- reactive({
    maxcut <- 0
    if (is.null(input$datafile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    else
    {
      obs<-read.delim(input$datafile$datapath)
      left <- 0
      if(input$disttype!="Perpendicular")
      {
        validate(
          need(!is.null(obs$radial),"Missing radial distance. Must be named radial."),
          need(!is.null(obs$angle),"Missing sighting angle. Must be named angle.")
        )
        obs$distance=obs$radial*sin(obs$angle*pi/180)
      } else 
        validate(
             need(!is.null(obs$distance),"Missing perpendicular distance. Must be named distance."))
      validate(
          need(!is.null(obs$object),"field in observations file named object is missing"),
          need(!is.null(obs$observer),"field in observations file named observer is missing"),
          need(!is.null(obs$detected),"field in observations file named detected is missing"),
          need(!is.null(obs$size) | input$obstype=="Single","field in observations file named size is missing")          
        )
      if (input$dobins == 0)
      {
        if(!is.null(obs$distbegin)&!is.null(obs$distend))
        {
          cutpoints=sort(unique(c(obs$distbegin,obs$distend)))
          updateCheckboxInput(session,"binned",value=TRUE)
          updateTextInput(session,"cutpoints",value=paste("c(",paste(cutpoints,collapse=","),")",sep=""))
          left <- cutpoints[1]
          maxcut <- cutpoints[length(cutpoints)]
        }        
      } else
      {
        if(input$binned)
        {
          input$dobins
          isolate({
            if(input$cutpoints!="")
            {
              cutpoints=eval(parse(text=input$cutpoints))
              cutpoints=cutpoints[order(cutpoints)]
              if(any(obs$distance>cutpoints[length(cutpoints)]))
              {
                width <- cutpoints[length(cutpoints)]
                obs<- obs[obs$distance<=cutpoints[length(cutpoints)],]
              }
              if(any(obs$distance<cutpoints[1]))
              {
                left <- cutpoints[1]
                obs<- obs[obs$distance>=left,]
              } 
              int=cut(obs$distance,cutpoints,right=FALSE)
              obs$distbegin=cutpoints[as.numeric(int)]
              obs$distend=cutpoints[as.numeric(int)+1]
              maxcut <- cutpoints[length(cutpoints)]
            }})
        }        
      }
      if(input$width==0)
        width<-max(obs$distance)
      else
      {
        width<-input$width
        width <- max(width, maxcut)
      }
      obs <- obs[obs$distance<width,]          
      updateNumericInput(session,"width",value=width)        
    }
    if(input$obstype=="Single" & !is.null(obs$size))obs$size <- NULL
    return(list(obs=obs,width=width,left=left))
  })
  # Render to data table
  output$datatable <- renderDataTable(data()$obs)
  # Load observation data from file; update width with max distance
  regions <- reactive({
    if (is.null(input$regionfile)) {
      # User has not uploaded a file yet
      return(NULL)
    } else
    {
      regions <-read.delim(input$regionfile$datapath)
      validate(
        need(!is.null(regions$Region.Label),"field in regions file named Region.Label is missing"),
        need(!is.null(regions$Area),"field in regions file named Area is missing")
      )
      return(regions)
    }
  })
  samples <- reactive({
    if (is.null(input$samplefile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    else
    {
      samples <-read.delim(input$samplefile$datapath)
      validate(
        need(!is.null(samples$Region.Label),"field in samples file named Region.Label is missing"),
        need(!is.null(samples$Sample.Label),"field in samples file named Sample.Label is missing"),
        need(!is.null(samples$Effort),"field in samples file named Effort is missing")
      )  
      return(samples)
    }})  
  obs <- reactive({
    if (is.null(input$linkfile)) {
      return(NULL)
    }
    else
    {
      obs <-read.delim(input$linkfile$datapath)
      validate(
        need(!is.null(obs$Region.Label),"field in obs-link file named Region.Label is missing"),
        need(!is.null(obs$Sample.Label),"field in obs-link file named Sample.Label is missing"),
        need(!is.null(obs$object),"field in obs-link file named object is missing")
      )  
      return(obs)
    }})  
  # Once go has been pressed, fit model
  results <- reactive({
    if (input$goButton == 0)
       invisible()
    else
    {
      input$goButton
      validate(
        need(!input$binned | (!is.null(data()$obs$distbegin) & !is.null(data()$obs$distend)),
             "Missing distance intervals for binned data. Must be named distbegin and distend."))
      isolate({
        if(is.null(data()))invisible()
        df=data()$obs
        if(input$numobs=="Single")mr.formula=NULL
        if(input$numobs=="Double" & input$indeptype=="Full")
        {
          key <- NULL
          adj <- NULL
        }else
        {
          key<-switch(input$keyfct,
                      "Half-normal"="hn",
                      "Hazard rate"="hr",
                      "Uniform"="unif",
                      "Gamma"="gamma")
          adj<-switch(input$adjfct,
                      "Cosine"="cos",
                      "Simple Polynomial"="poly",
                      "Hermite Polynomial"="herm",
                      "None"=NULL)          
        }
        if(input$type=="Point")
          meta.data="point=TRUE"
        else
          meta.data="point=FALSE"
        meta.data <- paste(meta.data,",width=",data()$width,sep="")
        meta.data <- paste(meta.data,",left=",data()$left,sep="")
        if(input$binned)meta.data <- paste(meta.data,",binned=TRUE")
        if(input$binned)meta.data <- paste(meta.data,",breaks=c(", paste(input$cutpoints,collapse=","),")",sep="")
        method="ds"
        if(input$numobs=="Double")
        {
          method=ifelse(input$indeptype=="Point","",".fi")
          method=paste(ifelse(input$config=="Independent","io",ifelse(input$config=="Trial","trial","rem")),method,sep="")
        }
        ddfstring <- ddf_string(key,adj,input$order,input$scale.formula,meta.data,method,input$mr.formula)
        return(list(model=eval(parse(text=ddfstring)),ddfstring=ddfstring,method=method))
      })
    }
  })
  # Reactive to number of bins, compute gof
  gof <- reactive({
    if(is.null(data()))invisible()
    isolate({width<-data()$width})
    if(!input$binned)
      bins<-seq(0, width, length.out = input$bins + 1)
    else
      bins <- eval(parse(text=input$cutpoints))
    gof <- ddf.gof(results()$model,qq=FALSE,breaks=bins)
    qqresults <- qqplot.ddf(results()$model,plot=FALSE)
    list(bins=bins,gof=gof,qqresults=qqresults)
  })
  # Reactive to number of bins, produce and output plot
  output$distPlot <- renderPlot({
    if (input$goButton == 0) return()
    else
    {
        if(is.null(data()))invisible()
        if(input$numobs=="Double")
        {
          if(input$config=="Independent")
            par(mfrow=c(2,3))
          else
            if(input$config=="Removal")
              par(mfrow=c(2,2))
            else
              par(mfrow=c(1,2))          
        }
        plot(results()$model,breaks=gof()$bins)
    }
  })
  # Once go button pushed output gof chi-square table
  output$gofvalues <-   reactive({
    if (input$goButton == 0)
      return("No model fit yet")
    else
    {
      if(is.null(data()))invisible()
        paste(capture.output(print(gof()$gof)),collapse="\n")
    }
   
  })
  output$qqplot <- renderPlot({
    if (input$goButton == 0 | input$binned)
      return()
    else
    { 
       plot(gof()$qqresults$edf[,2],gof()$qqresults$cdf,xlab="Empirical cdf",ylab="Fitted cdf",
         xlim=c(0,1),ylim=c(0,1))
       abline(0,1)
    }
  })  
  output$gof <- renderDataTable(
  {
  if (input$goButton == 0)
    return()
  else
  {
    if(is.null(data()))invisible()
    t(chitable(gof()$gof$chisquare$chi1$observed, gof()$gof$chisquare$chi1$expected))                
  }
})
# Once go button pushed output gof chi-square P-value
output$qqstats <- renderText( 
{
  if (input$goButton == 0)
    return("")
  else
  {
    if(is.null(data()))invisible()
    if(input$binned)
       "No test for binned data."
    else
       paste("The Kolmogorov-Smirnov test statistic is ",sprintf("%.3f",gof()$qqresults$ks$Dn,digits=3)," (P = ", sprintf("%.3f",gof()$qqresults$ks$p),"). The Cramer-von Mises test statistic (unweighted) is ", sprintf("%.3f",gof()$qqresults$CvM$W), " (P = ",sprintf("%.3f",gof()$qqresults$CvM$p),")",sep="")
 }})
output$P <- renderText( 
{
  if (input$goButton == 0)
    return("")
  else
  {
    if(is.null(data()))invisible()
    if(!is.na(gof()$gof$chisquare$chi1$p)){
      paste("\nP =",format(gof()$gof$chisquare$chi1$p,digits=5),
            " with ",gof()$gof$chisquare$chi1$df," degrees of freedom\n",sep="")
    }else{
      "\nNo degrees of freedom for test\n"
    }
  }})
output$fit <- reactive( 
{
  if (input$goButton == 0)
    return("No model fitted yet")
  else
  {
    input$goButton
    isolate({
      if(is.null(data()))invisible()
      paste(capture.output(print(results()$model)),collapse="\n")
    })
  }
})
output$abundance <- renderText( 
{
  if (input$goButton == 0)
    return("No model fitted yet")
  else
  {
    input$goButton
#    isolate({
      validate(
        need(!is.null(obs()),"Missing observation links file."),
        need(!is.null(samples()),"Missing samples file"),
        need(!is.null(regions()),"Missing regions file")
      )  
      paste(capture.output(print(dht(results()$model,regions(),samples(),obs()))),collapse="\n")
  #})
  }
})
observe({
  if(input$stopButton==0)
    return()
  isolate({
    obs <- data()$obs
    rownames(obs) <- NULL
    binstring <- paste("c(",paste(gof()$bins,collapse=","),")",sep="")
    stopApp(list(data=c(paste(names(obs),collapse=","),apply(obs,1, paste, collapse=",")),ddfstring=results()$ddfstring,binstring=binstring))
  })         
})
})

