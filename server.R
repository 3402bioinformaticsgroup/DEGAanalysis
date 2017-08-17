# server.R
#install.packages("kernlab")
library(shiny)
library(pwr)

shinyServer(function(input, output){

output$input.view<-renderTable({
        files<-input$NanoStringData
        if(is.null(files))
        return(NULL)
        data<-as.matrix(read.csv(files$datapath, header=T, row.names=1))
        as.matrix(colnames(data))
    })
# plot power
output$plotPower<-renderPlot({
    library(pwr)
    files<-input$NanoStringData
    if(is.null(files))
    return(NULL)
    data<-as.matrix(read.csv(files$datapath, header=T, row.names=1))
    
    ref.file<-input$PanelReference   # read in panel reference
    if(is.null(ref.file))
    return(NULL)
    ref<-as.matrix(read.csv(ref.file$datapath, header=T))
    CodeClass<-ref[,1]
    logFC<-as.numeric(input$LogFC.cutoff) # logFC.cutoff
    FC<-2^logFC                           # FC.cutoff
    positive<-data[grep("Positive", CodeClass),]   #not necessary here
    data<-data[grep("Endogenous", CodeClass),]     #Cut data, leave only probes
    table<-matrix(rep(0,dim(data)[1]*12), ncol=12)
    rownames(table)<-rownames(data)
    colnames(table)<-c("mean.1", "mean.2", "sd.1", "sd.2",
    "sample size(n)","effective size (d)","power",
    "Log2FC", "p-value", "sample size(N,power=0.8,d=defined)",
    "CV.1","CV.2") # CVs are in absolute level
    bio.condition1<-as.numeric(unlist(strsplit(input$Condition1,","))) # define condition 1
    bio.condition2<-as.numeric(unlist(strsplit(input$Condition2,","))) # define condition 2
    
    for (i in 1:dim(table)[1]){
        table[i,1] <- mean(data[i,bio.condition1]) #mean.1
        table[i,2] <- mean(data[i,bio.condition2]) #mean.2
        table[i,3] <- sd(data[i,bio.condition1])   #sd.1
        table[i,4] <- sd(data[i,bio.condition2])   #sd.2
        table[i,5] <- length(data[i,bio.condition1])+length(data[i,bio.condition2]) # n
        table[i,6] <- abs(table[i,1]-table[i,2])/ (((table[i,3])^2 + (table[i,4])^2)/2)^0.5 # d
        if(table[i,6]!="NaN"){
            if (table[i,6]>0.001){
            if (table[i,6]<5){table[i,7] <- pwr.t.test(n = table[i,5], d = table[i,6], sig.level = 0.05, type ="two.sample")$power}
          }
        }
        table[i,8] <- log(table[i,2]+0.01,2)-log(table[i,1]+0.01,2)    # log2FC
        if(table[i,6]!="NaN"){
        table[i,9] <- t.test(data[i,bio.condition2],data[i,bio.condition1])[[3]] #p-value
        }
        if(table[i,6]!="NaN"){
            if (table[i,6]>0.001){
            if (table[i,6]<5){table[i,10] <- pwr.t.test(d = table[i,6], power=0.8, sig.level = 0.05, type ="two.sample")$n}
          }
        }
        if(table[i,3]>0){table[i,11]<-table[i,3]/table[i,1]}
        if(table[i,4]>0){table[i,12]<-table[i,4]/table[i,2]}
    }
    mean.cv<-mean(c(mean(table[,11]), mean(table[,12])))
    mean.mean<-mean(c(mean(table[,1]), mean(table[,2])))
    # range of logFC
    r <- seq(0.3, 3, 0.1)
    nr <- length (r)
    # power values
    p<- seq(0.4, 0.9, 0.1)
    np<-length(p)
    # obtain sample sizes
    samsize <- array(numeric(nr*np), dim=c(nr,np))
    for (i in 1:np){
        for (j in 1:nr){     #abs(table[i,1]-table[i,2])/ (((table[i,3])^2 + (table[i,4])^2)/2)^0.5
            
            tryCatch({
                
            effective.size<-abs(mean.mean*2^r[j]-mean.mean)/((((mean.mean*2^r[j])*mean.cv)^2 + (mean.mean*mean.cv)^2)/2)^0.5
            #print(effective.size)
            if(effective.size!="NaN"){
            if(effective.size > 0.001){
                if (effective.size < 5){
            result <- pwr.t.test(n = NULL, d = effective.size,sig.level = .05, power = p[i],type ="two.sample")
                    samsize[j,i] <- ceiling(result$n)
                }
               }
              }
          
          
            }, error=function(e){})
            
        }
    }
    # set up graph
    xrange <- range(r)
    yrange <- round(range(samsize))
    colors <- rainbow(length(p))
    plot(xrange, yrange, type="n",
    xlab="logFC",
    ylab="Sample Size (n)" )
    
    # add power curves
    for (i in 1:np){
        points(r, samsize[,i], type="l", lwd=2, col=colors[i])
    }
    
    # add annotation (grid lines, title, legend)
    #abline(v=0, h=seq(0,yrange[2],50), lty=2, col="grey89")
    #abline(h=0, v=seq(xrange[1],xrange[2],.02), lty=2,col="grey89")
    title("Sample Size Estimation\n
    Sig=0.05")
    legend("topright", title="Power", as.character(p),
    fill=colors)

    
    })
# plot DEG
output$plotDEG<-renderPlot({
    library(pwr)
    files<-input$NanoStringData
    if(is.null(files))
    return(NULL)
    data<-as.matrix(read.csv(files$datapath, header=T, row.names=1))
    
    ref.file<-input$PanelReference   # read in panel reference
    if(is.null(ref.file))
    return(NULL)
    ref<-as.matrix(read.csv(ref.file$datapath, header=T))
    CodeClass<-ref[,1]
    logFC<-as.numeric(input$LogFC.cutoff) # logFC.cutoff
    FC<-2^logFC                           # FC.cutoff
    positive<-data[grep("Positive", CodeClass),]   #not necessary here
    data<-data[grep("Endogenous", CodeClass),]     #Cut data, leave only probes
    table<-matrix(rep(0,dim(data)[1]*10), ncol=10)
    rownames(table)<-rownames(data)
    colnames(table)<-c("mean.1", "mean.2", "sd.1", "sd.2",
                       "sample size(n)","effective size (d)","power",
                        "Log2FC", "p-value", "sample size(N,power=0.8,d=defined)")
    bio.condition1<-as.numeric(unlist(strsplit(input$Condition1,","))) # define condition 1
    bio.condition2<-as.numeric(unlist(strsplit(input$Condition2,","))) # define condition 2
    
    for (i in 1:dim(table)[1]){
        table[i,1] <- mean(data[i,bio.condition1]) #mean.1
        table[i,2] <- mean(data[i,bio.condition2]) #mean.2
        table[i,3] <- sd(data[i,bio.condition1])   #sd.1
        table[i,4] <- sd(data[i,bio.condition2])   #sd.2
        table[i,5] <- length(data[i,bio.condition1])+length(data[i,bio.condition2]) # n
        table[i,6] <- abs(table[i,1]-table[i,2])/ (((table[i,3])^2 + (table[i,4])^2)/2)^0.5 # d
        if(table[i,6]!="NaN"){
        if (table[i,6]>0.001){
            if (table[i,6]<5){table[i,7] <- pwr.t.test(n = table[i,5], d = table[i,6], sig.level = 0.05, type ="two.sample")$power}
        }
        }
        table[i,8] <- log(table[i,2]+0.01,2)-log(table[i,1]+0.01,2)    # log2FC
        if(table[i,6]!="NaN"){
        table[i,9] <- t.test(data[i,bio.condition2],data[i,bio.condition1])[[3]] #p-value
        }
        if(table[i,6]!="NaN"){
        if (table[i,6]>0.001){
            if (table[i,6]<5){table[i,10] <- pwr.t.test(d = table[i,6], power=0.8, sig.level = 0.05, type ="two.sample")$n}
        }
        }
    }
    ToPlot<-table
    color<-rep("gray", dim(ToPlot)[1])
    for (i in 1:dim(ToPlot)[1]){
        if(ToPlot[i,9]<0.05){color[i]<-"red"}
    }
    plot(log(((ToPlot[,1]+ToPlot[,2])/2)+0.01,2), ToPlot[,8], col=color,pch=16, main="Differentially Expressed Genes", xlab="LogMean", ylab="LogFC", ylim=c(min(ToPlot[,8]), max(ToPlot[,8])))
    abline(h=c(logFC, -1*logFC), lty=2, col="red")
    
    })
# Download
output$downloadData <- downloadHandler(
filename = function() { paste(input$outputfile, '.csv', sep='') },
content = function(file) {
    library(pwr)
    files<-input$NanoStringData
    if(is.null(files))
    return(NULL)
    data<-as.matrix(read.csv(files$datapath, header=T, row.names=1))
    
    ref.file<-input$PanelReference   # read in panel reference
    if(is.null(ref.file))
    return(NULL)
    ref<-as.matrix(read.csv(ref.file$datapath, header=T))
    CodeClass<-ref[,1]
    logFC<-as.numeric(input$LogFC.cutoff) # logFC.cutoff
    FC<-2^logFC                           # FC.cutoff
    positive<-data[grep("Positive", CodeClass),]   #not necessary here
    data<-data[grep("Endogenous", CodeClass),]     #Cut data, leave only probes
    table<-matrix(rep(0,dim(data)[1]*10), ncol=10)
    rownames(table)<-rownames(data)
    colnames(table)<-c("mean.1", "mean.2", "sd.1", "sd.2",
    "sample size(n)","effective size (d)","power",
    "Log2FC", "p-value", "sample size(N,power=0.8,d=defined)")
    bio.condition1<-as.numeric(unlist(strsplit(input$Condition1,","))) # define condition 1
    bio.condition2<-as.numeric(unlist(strsplit(input$Condition2,","))) # define condition 2
    
    for (i in 1:dim(table)[1]){
        table[i,1] <- mean(data[i,bio.condition1]) #mean.1
        table[i,2] <- mean(data[i,bio.condition2]) #mean.2
        table[i,3] <- sd(data[i,bio.condition1])   #sd.1
        table[i,4] <- sd(data[i,bio.condition2])   #sd.2
        table[i,5] <- length(data[i,bio.condition1])+length(data[i,bio.condition2]) # n
        table[i,6] <- abs(table[i,1]-table[i,2])/ (((table[i,3])^2 + (table[i,4])^2)/2)^0.5 # d
        if(table[i,6]!="NaN"){
        if (table[i,6]>0.001){
            if (table[i,6]<5){table[i,7] <- pwr.t.test(n = table[i,5], d = table[i,6], sig.level = 0.05, type ="two.sample")$power}
        }
        }
        table[i,8] <- log(table[i,2]+0.01,2)-log(table[i,1]+0.01,2)    # log2FC
        if(table[i,6]!="NaN"){
        table[i,9] <- t.test(data[i,bio.condition2],data[i,bio.condition1])[[3]] #p-value
        }
        if(table[i,6]!="NaN"){
        if (table[i,6]>0.001){
            if (table[i,6]<5){table[i,10] <- pwr.t.test(d = table[i,6], power=0.8, sig.level = 0.05, type ="two.sample")$n}
        }
        }
    }
     write.csv(table,file)

      })

})