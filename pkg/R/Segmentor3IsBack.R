Segmentor<- function(data=numeric(), file=character(), model=1, Kmax = 15, param = 0.15) UseMethod("Segmentor")

Segmentor.default <-function(data=numeric(), file=character(), model=1, Kmax = 15, param = 0.15)
{
  if ((model!=1)&(model!=2)&(model!=3)&(model!=4))
    stop("Choose model=1 (Poisson), 2 (normal), 3 (Negative Binomial) or 4 (Normal-Variance)")
  if ((length(data)==0)&(length(file)==0))
    stop("Give me a vector of data to segment or a file to open")
  if ((length(data)>0)&(length(file)>0))
    stop("Give me either a vector of data or a file to open but not both")

  breaks=matrix(0,nrow=Kmax,ncol=Kmax)
  breaks = as.vector(breaks)
  parameters=matrix(0,nrow=Kmax,ncol=Kmax)
  parameters=as.vector(parameters)
  likelihood=rep(0, Kmax)
  likelihood=as.vector(likelihood)
  n = 0    
  if (length(file)>0)
  {
	data<-read.table(file)
	data<-data[,1]
  }
  n = length(data)
  if ((model==4) & (param==0.15))
	param = mean(data)
  if (model==1)
	Rep<-.C("SegmentPoisson", Size = as.integer(n),KMax = as.integer(Kmax), Data = as.integer(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), PACKAGE="Segmentor3IsBack") else if (model==3)
	Rep<-.C("SegmentBinNeg", Size = as.integer(n),KMax = as.integer(Kmax), theta = as.double(param), Data = as.integer(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), PACKAGE="Segmentor3IsBack") else if (model==2)
	Rep<-.C("SegmentNormal", Size = as.integer(n),KMax = as.integer(Kmax), Data = as.double(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), PACKAGE="Segmentor3IsBack") else if (model==4)
	Rep<-.C("SegmentVariance", Size = as.integer(n),KMax = as.integer(Kmax), theta = as.double(param), Data = as.double(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), PACKAGE="Segmentor3IsBack")

  breaks=matrix(Rep$Breakpoints,ncol=Kmax)
  breaks = t(breaks)
  parameters = matrix(Rep$Parameters,ncol=Kmax)
  parameters = t(parameters)
  rownames(breaks)<-c("1 segment", paste(2:Kmax, "segments"))
  rownames(parameters)<-c("1 segment", paste(2:Kmax, "segments"))
  colnames(breaks)<-c(paste(1:Kmax, "th break",sep=""))
  colnames(parameters)<-c(paste(1:Kmax, "th parameter",sep=""))
  likelihood=matrix(Rep$Likelihood,ncol=1)
  rownames(likelihood)<-c("1 segment", paste(2:Kmax, "segments"))
  if (model==1) 
  {
      model.dist="Poisson"
      Segmentor.res=list(model=model.dist,breaks=breaks,parameters=parameters,likelihood=likelihood)
  }
  if (model==2) 
  {
    variance <- rep(0,Kmax)
    variance[1] = variance[1] + breaks[1,1]*var(data[1:breaks[1,1]])/n
    for (j in 2:Kmax)
    {
	for (i in 2:j)
	{
	  var = (breaks[j,i]-breaks[j,i-1])*var(data[(breaks[j,i-1]+1):breaks[j,i]])/n
	  variance[j] = variance[j]+var
	}
	variance[j] = variance[j] + breaks[j,1]*var(data[1:breaks[j,1]])/n
    }
    likelihood = likelihood /(2*variance) + n/2*log(2*pi*variance)
    model.dist="normal"
    Segmentor.res=list(model=model.dist,breaks=breaks,parameters=parameters,likelihood=likelihood)
  }
  if (model==3) 
  {
      model.dist="Negative binomial"
      Segmentor.res=list(model=model.dist,breaks=breaks,parameters=parameters,likelihood=likelihood)
  }
  if (model==4) 
  {
      likelihood = likelihood + n/2*log(2*pi)
      model.dist="Variance Segmentation"
      Segmentor.res=list(model=model.dist,breaks=breaks,parameters=parameters,likelihood=likelihood)
  }
  class(Segmentor.res) <- "Segmentor"
  Segmentor.res

}

print.Segmentor <-function(x,...)
{
  cat("\n Model used for the segmentation: \n")
  print(x$model)
  
  cat("\n Table of optimal breakpoints: \n")
  print(x$breaks)

  cat("\n Table of negative log-likelihood for each optimal segmentation: \n")
  print(x$likelihood)
  
  if (is.element("parameters",names(x)))
  {
    cat("\n Table of parameters: ")
    if (x$model=="Poisson")
      cat("(mean of the signal in each segment) \n")
    if (x$model=="normal")
      cat("(mean of the signal in each segment) \n")
    if (x$model=="Negative binomial")
      cat("(success-probability of the signal in each segment) \n")
    if (x$model=="Variance Segmentation")
      cat("(Variance of the signal in each segment) \n")
    print(x$parameters)
  }
}

