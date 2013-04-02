Segmentor<- function(data=numeric(), model=1, Kmax = 15, phi = numeric(), m = numeric()) UseMethod("Segmentor")

Segmentor.default <-function(data=numeric(), model=1, Kmax = 15, phi = numeric(), m = numeric())
{
  if ((model!=1)&(model!=2)&(model!=3)&(model!=4))
    stop("Choose model=1 (Poisson), 2 (normal), 3 (Negative Binomial) or 4 (Normal-Variance)")
  if (length(data)==0)
    stop("Give me a vector of data to segment")

	n = length(data)
  breaks=matrix(0,nrow=Kmax,ncol=Kmax)
  breaks = as.vector(breaks)
  parameters=matrix(0,nrow=Kmax,ncol=Kmax)
  parameters=as.vector(parameters)
  cost=matrix(0,nrow=Kmax,ncol=n)
  cost=as.vector(cost)
  likelihood=rep(0, Kmax)
  likelihood=as.vector(likelihood)
  if ((model==4) & (length(m)==0))
	m = mean(data)
  if ((model==3) & (length(phi)==0))
  {
    h<-15
		Xcum = cumsum(data)
		X2cum = cumsum(data^2)
		M = (Xcum[h:n] - c(0, Xcum[1:(n-h)])) / h
		S2 = (X2cum[h:n] - c(0, X2cum[1:(n-h)])) / (h-1) - h/(h-1)*M^2
		K = M^2 / (S2-M)
		phi = median(K[!is.na(K)])
		while ((phi<0)&(h<(n/2)))
		{
			h<-2*h
			M = (Xcum[h:n] - c(0, Xcum[1:(n-h)])) / h
			S2 = (X2cum[h:n] - c(0, X2cum[1:(n-h)])) / (h-1) - h/(h-1)*M^2
			K = M^2 / (S2-M)
			phi = median(K[!is.na(K)])   	
		}
  }
  if (model==1)
	Rep<-.C("SegmentPoisson", Size = as.integer(n),KMax = as.integer(Kmax), Data = as.integer(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), PACKAGE="Segmentor3IsBack") else if (model==3)
	Rep<-.C("SegmentBinNeg", Size = as.integer(n),KMax = as.integer(Kmax), theta = as.double(phi), Data = as.integer(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), PACKAGE="Segmentor3IsBack") else if (model==2)
	Rep<-.C("SegmentNormal", Size = as.integer(n),KMax = as.integer(Kmax), Data = as.double(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), PACKAGE="Segmentor3IsBack") else if (model==4)
	Rep<-.C("SegmentVariance", Size = as.integer(n),KMax = as.integer(Kmax), mu = as.double(m), Data = as.double(data), Breakpoints = as.integer(breaks), Parameters = as.double(parameters), Likelihood = as.double(likelihood), Cost = as.double(cost), PACKAGE="Segmentor3IsBack")

  breaks=matrix(Rep$Breakpoints,ncol=Kmax)
  breaks = t(breaks)
  parameters = matrix(Rep$Parameters,ncol=Kmax)
  parameters = t(parameters)
  cost = matrix(Rep$Cost,ncol=Kmax)
  cost = t(cost)
  rownames(breaks)<-c("1 segment", paste(2:Kmax, "segments"))
  rownames(parameters)<-c("1 segment", paste(2:Kmax, "segments"))
  colnames(breaks)<-c(paste(1:Kmax, "th break",sep=""))
  colnames(parameters)<-c(paste(1:Kmax, "th parameter",sep=""))
  likelihood=matrix(Rep$Likelihood,ncol=1)
  rownames(likelihood)<-c("1 segment", paste(2:Kmax, "segments"))
  if (model==1) 
  {
      model.dist="Poisson"
      likelihood=likelihood+sum(lgamma(data+1))
      Segmentor.res=new("Segmentor",model=model.dist,breaks=breaks,parameters=parameters,likelihood=likelihood,Kmax=Kmax,Cost=cost)
  }
  if (model==2) 
  {
		data1<-data[-(1:3)]
		data2<-data[-c(1,2,n)]
		data3<-data[-c(1,n-1,n)]
		data4<-data[-c(n-2,n-1,n)]
		d<-c(0.1942, 0.2809, 0.3832, -0.8582)
		v2<-d[1]*data1+d[2]*data2+d[3]*data3+d[4]*data4
		v2<-v2*v2
		var<-sum(v2)/(n-3)
    likelihood = likelihood /(2*var) + n/2*log(2*pi*var)
    model.dist="Normal"
    Segmentor.res=new("Segmentor",model=model.dist,breaks=breaks,parameters=parameters,likelihood=likelihood,Kmax=Kmax,Cost=cost)
  }
  if (model==3) 
  {
      model.dist="Negative binomial"
      Segmentor.res=new("Segmentor",model=model.dist,breaks=breaks, overdispersion=phi,parameters=parameters,likelihood=likelihood, Kmax=Kmax,Cost=cost)
  }
  if (model==4) 
  {
      likelihood = likelihood + n/2*log(2*pi)
      model.dist="Variance Segmentation"
      Segmentor.res=new("Segmentor",model=model.dist,breaks=breaks,mean=m,parameters=parameters, likelihood=likelihood,Kmax=Kmax,Cost=cost)
  }
  Segmentor.res

}

SelectModel <-function(x,penalty="oracle",seuil=n/log(n),keep=FALSE)
{
	if ((penalty!='BIC') & (penalty!='mBIC') & (penalty!='AIC') & (penalty!='oracle'))
		stop("penalty must be BIC, mBIC, AIC or oracle")
	if (class(x)!="Segmentor")
		stop("x must be an object of class Segmentor returned by the Segmentor function")
	n<-getBreaks(x)[1,1]
	Kmax<-getKmax(x)
	sizenr<-function(k) {	vec<-sum(log(diff(c(1,getBreaks(x)[k,1:k]))))}
	saut<-function(Lv, pen,Kseq,seuil=sqrt(n)/log(n),biggest=TRUE)
	{
		J=-Lv;Kmax=length(J); k=1;kv=c();dv=c();pv=c();dmax=1
		while (k<Kmax) {
				pk=(J[(k+1):Kmax]-J[k])/(pen[k]-pen[(k+1):Kmax])
				pm=max(pk); dm=which.max(pk); dv=c(dv,dm); kv=c(kv,k); pv=c(pv,pm)
				if (dm>dmax){  
				  dmax=dm; kmax=k; pmax=pm  
				  }
				k=k+dm
		 } 
		if (biggest)
		{
			pv=c(pv,0); kv=c(kv,Kmax); dv=diff(kv); dmax=max(dv); rt=max(dv); rt=which(dv==rt)
			pmax=pv[rt[length(rt)]]
			alpha=2*pmax
			km=kv[alpha>=pv]; Kh =Kseq[km[1]] 
			return(c(Kh,alpha))
		} else
		{
			paux<-pv[which(kv<=seuil)]
			alpha<-2*min(paux)
			km=kv[alpha>=pv];	Kh =Kseq[km[1]] 
			return(c(Kh,alpha))	
		}
	}
	if (getModel(x)=="Poisson")
	{
		if(penalty=='mBIC')
			K<-which.min(crit<-getLikelihood(x)+0.5*sapply(1:Kmax,sizenr)+(1:Kmax-0.5)*log(n))
		if (penalty=='BIC')
			K<-which.min(crit<-getLikelihood(x)+1:K*log(n))
		if (penalty=='AIC')
			K<-which.min(crit<-getLikelihood(x)+1:K*2)
		if (penalty=='oracle')
		{
			Kseq=1:Kmax
			pen=Kseq*(1+4*sqrt(1.1+log(n/Kseq)))*(1+4*sqrt(1.1+log(n/Kseq)))
			r1=saut(-getLikelihood(x)[Kseq],pen,Kseq)
			crit1<-getLikelihood(x)[Kseq]+r1[2]*pen
			r2=saut(-getLikelihood(x)[Kseq],pen,Kseq,seuil,biggest=FALSE)
			crit2<-getLikelihood(x)[Kseq]+r2[2]*pen
			crit<-cbind(crit1,crit2)
			K<-c(r1[1],r2[1])
		}	
	} else if (getModel(x)=="Negative binomial")
	{
		if(penalty=='mBIC')
			stop("no mBIC for Negative Binomial model")
		if (penalty=='BIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*log(n))
		if (penalty=='AIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*2)	
		if (penalty=='oracle')
		{
			Kseq=1:Kmax
			pen=Kseq*(1+4*sqrt(1.1+log(n/Kseq)))*(1+4*sqrt(1.1+log(n/Kseq)))
			r1=saut(-getLikelihood(x)[Kseq],pen,Kseq)
			crit1<-getLikelihood(x)[Kseq]+r1[2]*pen
			r2=saut(-getLikelihood(x)[Kseq],pen,Kseq,seuil,biggest=FALSE)
			crit2<-getLikelihood(x)[Kseq]+r2[2]*pen
			crit<-cbind(crit1,crit2)
			K<-c(r1[1],r2[1])
		}	
	} else if (getModel(x)=='Normal')
	{
		if(penalty=='mBIC')
			K<-which.min(getLikelihood(x)+0.5*sapply(1:Kmax,sizenr)+(1:Kmax-0.5)*log(n))
		if (penalty=='BIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*log(n))
		if (penalty=='AIC')
			K<-which.min(getLikelihood(x)+((1:Kmax)+1)*2)		
		if (penalty=='oracle')
		{
			Kseq=1:Kmax
			pen=Kseq*(2*log(n/Kseq)+5)
			r1=saut(-getLikelihood(x)[Kseq],pen,Kseq)
			crit1<-getLikelihood(x)[Kseq]+r1[2]*pen
			r2=saut(-getLikelihood(x)[Kseq],pen,Kseq,seuil,biggest=FALSE)
			crit2<-getLikelihood(x)[Kseq]+r2[2]*pen
			crit<-cbind(crit1,crit2)
			K<-c(r1[1],r2[1])
		}	
	} else
	{
		if(penalty=='mBIC')
			stop("no mBIC for Variance model")
		if (penalty=='BIC')
			K<-which.min(getLikelihood(x)+((1:K)+1)*log(n))
		if (penalty=='AIC')
			K<-which.min(getLikelihood(x)+((1:K)+1)*2)		
		if(penalty=='oracle')
			stop("no oracle penalty for Variance model")
	}
	if (keep)
		res<-list(K=K,criterion=crit)
	else res<-K
	return(res)
}



print.Segmentor <-function(x,...)
{
  cat("\n Model used for the segmentation: \n")
  print(getModel(x))
  
  cat("\n Table of optimal breakpoints: \n")
  print(getBreaks(x))

  cat("\n Table of negative log-likelihood for each optimal segmentation: \n")
  print(getLikelihood(x))

  if (is.element("overdispersion",names(x)))
  {
    cat("\n Value of the overdispersion used for the segmentation: \n")
    print(getOverdispersion(x))
  }  
  if (is.element("mean",names(x)))
  {
    cat("\n Value of the mean used for the segmentation: \n")
    print(getMean(x))
  }  
  if (is.element("parameters",names(x)))
  {
    cat("\n Table of parameters: ")
    if (getModel(x)=="Poisson")
      cat("(mean of the signal in each segment) \n")
    if (getModel(x)=="Normal")
      cat("(mean of the signal in each segment) \n")
    if (getModel(x)=="Negative binomial")
      cat("(success-probability of the signal in each segment) \n")
    if (getModel(x)=="Variance Segmentation")
      cat("(Variance of the signal in each segment) \n")
    print(getParameters(x))
  }
}

