#' fits the model of a Lee and Carter (1992, JASA) on a matrix of log mortality data. 
#' 
#' 	@param logm A list with entries: logm (a matrix of log mortality data), cn (country name, one of: USA, FRA, DNK, NED), gen (gender, one of: female, male, total), startyear (start of the sample), endyear (end of the sample), startage , endage. 
#'   @param trend Should the data be detrended prior to the singular value decomposition. Default = TRUE. 
#'   @param ols  Estimate the loadings on the stochastic trend by OLS. Default = FALSE. 
#'   @return The initial list augmented with the parameters of the Lee Carter model. 
leecarter<-function(lcdat,trend=TRUE,ktols=FALSE){

	logm	<-lcdat$logm
	model	<-lcdat

# Extracting the deterministic loadings	
	if(!trend){
		lm.det	<-lm(logm~1)
		alpha	<-lm.det$coeff[1,]}

	if(trend){
		lt<-1:nrow(logm)#linear trend
		lt<-lt-mean(lt) #demeaning the trend
		lm.det<-lm(logm~1+lt)
		
		alpha<-lm.det$coeff[1,]
		gamma<-lm.det$coeff[2,]
		model$gamma<-gamma}


#getting the demeaned mortality rates
	logmdm<-lm.det$residuals


#Lee Carter kt and beta + identification 
#singular value decomposition (estimation of the LC trend)
	singv	<-svd(logmdm,1,1)
	kt		<- matrix(singv$d[1]*sum(singv$v)*singv$u,ncol=1)
	model$kt<- kt

	# Computing beta from the eigenvalues
	if(!ktols){
		beta <- matrix(singv$v/sum(singv$v),nrow=1)	
		btci <- NULL}
	# Computing beta by ols, and appending std.errs
	if(ktols){
		ols  <- lm(logmdm ~ kt - 1)
		xxi <- solve( t(kt)%*%kt )
		ssq <- diag(t(residuals(ols))%*%residuals(ols))/(length(kt)-1)
		beta <- coef(ols)[1,]
		btci <- matrix(1.96*sqrt(xxi*ssq/length(kt)),nrow=1)
	}
		
	model$beta <- matrix(beta,nrow=1)
	model$btci <- btci

	names(model$beta) <- model$minage:model$maxage
	if(!is.null(model$btci))names(model$btci) <- model$minage:model$maxage
	names(model$kt) <- model$startyear:model$endyear

# Computing the fitted series under the different scenarios for the model.	
	mfit <- matrix(1,nrow=nrow(logm),ncol=1) %*% matrix(alpha,nrow=1)

	 mfit <- mfit + kt %*% beta 
	if(trend) mfit <- mfit + matrix(lt,ncol=1) %*% matrix(gamma,nrow=1)

# computing R2:
	res<-logm - mfit
	RSS<-colSums(res^2)
	TSS<-colSums((logm-matrix(1,nrow(logm),1)%*%t(as.matrix(colMeans(logm))))^2)
	colR2 <- 1-RSS/TSS

# Creating the output list	
	model$colR2<-colR2
	model$fitted<-mfit
	model$res<-res
	model$alpha<-alpha
	model$lmort<-logm
	model$lmortdm<-logmdm
	model$trend<-trend

	names(model$alpha)<-model$minage:model$maxage
	if(trend) names(model$gamma)<-model$minage:model$maxage
	
	return(model)
	}

