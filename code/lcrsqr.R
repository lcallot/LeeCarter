# A function to get the R2 for a given smpl
lcrsqr <- function(smpl){
	
  # estimating the 3 lee carter models.
  LC	 <- est.lc(smpl)

  det.res	<- LC$DLC$lmortdm # demeaned detrended.
  dem.res	<- LC$CLC$lmortdm # demeaned

  # Saving residuals from each model
  clc.res <- LC$CLC$res
  dlc.res <- LC$DLC$res
  tlc.res <- det.res

  rsqr=c(1-sum(clc.res^2)/sum(dem.res^2)
            ,1-sum(dlc.res^2)/sum(dem.res^2)
            ,1-sum(tlc.res^2)/sum(dem.res^2)
            ,1-sum(clc.res^2)/sum(det.res^2)
            ,1-sum(dlc.res^2)/sum(det.res^2)
  		  )
return(rsqr)
}



est.lc <- function(smpl){
	
	LC <- list()
	LC$CLC	<- leecarter(smpl,trend=FALSE)
	LC$DLC	<- leecarter(smpl,trend=TRUE,ktols=TRUE)
	LC$DET	<- leecarter(smpl,trend=TRUE,ktols=TRUE)

	return(LC)
}

