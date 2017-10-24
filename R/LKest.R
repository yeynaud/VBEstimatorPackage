#' Linf and K estimator for a Von Bertalanffy Age/Length relationship
#'
#' This function allows you to estimates the growth rate \code{K} and the theoretical maximum size \code{Linf} from a dataset made of \code{Age} and \code{length} measurements using the method described in Kimura D.K.,(1980) Likelihood methods for the von Bertalanffy growth curve. Fishery Bulletin 77: 765–776.

#' @param Length a vector comtaining the length measurements
#' @param Age  a vector comtaining the age measurements
#' @param plotVBGF  if \code{TRUE} plot of the VB model, \code{FALSE} does not provide a plot (default to \code{TRUE})
#' @param t0_estim  the value of t0 for the model (default to 0)
#' @param na.rm if \code{TRUE} remove NA in the first and second arguments, \code{FALSE} else (default to \code{TRUE})
#' @param boot if \code{TRUE} plot 1000 parameters estimations value issued through a boostrap procedure (the corresponding 95 percent confidence interval is plotted as a dotted line) - default to \code{FALSE}
#' @param boundary if \code{TRUE}, value of K bounded between 0 and 1 - default to \code{FALSE}
#' @return A list containing (i) the Von Bertalanffy growth curve fitted values, (ii) the Von Bertalanffy growth curve parameters , (iii) the confidence interval fitted values obtained using the appraoch described in Kimura (1980), (iv) the estimated Von Bertalanffy parameters using the approach described in Kimura (1980), and (v) the r-squared value.
#' @keywords Fish, Von Bertalanffy, Linf, K, estination, Kimura, boostrap
#' @export
#' @examples
#' no example so far, chill


LKest<-function(Length,Age,plotVBGF=TRUE,t0_estim=0,na.rm=TRUE,boot=FALSE,boundary=FALSE){


if(na.rm){
  bueno=!is.na(Length)&!is.na(Age)
  Length=Length[bueno]
  Age=Age[bueno]

}

t0=t0_estim


vbgf<-function(Linf, k, t0, Age) Linf*(1- exp(-k*(Age-t0)))

vbgf_ls<-function(x) {

  Linf <- x[1]
  k <- x[2]

  temp=(Linf*(1- exp(-k*(Age-t0))))
  return(sum((Length-temp)^2, na.rm=TRUE))

}


temp_output<-optim(c(max(Length),runif(1)),vbgf_ls,method="Nelder-Mead")
if(boundary){
temp_output <-optim(c(temp_output$par[1],temp_output$par[2]),vbgf_ls,method="BFGS")
}
if(boundary==FALSE){
temp_output <-optim(c(temp_output$par[1],temp_output$par[2]),vbgf_ls,method="L-BFGS-B",upper=c(10*max(Length),1),lower=c(0,0))
}

best_VBGF=vector()
best_VBGF[1] = temp_output$par[2]
best_VBGF[2] = temp_output$par[1]


reg=cbind(seq(0,max(Age,na.rm=T),1/500),vbgf(best_VBGF[2],best_VBGF[1], t0, seq(0,max(Age,na.rm=T),1/500)))
reg_exact=cbind(Age,vbgf(best_VBGF[2],best_VBGF[1], 0,Age))
SSres=sum((Length-reg_exact[,2])^2)
SStot=sum((Length-mean(Length))^2)






# The Lima bean
##############################################################################################
N=length(Age)

S_hats=sum((Length - (best_VBGF[2] * (1 - exp(-best_VBGF[1] * (Age - t0)))))^2)


F_temp = qf(0.95,3,N-3)

c_q = S_hats * (1 + ((3/(N-3)) * F_temp))

K_temp=best_VBGF[1]

# Upper left
lima_bean_ul=matrix(NA,1,2)
repeat{

       A= sum((1 - exp(-K_temp * (Age - t0)))^2)
       B=-2*sum(Length * (1 - exp(-K_temp * (Age - t0))))
       C=sum(Length^2)-c_q

      delta=B^2-4*A*C
      if(delta<0|K_temp>10*best_VBGF[1]|K_temp<0){
        break
      }
      Linf_temp =(-B+sqrt(delta))/(2*A)
      lima_bean_ul = rbind(lima_bean_ul,c(K_temp,Linf_temp))
      K_temp=K_temp-0.001
      #Lima_l=
       }
K_temp=best_VBGF[1]
lima_bean_ul=lima_bean_ul[-1,]

# Lower left
lima_bean_ll=matrix(NA,1,2)
repeat{

  A= sum((1 - exp(-K_temp * (Age - t0)))^2)
  B=-2*sum(Length * (1 - exp(-K_temp * (Age - t0))))
  C=sum(Length^2)-c_q

  delta=B^2-4*A*C
  if(delta<0|K_temp>10*best_VBGF[1]|K_temp<0){
    break
  }
  Linf_temp =(-B-sqrt(delta))/(2*A)
  lima_bean_ll = rbind(lima_bean_ll,c(K_temp,Linf_temp))
  K_temp=K_temp-0.001
  #Lima_l=
}
K_temp=best_VBGF[1]
lima_bean_ll=lima_bean_ll[-1,]

# Lower right
lima_bean_lr=matrix(NA,1,2)
repeat{

  A= sum((1 - exp(-K_temp * (Age - t0)))^2)
  B=-2*sum(Length * (1 - exp(-K_temp * (Age - t0))))
  C=sum(Length^2)-c_q

  delta=B^2-4*A*C
  if(delta<0|K_temp>10*best_VBGF[1]|K_temp<0){
    break
  }
  Linf_temp =(-B-sqrt(delta))/(2*A)
  lima_bean_lr = rbind(lima_bean_lr,c(K_temp,Linf_temp))
  K_temp=K_temp+0.001
  #Lima_l=
}
K_temp=best_VBGF[1]
lima_bean_lr=lima_bean_lr[-1,]
# Upper right
lima_bean_ur=matrix(NA,1,2)
repeat{

  A= sum((1 - exp(-K_temp * (Age - t0)))^2)
  B=-2*sum(Length * (1 - exp(-K_temp * (Age - t0))))
  C=sum(Length^2)-c_q

  delta=B^2-4*A*C
  if(delta<0|K_temp>10*best_VBGF[1]|K_temp<0){
    break
  }
  Linf_temp =(-B+sqrt(delta))/(2*A)
  lima_bean_ur = rbind(lima_bean_ur,c(K_temp,Linf_temp))
  K_temp=K_temp+0.001
  #Lima_l=
}
lima_bean_ur=lima_bean_ur[-1,]

lima_bean=rbind(lima_bean_ur,lima_bean_lr[nrow(lima_bean_lr):1,],lima_bean_ll,lima_bean_ul[nrow(lima_bean_ul):1,])

forplot=lima_bean

Length_ori=Length
Age_ori=Age
######

if(boot){
estim_boot=matrix(NA,5000,2)
pb <- txtProgressBar(min = 0, max = 5000, style = 3)
for(boot in 1:5000){
  setTxtProgressBar(pb, boot)
  mix=sample(length(Age),length(Age),replace=T)
  Age=Age_ori[mix]
   Length=Length_ori[mix]
  #temp_output <- optim(c(5,0.5),vbgf_ls,method="Nelder-Mead")
  #temp_output <- optim(c(temp_output$par[1],temp_output$par[2]),vbgf_ls,method="L-BFGS-B",lower=c(0,0),upper=c(1000,10.0))
  temp_output<-optim(c(best_VBGF[2],best_VBGF[1]),vbgf_ls,method="Nelder-Mead")
  temp_output <-optim(c(temp_output$par[1],temp_output$par[2]),vbgf_ls,method="BFGS")
  #temp_output <- optim(c(temp_output$par[1],temp_output$par[2]),vbgf_ls,method="Nelder-Mead")
  #temp_output <-optim(c(temp_output$par[1],temp_output$par[2]),vbgf_ls,method="L-BFGS-B",lower=c(NA,0),upper=c(NA,1))
  #temp_output <- nls(Length ~ vbgf(Linf, k, t0, Age), start=list(Linf=max(Length), k=0.5))

  estim_boot[boot,]=temp_output$par
  #estim_boot[boot,]=c(summary(temp_output)$coefficients[1],summary(temp_output)$coefficients[2])
}
close(pb)
}
######


full_res=list()

full_res[[1]]=reg[which(reg[,1]>=min(Age_ori)&which(reg[,1]<=max(Age_ori))),]
full_res[[2]]=paste('K = ',best_VBGF[1],', Linf = ',best_VBGF[2],sep='')
full_res[[3]]=forplot
full_res[[4]]=c(best_VBGF[1],best_VBGF[2])
full_res[[5]]=c(1-(SSres/SStot))


names(full_res)=c('Von Bertalanffy growth curve fitted values','Von Bertalanffy growth curve parameters ','Lima beans fitted values','Von Bertalanffy parameters','R2')
if(plotVBGF){
  if(boot){
  par(mfrow=c(1,3))}else{par(mfrow=c(1,2))}

  plot(forplot,lwd=3,type='l',las=1,ylab='L inf',xlab='K',main='Parameter estimation,Kimura(1980)')
  points(best_VBGF[1],best_VBGF[2],lwd=3,pch=3,col='red')

  if(boot){
library(MASS)
boot_estimation_contour=kde2d(estim_boot[,2][sort(estim_boot[,2],index.return=T)$ix],estim_boot[,1][sort(estim_boot[,2],index.return=T)$ix])
plot(estim_boot[,2],estim_boot[,1],lwd=3,col='grey',pch=16,main='Parameter estimation, bootstrap',las=1,ylab='L inf',xlab='K',xlim=c(min(boot_estimation_contour$x,na.rm=T),max(boot_estimation_contour$x,na.rm=T)),ylim=c(min(boot_estimation_contour$y,na.rm=T),max(boot_estimation_contour$y,na.rm=T)))
lines(forplot,lwd=3,type='l',las=1,ylab='L inf',xlab='K',main='Parameter estimation, Kimura(1980)',lty=3,col='red')
contour(boot_estimation_contour,levels=.05,col='blue',add=T,lwd=2,labels='95%',labcex=1)
points(mean(estim_boot[,1]),mean(estim_boot[,2]),lwd=3,pch=3,col='blue')
points(best_VBGF[1],best_VBGF[2],lwd=3,pch=3,col='red')
}

plot(Age_ori,Length_ori,pch=16,col='grey',las=1,main='Growth curve',xlab='Age',ylab='Length',xlim=c(0,max(Age_ori)),ylim=c(0,max(Length_ori)))
lines(reg[which(reg[,1]>=min(Age_ori)&which(reg[,1]<=max(Age_ori))),],lty=1,lwd=3)
}else{cat('No VB plot requested',fill=T)}
par(mfrow=c(1,1))
return(full_res)
}



