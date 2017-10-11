#' A function to create fake Length/Age relationship for fish
#'
#' This function allows you create fake data, with added noise.
#' @param n Number of fish captured
#' @param maxA Age of the oldest fish captured (years)
#' @param noiseK  Standard deviation of the population growth rate (mm/years)
#' @param noiseL  Standard deviation of the population maximum size (mm)
#' @param K  Mean of the population growth rate (mm/years)
#' @param Linf  Mean of the population maximum size (mm)
#' @return A dataframe containing the simulated Age(year) and Length(mm) along with a plot presenting Length(mm) as a function of age(year).  
#' @keywords Fish, Von Bertalanffy, simulation
#' @export
#' @examples
#' function(12,32,0,01,10,0.5,1200)

# This function creates fake data for fish age/length relationship, enjoy !
fake_fish_creator=function(n,maxA,noiseK,noiseL,K,Linf){
captured_fish=sample(maxA,n)
K=rnorm(length(captured_fish),mean = K, sd=noiseK)
Linf=rnorm(length(captured_fish),mean = Linf, sd=noiseL)

L=vector()
for(i in 1:length(captured_fish)){
L[i]=Linf[i]*(1-exp(-K[i]*captured_fish[i]))
}
plot(captured_fish,L,'p',las=1,xlab='Age (years)',ylab='Length (cm)',pch=16)
return(data.frame(Age=captured_fish,Length=L))
}

