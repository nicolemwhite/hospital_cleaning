#R functions required for fitting FP1 and FP2 - scaling, centering and transformations


compute_scale_x<-function(x){
  lrange.x<-log10(max(x)-min(x))
  scale.x<-10^(sign(lrange.x)*floor(abs(lrange.x)))
  return(scale.x)
}

fp.scale <- function(x, scaling = TRUE)
{
  scale <- 1; shift <- 0
  if(scaling) {
    if(min(x) <= 0) {
      z <- diff(sort(x))
      shift <- min(z[z > 0]) - min(x)
      shift <- ceiling(shift*10)/10
    } 
    #
    #	range <- median(x+shift)
    range <- mean(x+shift)
    #    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
    scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
  }
  #
  return(list(shift=shift, scale=scale))
}

#transformation of (scaled) continuous covariate for a single p
FP_transform_x_1<-function(x_scaled,p,shift){

  #compute center of x (center)
  #assumes that x has already been scaled
  center<-mean(x_scaled+shift)

  #p=0 case
  if(p==0){
    y<-log(x_scaled+shift) #-log(center)
  }
  #p!=0 case
  if(p!=0){
    y<-((x_scaled+shift)^p) #-(center^p)
  }
#  y<-ifelse(p!=0,x_scaled^p,log(x_scaled))
  
  return(y)
} 

#transformation of (scaled) continuous covariate for repeated powers case (only computes second term)
FP_transform_x_2<-function(x_scaled,p){
  
  center<-mean(x_scaled)
  
  #p=0 case
  if(p==0){
    y<-(log(x_scaled)^2)-(log(center)^2)
  }
  
  #p!=0 case
  if(p!=0){
    y<-(x_scaled^p)*log(x_scaled)-(center^p)*log(center)
  }
  return(y)
}


#function to extract all AICs given list object
extract_AIC_FP<-function(glmer.object){
  out<-sapply(1:length(glmer.object),function(m) extractAIC(glmer.object[[m]])[2])
  return(out)
}