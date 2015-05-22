cosvarcar=function(aside=50, bside=50, cside=50, subsets=1){
  temp=sort(c(aside,bside),decreasing = TRUE)
  aside=temp[1];bside=temp[2]
  x=((aside/bside)-1.0)**0.5
  return(cv=(1.0-0.03*x)*(219.7-52.4*log10(aside*bside*291.0) + 3.21*(log10(aside*bside*291.0))^2)/sqrt(subsets*cside/291.0))
}

cosvarsph=function(long = c(129, 141), lat = c(-2, 3), zmax=1, zmin=0, subsets=1, inunit='deg'){
  if(inunit %in% c('deg','amin','asec','rad','sex')==FALSE){stop('inunit must be one of deg, amin, asec, rad or sex')}
  if(length(long)==1){long=c(0,long)}
  if(length(lat)==1){lat=c(0,lat)}
  fullsky=129600/pi
  if(inunit=='sex'){long=hms2deg(long,sep=sep);lat=dms2deg(lat,sep=sep)}
  if(inunit=='amin'){long=long/60;lat=lat/60}
  if(inunit=='asec'){long=long/3600;lat=lat/3600}
  if(inunit=='rad'){long=long*180/pi;lat=lat*180/pi}
  CoDistLow = cosdistCoDist(z=zmin,H0=70,OmegaM=0.3)     
  CoDistHigh = cosdistCoDist(z=zmax,H0=70,OmegaM=0.3)
  cside=CoDistHigh-CoDistLow
  aside=cos(mean(lat)*pi/180)*(abs(diff(long))/360)*(CoDistLow+cside/2)
  bside=(abs(diff(long))/180)*(CoDistLow+cside/2)
  return(cosvarcar(aside=aside, bside=bside, cside=cside, subsets=subsets))
}