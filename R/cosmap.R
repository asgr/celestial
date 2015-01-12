cosmapfunc=function(cosparamx='z', cosparamy='CoDist', H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange=c(0,20), step='z', res=100){
  paramlistx=c('z', 'a', 'CoDist', 'LumDist', 'CoDistTran', 'DistMod', 'CoVol', 'UniAgeAtz','TravelTime')
  paramlisty=c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime')
  if(! cosparamx %in% paramlistx){stop('cosparamx is not an allowed cosmological parameter, see help options.')}
  if(! cosparamy %in% paramlisty){stop('cosparamy is not an allowed cosmological parameter, see help options.')}
  age=FALSE
  if(cosparamx %in% c('UniAgeAtz','TravelTime')){age=TRUE}
  if(cosparamy %in% c('UniAgeAtz','TravelTime')){age=TRUE}
  if(step=='z'){zvals=seq(zrange[1],zrange[2],len=res)}
  if(step=='logz'){
    zrangelog=log10(1+zrange)
    zvalslog=seq(zrangelog[1],zrangelog[2],len=res)
    zvals=10^zvalslog-1
  }
  if(step=='a'){
    avals=seq(1/(1+zrange[1]),1/(1+zrange[2]),len=res)
    zvals=1/avals-1
  }
  temp=cosdist(zvals, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, age=age)[,c(cosparamx,cosparamy)]
  return=approxfun(temp[,1],temp[,2])
}

cosmapval=function(val=1, cosparamx='z', cosparamy='CoDist', H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange=c(0,20), step='z', res=100){
  temp=cosmapfunc(cosparamx=cosparamx, cosparamy=cosparamy, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, zrange=zrange, res=res, step=step)
  return(temp(val))
}