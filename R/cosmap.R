cosmapfunc=function(cosparamx='CoVol', cosparamy='z', H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange=c(0,20), step='z', res=100){
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

cosmapval=function(val=50, cosparam='CoVol', H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange=c(0,100), res=10, iter=12, age=FALSE){
  temp=function(val, cosparam, H0, OmegaM, OmegaL, zlo, zhi, res, iter ,age){
    if(cosparam=='DistMod' & zlo==0){zlo=1e-5}
    zrangetemp=c(zlo, zhi)
    for(i in 1:iter){
    tempz=cosmapfunc(cosparamx=cosparam, cosparamy='z', H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, zrange=zrangetemp, step='z', res=res)
    currentz=tempz(val)
    zlonew=max(zrangetemp[1],currentz-(zrangetemp[2]-zrangetemp[1])/res)
    zhinew=currentz+(zrangetemp[2]-zrangetemp[1])/res
    zrangetemp=c(zlonew,zhinew)
    }
    out=cosdist(currentz, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, age = age)
    if(age==FALSE & cosparam %in% c('UniAgeAtz','TravelTime')){
      error=NA
    }else{
      error=abs(val-out[1,cosparam])
      if(error>0){error=error/out[1,cosparam]}
    }
    out=as.numeric(out)
    if(age){outfin=c(z=out[1], a=out[2], CoDist=out[3], LumDist=out[4], AngDist=out[5], CoDistTran=out[6], DistMod=out[7], AngSize=out[8], CoVol=out[9],HubTime=out[10], UniAgeNow=out[11], UniAgeAtz=out[12], TravelTime=out[13], error = error)}
    if(age==FALSE){outfin=c(z=out[1], a=out[2], CoDist=out[3], LumDist=out[4], AngDist=out[5], CoDistTran=out[6],  DistMod=out[7], AngSize=out[8], CoVol=out[9],error = error)}
    return=outfin
  }
  return(as.data.frame(t(Vectorize(temp)(val = val, cosparam = cosparam, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, zlo = zrange[1], zhi = zrange[2], res = res, iter = iter, age = age))))
}