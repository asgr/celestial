cosmapfunc=function(cosparamx='CoVol', cosparamy='z', H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange=c(0,20), step='z', res=100){
  
  paramlistx=c('z', 'a', 'CoDist', 'LumDist', 'CoDistTran', 'DistMod', 'CoVol', 'UniAgeAtz','TravelTime','H','OmegaM','OmegaL','OmegaK','Factor','Rate','RhoCrit')
  paramlisty=c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime','H','OmegaM','OmegaL','OmegaK','Factor','Rate','RhoCrit')
  if(! cosparamx %in% paramlistx){stop('cosparamx is not an allowed cosmological parameter, see help options.')}
  if(! cosparamy %in% paramlisty){stop('cosparamy is not an allowed cosmological parameter, see help options.')}
  #age=FALSE
  #if(cosparamx %in% c('UniAgeAtz','TravelTime')){age=TRUE}
  #if(cosparamy %in% c('UniAgeAtz','TravelTime')){age=TRUE}
  if(cosparamx %in% c('z', 'a', 'CoDist', 'LumDist', 'CoDistTran', 'DistMod', 'CoVol', 'UniAgeAtz','TravelTime')){
    pre_x='cosdist'
  }else{
    pre_x='cosgrow'
  }
  if(cosparamy %in% c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime')){
    pre_y='cosdist'
  }else{
    pre_y='cosgrow'
  }
  if(step=='z'){
    zvals=seq(zrange[1],zrange[2],len=res)
  }
  if(step=='logz'){
    zrangelog=log10(1+zrange)
    zvalslog=seq(zrangelog[1],zrangelog[2],len=res)
    zvals=10^zvalslog-1
  }
  if(step=='a'){
    avals=seq(1/(1+zrange[1]),1/(1+zrange[2]),len=res)
    zvals=1/avals-1
  }
  if(cosparamx %in% c('CoDist', 'LumDist', 'CoDistTran', 'DistMod', 'CoVol', 'UniAgeAtz','TravelTime','H','RhoCrit')){
    combxparams=list(z=zvals, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL)
  }
  if(cosparamx %in% c('z', 'a')){
    combxparams=list(z=zvals)
  }
  if(cosparamx %in% c('OmegaM', 'OmegaL', 'OmegaK','Factor', 'Rate')){
    combxparams=list(z=zvals, OmegaM = OmegaM, OmegaL = OmegaL)
  }
  
  if(cosparamy %in% c('CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime','H','RhoCrit')){
    combyparams=list(z=zvals, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL)
  }
  if(cosparamy %in% c('z', 'a')){
    combyparams=list(z=zvals)
  }
  if(cosparamy %in% c('OmegaM', 'OmegaL', 'OmegaK','Factor', 'Rate')){
    combyparams=list(z=zvals, OmegaM = OmegaM, OmegaL = OmegaL)
  }
  
  tempx=do.call(paste(pre_x,cosparamx,sep=''),combxparams)
  tempy=do.call(paste(pre_y,cosparamy,sep=''),combyparams)
  return=approxfun(tempx,tempy)
}

cosmapval=function(val=50, cosparam='CoVol', H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange=c(0,100), res=10, iter=12, out='cos'){
  temp=function(val, cosparam, H0, OmegaM, OmegaL, zlo, zhi, res, iter, out){
    if(cosparam=='DistMod' & zlo==0){zlo=1e-5}
    zrangetemp=c(zlo, zhi)
    for(i in 1:iter){
    tempz=cosmapfunc(cosparamx=cosparam, cosparamy='z', H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, zrange=zrangetemp, step='z', res=res)
    currentz=tempz(val)
    zlonew=max(zrangetemp[1],currentz-(zrangetemp[2]-zrangetemp[1])/res)
    zhinew=currentz+(zrangetemp[2]-zrangetemp[1])/res
    zrangetemp=c(zlonew,zhinew)
    }
    if(out=='cos'){
      outdist=unlist(cosdist(currentz, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, age = TRUE, error=T))
      outgrow=unlist(cosgrow(currentz, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL))
      output=c(outdist,outgrow[3:9])
      Error=abs(val-output[[cosparam]])
      if(Error>0){Error=Error/output[[cosparam]]}
      output=c(output,MapError=Error)
    }
    if(out=='z'){
      output=currentz
    }
    return=output
  }
  if(out=='cos'){
    output=as.data.frame(t(Vectorize(temp)(val = val, cosparam = cosparam, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, zlo = zrange[1], zhi = zrange[2], res = res, iter = iter, out = out)))
  }
  if(out=='z'){
    output=Vectorize(temp)(val = val, cosparam = cosparam, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, zlo = zrange[1], zhi = zrange[2], res = res, iter = iter, out = out)
  }
  return(output)
}