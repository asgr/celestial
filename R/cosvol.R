cosvol=function(area=60,zmax=1,zmin=0,H0=100,OmegaM=0.3,OmegaL=1-OmegaM,inunit='deg2'){
  if(inunit=='amin2'){area=area/3600}
  if(inunit=='asec2'){area=area/12960000}
  if(inunit=='rad2' | inunit=='sr'){area=area*(180/pi)^2}
  vols=cosdistCoVol(z=c(zmin,zmax), H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL)
  totalvol=as.numeric((vols[2]-vols[1])*area*pi/129600)
  return(totalvol)
}