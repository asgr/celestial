cosNFW=function(Rad=0, Rho0=2.412e15, Rs=0.03253){
  return(Rho0/((Rad/Rs)*(1+Rad/Rs)^2))
}

cosNFWmass_c=function(Rho0=2.412e15, Rs=0.03253, c=5){
  Rmax=c*Rs
  return(4*pi*Rho0*Rs^3*(log((Rs+Rmax)/Rs)-Rmax/(Rs+Rmax)))
}

cosNFWmass_Rmax=function(Rho0=2.412e15, Rs=0.03253, Rmax=0.16265){
  return(4*pi*Rho0*Rs^3*(log((Rs+Rmax)/Rs)-Rmax/(Rs+Rmax)))
}

cosNFWsigma=function(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM, Rho = "crit", DeltaVir = 200, ref){
  Mvir=coshaloRvirToMvir(Rvir=Rs*c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Rho=Rho, DeltaVir=DeltaVir, ref=ref)
  Rho0=Mvir/cosNFWmass_c(Rho0=1, Rs=Rs, c=c)
  x=Rad/Rs
  temp=rep(NA, length(x))
  lowx=which(x<1)
  onex=which(x==1)
  hix=which(x>1)
  temp[lowx]=(2*Rs*Rho0/(x[lowx]^2-1))*(1-(2/sqrt(1-x[lowx]^2))*atanh(sqrt((1-x[lowx])/(1+x[lowx]))))
  temp[onex]=(2*Rs*Rho0)/3
  temp[hix]=(2*Rs*Rho0/(x[hix]^2-1))*(1-(2/sqrt(x[hix]^2-1))*atan(sqrt((x[hix]-1)/(1+x[hix]))))
  return(temp)
}

cosNFWsigma_mean=function(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM, Rho = "crit", DeltaVir = 200, ref){
  Mvir=coshaloRvirToMvir(Rvir=Rs*c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Rho=Rho, DeltaVir=DeltaVir, ref=ref)
  Rho0=Mvir/cosNFWmass_c(Rho0=1, Rs=Rs, c=c)
  x=Rad/Rs
  temp=rep(NA, length(x))
  lowx=which(x<1)
  onex=which(x==1)
  hix=which(x>1)
  temp[lowx]=(4/x[lowx]^2)*Rs*Rho0*((2/sqrt(1-x[lowx]^2))*atanh(sqrt((1-x[lowx])/(1+x[lowx])))+log(x[lowx]/2))
  temp[onex]=4*Rs*Rho0*(1+log(0.5))
  temp[hix]=(4/x[hix]^2)*Rs*Rho0*((2/sqrt(x[hix]^2-1))*atan(sqrt((x[hix]-1)/(1+x[hix])))+log(x[hix]/2))
  return(temp)
}

cosNFWgamma=function(Rad=0.03253, Rs=0.03253, c=5, SigmaC=1, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM, Rho = "crit", DeltaVir = 200, ref){
  temp=(cosNFWsigma_mean(Rad=Rad, Rs=Rs, c=c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Rho=Rho, DeltaVir=DeltaVir)-cosNFWsigma(Rad=Rad, Rs=Rs, c=c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Rho=Rho, DeltaVir=DeltaVir))/SigmaC
  return(temp)
}
