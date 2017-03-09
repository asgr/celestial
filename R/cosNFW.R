cosNFW=function(Rad=0, Rho0=2.412e15, Rs=0.03253){
  return(Rho0/((Rad/Rs)*(1+Rad/Rs)^2))
}

cosNFWmass_c=function(Rho0=2.412e15, Rs=0.03253, c=5, Munit=1, Lunit=1e6){
  Rs=Rs*Lunit/1e6
  Rmax=c*Rs
  return(4*pi*Rho0*Rs^3*(log((Rs+Rmax)/Rs)-Rmax/(Rs+Rmax))/Munit)
}

cosNFWmass_Rmax=function(Rho0=2.412e15, Rs=0.03253, Rmax=0.16265, Munit=1, Lunit=1e6){
  Rs=Rs*Lunit/1e6
  Rmax=Rmax*Lunit/1e6
  return(4*pi*Rho0*Rs^3*(log((Rs+Rmax)/Rs)-Rmax/(Rs+Rmax))/Munit)
}

cosNFWsigma=function(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR,  OmegaR=0, Rho = "crit", DeltaVir = 200, Munit = 1, Lunit = 1e+06, Vunit = 1e3, ref){
  Mvir=coshaloRvirToMvir(Rvir=Rs*c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Rho=Rho, DeltaVir=DeltaVir, Munit=Munit, Lunit=Lunit, Vunit=Vunit, ref=ref)
  Rho0=Mvir/cosNFWmass_c(Rho0=1, Rs=Rs, c=c, Munit=Munit, Lunit=Lunit)
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

cosNFWsigma_mean=function(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR,  OmegaR=0, Rho = "crit", DeltaVir = 200, Munit = 1, Lunit = 1e+06, Vunit = 1e3, ref){
  Mvir=coshaloRvirToMvir(Rvir=Rs*c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Rho=Rho, DeltaVir=DeltaVir, Munit=Munit, Lunit=Lunit, Vunit=Vunit, ref=ref)
  Rho0=Mvir/cosNFWmass_c(Rho0=1, Rs=Rs, c=c, Munit=Munit, Lunit=Lunit)
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

cosNFWgamma=function(Rad=0.03253, Rs=0.03253, c=5, SigmaC=1, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR,  OmegaR=0, Rho = "crit", DeltaVir = 200, Munit = 1, Lunit = 1e+06, Vunit = 1e3, ref){
  temp=(cosNFWsigma_mean(Rad=Rad, Rs=Rs, c=c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL,  OmegaR=OmegaR, Rho=Rho, DeltaVir=DeltaVir, Munit=Munit, Lunit=Lunit, Vunit=Vunit, ref=ref)-cosNFWsigma(Rad=Rad, Rs=Rs, c=c, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL,  OmegaR=OmegaR, Rho=Rho, DeltaVir=DeltaVir, Munit=Munit, Lunit=Lunit, Vunit=Vunit, ref=ref))/SigmaC
  return(temp)
}

cosNFWduffym2c=function(M=2e12, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR,  OmegaR=0, Rho = "crit", A=6.71, B=-0.091, C=-0.44, Munit = 1, ref){
  M=M*Munit
  if(Rho=='mean'){M=M/cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL,  OmegaR=OmegaR, ref=ref)}
  return(A*(M/2e12/(H0/100))^B*(1+z)^C)
}

cosNFWvesc=function(Rad=0.16264, Mvir=1e12, c=5, f=Inf, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1 - 
    OmegaM - OmegaR, OmegaR = 0, Rho = "crit", Dist = "Co", DeltaVir = 200, Munit = 1, Lunit = 1e+06, Vunit = 1e3, ref){
  G = 6.67384e-11
  msol_to_kg = 1.98892e+30
  pc_to_m = 3.08568e+16
  g = G * msol_to_kg/(pc_to_m)
  g = g * Munit/(Lunit * Vunit^2)
  
  R200=coshaloMvirToRvir(Mvir=Mvir, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Rho=Rho, Dist=Dist, DeltaVir=200, Munit = Munit, Lunit=Lunit, Vunit=Vunit, ref=ref)
  vcircR200=sqrt(g*Mvir/R200)
  
  x=Rad/R200
	g=log(1+c)-c/(1+c)
	vesc=rep(0,length(x))
	vesc[x<=f]=sqrt((vcircR200^2/g) * (log(1+c*x[x<=f])/x[x<=f]-c/(1+f*c)))
	vesc[x>f]=vcircR200/sqrt(x[x>f])
	return(vesc*sqrt(2))
}

cosNFWvcirc=function(Rad=0.16264, Mvir=1e12, c=5, f=Inf, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1 - 
    OmegaM - OmegaR, OmegaR = 0, Rho = "crit", Dist = "Co", DeltaVir = 200, Munit = 1, Lunit = 1e+06, Vunit = 1e3, ref){
  G = 6.67384e-11
  msol_to_kg = 1.98892e+30
  pc_to_m = 3.08568e+16
  g = G * msol_to_kg/(pc_to_m)
  g = g * Munit/(Lunit * Vunit^2)
  
  R200=coshaloMvirToRvir(Mvir=Mvir, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Rho=Rho, Dist=Dist, DeltaVir=200, Munit = Munit, Lunit=Lunit, Vunit=Vunit, ref=ref)
  Rmax=coshaloMvirToRvir(Mvir=Mvir, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Rho=Rho, Dist=Dist, DeltaVir=DeltaVir, Munit = Munit, Lunit=Lunit, Vunit=Vunit, ref=ref)
  Rho0=Mvir/cosNFWmass_Rmax(Rho0=1, Rs=R200/c, Rmax=Rmax, Munit=Munit, Lunit=Lunit)
  MassCont=cosNFWmass_Rmax(Rho0=Rho0, Rs=R200/c, Rmax=Rad, Munit=Munit, Lunit=Lunit)
  MassContf=cosNFWmass_Rmax(Rho0=Rho0, Rs=R200/c, Rmax=f*Rmax, Munit=Munit, Lunit=Lunit)
  MassCont[Rad/R200>f]=MassContf
  
  vcirc=sqrt(g*MassCont/Rad)
  vcirc[Rad==0]=0
  return(vcirc)
}
