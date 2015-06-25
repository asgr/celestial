cosNFW=function(Rad=0, Rho0=2.412e5, Rs=0.03253){
  Rho0=Rho0*1e10
  return(Rho0/((Rad/Rs)*(1+Rad/Rs)^2))
}

cosNFWmass_c=function(Rho0=2.412e5, Rs=0.03253, c=5){
  Rho0=Rho0*1e10
  Rmax=c*Rs
  return(4*pi*Rho0*Rs^3*(log((Rs+Rmax)/Rs)-Rmax/(Rs+Rmax)))
}

cosNFWmass_Rmax=function(Rho0=2.412e5, Rs=0.03253, Rmax=0.1626){
  Rho0=Rho0*1e10
  return(4*pi*Rho0*Rs^3*(log((Rs+Rmax)/Rs)-Rmax/(Rs+Rmax)))
}