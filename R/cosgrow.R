cosgrow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  OmegaSum=OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL
  Hz=H0*sqrt(OmegaSum)
  OmegaMAtz=(OmegaM*(1+z)^3)/OmegaSum
  OmegaLAtz=OmegaL/OmegaSum
  OmegaKAtz=(OmegaK*(1+z)^2)/OmegaSum
  Factor=(5*OmegaMAtz/2)/(OmegaMAtz^(4/7)-OmegaLAtz+(1+0.5*OmegaMAtz)*(1+OmegaLAtz/70))
  Rate=OmegaMAtz^(4/7)+(1+OmegaMAtz/2)*(OmegaLAtz/70)
  G=6.67384e-11 # m^3 kg^-1 s^-2
  Hub2=H0*sqrt(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL)
  km2m=1000
  Mpc2m=3.08567758e22
  Msol2kg=1.9891e30 # kg
  RhoCrit=(3*Hub2)/(8*pi*G)*(km2m^2)*Mpc2m/Msol2kg #MsolperMpc3
  return(data.frame(z=z, a=1/(1+z), H=Hz, OmegaM=OmegaMAtz, OmegaL=OmegaLAtz, OmegaK=OmegaKAtz, Factor=Factor, Rate=Rate, RhoCrit=RhoCrit))
}

cosgrowz=function(z = 1){
  return(z)
}

cosgrowa=function(z = 1){
  return(1/(1 + z))
}

cosgrowH=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  return(H0*sqrt(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaM=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  return((OmegaM*(1+z)^3)/(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaL=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  return(OmegaL/(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaK=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  return((OmegaK*(1+z)^2)/(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowFactor=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  Einva3=function(a, OmegaM, OmegaL, OmegaK){1/(a^3*(sqrt(OmegaM*a^(-3) + OmegaK*a^(-2) + OmegaL))^3)}
  temp=function(z, OmegaM, OmegaL, OmegaK){
    growthfactor=(5*OmegaM/2)*cosgrowH(z,H0=1,OmegaM=OmegaM,OmegaL=OmegaL)*(1+z)*integrate(Einva3,0,1/(1+z),OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    return=growthfactor
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosgrowRate=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttemp=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(-1 - OmegaMtemp/2 + OmegaLtemp + (5*OmegaMtemp)/(2*growthfacttemp))
}

cosgrowFactorApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return((5*OmegaMtemp/2)/(OmegaMtemp^(4/7)-OmegaLtemp+(1+0.5*OmegaMtemp)*(1+OmegaLtemp/70)))
}

cosgrowRateApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaMtemp=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  OmegaLtemp=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(OmegaMtemp^(4/7)+(1+OmegaMtemp/2)*(OmegaLtemp/70))
}

cosgrowRhoCrit=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  G=6.67384e-11 # m^3 kg^-1 s^-2
  Hub2=cosgrowH(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)^2 # (km/s / Mpc)^2
  km2m=1000
  Mpc2m=3.08567758e22
  Msol2kg=1.9891e30 # kg
  #rhocrit_kgperm3=(3*Hub2)/(8*pi*G) * ((km2m^2)/(Mpc2m^2)) #this is correct, should be ~9.2e-27 for H0=70 z=0
  #rhocrit_MsolperMpc3=rhocrit_kgperm3 * (Mpc2m^3)/Msol2kg #this is correct, should be ~2.8e11 for H0=100 z=0
  RhoCrit=(3*Hub2)/(8*pi*G) * (km2m^2)*Mpc2m/Msol2kg #compact form of the above, in units MsolperMpc3
  return(RhoCrit)
}
