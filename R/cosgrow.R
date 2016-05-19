.Einva3=function(a, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime){1/(a^3*(sqrt(OmegaR*a^(-4) + OmegaM*a^(-3) + OmegaK*a^(-2) + OmegaL*cosgrowRhoDE(z=1/a-1, w0=w0, wprime=wprime, rhoDE=1)))^3)}

cosgrow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.8, fSigma8=FALSE, Dist='Co', Mass='Msun', ref){
  z=as.numeric(z)
  if(!Dist %in% c('Co','Ang','m')){stop('Dist must be one of Co, Ang or m')}
  if(!Mass %in% c('Msun','kg')){stop('Mass must be one of Msun or kg')}
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  temp=function(z, H0, OmegaM, OmegaL, OmegaR, w0, wprime){
    OmegaK=1-OmegaM-OmegaL-OmegaR
    OmegaSum=OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)
    Hz=H0*sqrt(OmegaSum)
    CoVel=299792.458*integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    OmegaRatz=(OmegaR*(1+z)^4)/OmegaSum
    OmegaMatz=(OmegaM*(1+z)^3)/OmegaSum
    OmegaLatz=OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)/OmegaSum
    OmegaKatz=(OmegaK*(1+z)^2)/OmegaSum
    Factor=(5*OmegaM/2)*(Hz/H0)*(1+z)*integral(.Einva3, 0, 1/(1+z), OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, OmegaK=OmegaK, w0=w0, wprime=wprime)
    Factor0=(5*OmegaM/2)*integral(.Einva3, 0, 1, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, OmegaK=OmegaK, w0=w0, wprime=wprime)
    Sigma8atz=Sigma8*(Factor/Factor0)/(1+z)
    if(fSigma8==FALSE){
      Rate=-1 - OmegaMatz/2 + OmegaLatz + (5*OmegaMatz)/(2*Factor)
    }else{
      Rate=Sigma8atz*(-1 - OmegaMatz/2 + OmegaLatz + (5*OmegaMatz)/(2*Factor))
    }
    G=6.67384e-11 # m^3 kg^-1 s^-2
    km2m=1000
    Mpc2m=3.08567758e22
    Msol2kg=1.9891e30 # kg
    RhoCrit=(3*Hz^2)/(8*pi*G)*(km2m^2)/Mpc2m^2 #kg/m^3
    if(Mass=='Msun'){RhoCrit=RhoCrit/Msol2kg}
    if(Dist=='Ang'){RhoCrit=RhoCrit*Mpc2m^3}
    if(Dist=='Co'){RhoCrit=RhoCrit*Mpc2m^3/(1+z)^3}
    RhoMean=RhoCrit*OmegaMatz
    Decelq=OmegaMatz/2+OmegaRatz-OmegaLatz
  return=c(z=z, a=1/(1+z), H=Hz, CoVel=CoVel, OmegaM=OmegaMatz, OmegaL=OmegaLatz, OmegaR=OmegaRatz, OmegaK=OmegaKatz, Decelq=Decelq, Factor=Factor, Rate=Rate, Sigma8=Sigma8atz, RhoCrit=RhoCrit, RhoMean=RhoMean)
  }
  return(as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime))))
}

cosgrowz=function(z = 1){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(z)
}

cosgrowa=function(z = 1){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(1/(1 + z))
}

cosgrowH=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  return(H0*sqrt(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)))
}

cosgrowOmegaR=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  return((OmegaR*(1+z)^4)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)))
}

cosgrowOmegaM=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  return((OmegaM*(1+z)^3)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)))
}

cosgrowOmegaL=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  return(OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)))
}

cosgrowOmegaK=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  return((OmegaK*(1+z)^2)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)))
}

cosgrowDecelq=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  OmegaNorm=OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)
  OmegaMatz=(OmegaM*(1+z)^3)/OmegaNorm
  OmegaRatz=(OmegaR*(1+z)^4)/OmegaNorm
  OmegaLatz=OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)/OmegaNorm
  return(OmegaMatz/2+OmegaRatz-OmegaLatz)
}

cosgrowFactor=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  temp=function(z, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime){
    growthfactor=(5*OmegaM/2)*cosgrowH(z, H0=1, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)*(1+z)*integral(.Einva3,0,1/(1+z), OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, OmegaK=OmegaK, w0=w0, wprime=wprime)
    return=growthfactor
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR=OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosgrowRate=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.8, fSigma8=FALSE, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  OmegaLatz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  growthfacttemp=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  if(fSigma8==FALSE){
    Sigma8temp=1
  }else{
    Sigma8temp=cosgrowSigma8(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, Sigma8=Sigma8)
  }
  return(Sigma8temp*(-1 - OmegaMatz/2 + OmegaLatz + (5*OmegaMatz)/(2*growthfacttemp)))
}

cosgrowSigma8=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.8, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  OmegaLatz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  growthfacttempAt0=cosgrowFactor(z=0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  growthfacttempatz=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  return(Sigma8*(growthfacttempatz/growthfacttempAt0)/(1+z))
}

cosgrowFactorApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  OmegaLatz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  return((5*OmegaMatz/2)/(OmegaMatz^(4/7)-OmegaLatz+(1+0.5*OmegaMatz)*(1+OmegaLatz/70)))
}

cosgrowRateApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.8, fSigma8=FALSE, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  if(fSigma8==FALSE){
    Sigma8temp=1
  }else{
    Sigma8temp=cosgrowSigma8(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, Sigma8=Sigma8)
  }
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  OmegaLatz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  return(Sigma8temp*(OmegaMatz^(4/7)+(1+OmegaMatz/2)*(OmegaLatz/70)))
}

cosgrowSigma8Approx=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.8, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  OmegaLatz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  growthfacttempAt0=cosgrowFactorApprox(z=0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  growthfacttempatz=cosgrowFactorApprox(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  return(Sigma8*(growthfacttempatz/growthfacttempAt0)/(1+z))
}

cosgrowRhoCrit=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Dist='Co', Mass='Msun', ref){
  z=as.numeric(z)
  if(!Dist %in% c('Co','Ang','m')){stop('Dist must be one of Co, Ang or m')}
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  G=6.67384e-11 # m^3 kg^-1 s^-2
  km2m=1000
  Mpc2m=3.08567758e22
  Msol2kg=1.9891e30 # kg
  Hub2=cosgrowH(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)^2 # (km/s / Mpc)^2
  RhoCrit=(3*Hub2)/(8*pi*G)*(km2m^2)/Mpc2m^2 #MsolperMpc3
  if(Mass=='Msun'){RhoCrit=RhoCrit/Msol2kg}
  if(Dist=='Ang'){RhoCrit=RhoCrit*Mpc2m^3}
  if(Dist=='Co'){RhoCrit=RhoCrit*Mpc2m^3/(1+z)^3}
  return(RhoCrit)
}

cosgrowRhoMean=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Dist='Co', Mass='Msun', ref){
  z=as.numeric(z)
  if(!Dist %in% c('Co','Ang','m')){stop('Dist must be one of Co, Ang or m')}
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  G=6.67384e-11 # m^3 kg^-1 s^-2
  km2m=1000
  Mpc2m=3.08567758e22
  Msol2kg=1.9891e30 # kg
  Hub2=cosgrowH(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR)^2 # (km/s / Mpc)^2
  RhoCrit=(3*Hub2)/(8*pi*G)*(km2m^2)/Mpc2m^2 #MsolperMpc3
  if(Mass=='Msun'){RhoCrit=RhoCrit/Msol2kg}
  if(Dist=='Ang'){RhoCrit=RhoCrit*Mpc2m^3}
  if(Dist=='Co'){RhoCrit=RhoCrit*Mpc2m^3/(1+z)^3}
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
  RhoMean=RhoCrit*OmegaMatz
  return(RhoMean)
}

cosgrowEoSwDE=function(z=1, w0=-1, wprime=0){w0+2*wprime*(1-1/(1+z))}

cosgrowRhoDE=function(z=1, w0=-1, wprime=0, rhoDE=1){
  rhoDE*((1/(1+z))^(-(3+3*w0+6*wprime)))*exp(-6*wprime*(1-1/(1+z)))
}

cosgrowDeltaVir=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref){
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  if(OmegaM+OmegaL+OmegaR!=1){
    OmegaL=1-OmegaM
    OmegaR=0
    print(paste('Forcing Cosmology to OmegaM=',OmegaM,', OmegaL=',OmegaL,', OmegaR=0 for Bryan & Norman 1998 relation!',sep=''))
  }
  OmegaMatz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR)
  return(18*pi^2+82*(OmegaMatz-1)-39*(OmegaMatz-1)^2)
}

cosgrowCoVel=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    CoVel = 299792.458*integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return=CoVel
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosgrowPecVel=function(z=1, zob=1){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!all(is.finite(zob))){stop('All zob must be finite and numeric')}
  if(!all(zob> -1)){stop('All zob must be > -1')}
  zpec=(1+zob)/(1+z)-1
  PecVel=299792.458*((1+zpec)^2-1)/((1+zpec)^2+1)
  return(PecVel)
}

