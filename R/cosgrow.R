cosgrow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Sigma8=0.8, fSigma8=FALSE, Dist='Co', Mass='Msun', ref){
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
  temp=function(z, H0, OmegaM, OmegaL, OmegaR){
    #For parameters that can't use OmegaR
    OmegaK=1-OmegaM-OmegaL
    OmegaSum=OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL
    Hz=H0*sqrt(OmegaSum)
    OmegaMAtz=(OmegaM*(1+z)^3)/OmegaSum
    OmegaLAtz=OmegaL/OmegaSum
    OmegaKAtz=(OmegaK*(1+z)^2)/OmegaSum
    Einva3=function(a, OmegaM, OmegaL, OmegaK){1/(a^3*(sqrt(OmegaM*a^(-3) + OmegaK*a^(-2) + OmegaL))^3)}
    Factor=(5*OmegaM/2)*(Hz/H0)*(1+z)*integrate(Einva3,0,1/(1+z),OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    Factor0=(5*OmegaM/2)*integrate(Einva3,0,1,OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    Sigma8Atz=Sigma8*(Factor/Factor0)/(1+z)
    if(fSigma8==FALSE){
      Rate=-1 - OmegaMAtz/2 + OmegaLAtz + (5*OmegaMAtz)/(2*Factor)
    }else{
      Rate=Sigma8Atz*(-1 - OmegaMAtz/2 + OmegaLAtz + (5*OmegaMAtz)/(2*Factor))
    }
    #For parameters that can use OmegaR
    OmegaK=1-OmegaM-OmegaL-OmegaR
    OmegaSum=OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL
    Hz=H0*sqrt(OmegaSum)
    OmegaMAtz=(OmegaM*(1+z)^3)/OmegaSum
    OmegaLAtz=OmegaL/OmegaSum
    OmegaRAtz=(OmegaR*(1+z)^4)/OmegaSum
    OmegaKAtz=(OmegaK*(1+z)^2)/OmegaSum
    G=6.67384e-11 # m^3 kg^-1 s^-2
    km2m=1000
    Mpc2m=3.08567758e22
    Msol2kg=1.9891e30 # kg
    RhoCrit=(3*Hz^2)/(8*pi*G)*(km2m^2)/Mpc2m^2 #kg/m^3
    if(Mass=='Msun'){RhoCrit=RhoCrit/Msol2kg}
    if(Dist=='Ang'){RhoCrit=RhoCrit*Mpc2m^3}
    if(Dist=='Co'){RhoCrit=RhoCrit*Mpc2m^3/(1+z)^3}
    RhoMean=RhoCrit*OmegaMAtz
    
  return=c(z=z, a=1/(1+z), H=Hz, OmegaM=OmegaMAtz, OmegaL=OmegaLAtz, OmegaR=OmegaRAtz, OmegaK=OmegaKAtz, Factor=Factor, Rate=Rate, Sigma8=Sigma8Atz, RhoCrit=RhoCrit, RhoMean=RhoMean)
  }
  return(as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR))))
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

cosgrowH=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref){
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
  return(H0*sqrt(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaR=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref){
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
  return((OmegaR*(1+z)^4)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaM=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref){
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
  return((OmegaM*(1+z)^3)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaL=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref){
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
  return(OmegaL/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowOmegaK=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref){
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
  return((OmegaK*(1+z)^2)/(OmegaR*(1+z)^4 + OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL))
}

cosgrowFactor=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaK=1-OmegaM-OmegaL
  Einva3=function(a, OmegaM, OmegaL, OmegaK){1/(a^3*(sqrt(OmegaM*a^(-3) + OmegaK*a^(-2) + OmegaL))^3)}
  temp=function(z, OmegaM, OmegaL, OmegaK){
    growthfactor=(5*OmegaM/2)*cosgrowH(z,H0=1,OmegaM=OmegaM,OmegaL=OmegaL,OmegaR=0)*(1+z)*integrate(Einva3,0,1/(1+z),OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000L)$value
    return=growthfactor
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosgrowRate=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMAtz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  OmegaLAtz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  growthfacttemp=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  if(fSigma8==FALSE){
    Sigma8temp=1
  }else{
    Sigma8temp=cosgrowSigma8(z=z, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8)
  }
  return(Sigma8temp*(-1 - OmegaMAtz/2 + OmegaLAtz + (5*OmegaMAtz)/(2*growthfacttemp)))
}

cosgrowSigma8=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMAtz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  OmegaLAtz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  growthfacttempAt0=cosgrowFactor(z=0, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttempAtz=cosgrowFactor(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(Sigma8*(growthfacttempAtz/growthfacttempAt0)/(1+z))
}

cosgrowFactorApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
  }
  OmegaMAtz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  OmegaLAtz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  return((5*OmegaMAtz/2)/(OmegaMAtz^(4/7)-OmegaLAtz+(1+0.5*OmegaMAtz)*(1+OmegaLAtz/70)))
}

cosgrowRateApprox=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  if(fSigma8==FALSE){
    Sigma8temp=1
  }else{
    Sigma8temp=cosgrowSigma8(z=z, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8)
  }
  OmegaMAtz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  OmegaLAtz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  return(Sigma8temp*(OmegaMAtz^(4/7)+(1+OmegaMAtz/2)*(OmegaLAtz/70)))
}

cosgrowSigma8Approx=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  OmegaMAtz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  OmegaLAtz=cosgrowOmegaL(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=0)
  growthfacttempAt0=cosgrowFactorApprox(z=0, OmegaM=OmegaM, OmegaL=OmegaL)
  growthfacttempAtz=cosgrowFactorApprox(z=z, OmegaM=OmegaM, OmegaL=OmegaL)
  return(Sigma8*(growthfacttempAtz/growthfacttempAt0)/(1+z))
}

cosgrowRhoCrit=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Dist='Co', Mass='Msun', ref){
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
  return(RhoCrit)
}

cosgrowRhoMean=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Dist='Co', Mass='Msun', ref){
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
  OmegaMAtz=cosgrowOmegaM(z=z, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR)
  RhoMean=RhoCrit*OmegaMAtz
  return(RhoMean)
}
