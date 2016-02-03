.getcos=function(ref){
  cosref = NULL
  data('cosref',envir = environment())
  allownames=tolower(as.character(cosref[,'Ref']))
  if(tolower(ref) %in% allownames==FALSE){stop(paste('Provided ref name is not allowed, must be one of',paste(as.character(cosref[,'Ref']),sep='',collapse=', '),' (case insensitive). See ?cosref for details.'))}
  out=as.numeric(cosref[allownames==tolower(ref),])
  names(out)=colnames(cosref)
  return(out)
}

.Einv=function(z, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime){
  return(1/sqrt(OmegaR*(1+z)^4 + OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)))
}

.Einvz=function(z, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime){
  return(1/(sqrt(OmegaR*(1+z)^4 + OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1)) * (1 + z)))
}

cosdist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, age=FALSE, ref, error=FALSE){
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
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0=w0, wprime=wprime) {
    HubDist = (299792.458/H0)
    temp = suppressWarnings(integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
    CoDist = HubDist * temp
    if(error){
      if(z>0){
        RelError = 1e-8
      }else{
        RelError=0
      }
    }
    
    if(OmegaK==0){
      CoDistTran = CoDist
      CoVol = ((4/3) * pi * CoDist^3)/1e9
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
        CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asinh(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
      }
    if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
        CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asin(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
      }
    }
      
    a=1/(1+z)
    LumDist = (1+z)*CoDistTran
    AngDist = CoDistTran/(1+z)
    if(z>=0){DistMod = 5*log10(LumDist)+25}else{DistMod=NA}
    AngScale = AngDist*(pi/(180*60*60))*1000
      
    if (age) {
      HT = (3.08568025e+19/(H0*31556926))/1e9
      UniAge = HT*integral(.Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
      zAge = HT*integral(.Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    }
    if(error){
      if (age) {
        return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngScale = AngScale, CoVol = CoVol, HubTime = HT, UniAgeNow = UniAge, UniAgeAtz = UniAge - zAge, TravelTime = zAge, RelError=RelError)
      }
      else {
        return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngScale = AngScale, CoVol = CoVol, RelError=RelError)
      }
    }else{
      if (age) {
        return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngScale = AngScale, CoVol = CoVol, HubTime = HT, UniAgeNow = UniAge, UniAgeAtz = UniAge - zAge, TravelTime = zAge)
      }
      else {
        return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngScale = AngScale, CoVol = CoVol)
      }
    }
  }
  return(as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))))
}

cosdistz=function(z = 1){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(z)
}

cosdista=function(z = 1){
  z=as.numeric(z)
  if(!all(is.finite(z))){stop('All z must be finite and numeric')}
  if(!all(z> -1)){stop('All z must be > -1')}
  return(1/(1 + z))
}

cosdistCoDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return=CoDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistCoDistTran=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    return=CoDistTran
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistLumDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    LumDist = (1+z) * CoDistTran
    return=LumDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistAngDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    AngDist = CoDistTran / (1 + z)
    return=AngDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistDistMod=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    if(z>=0){DistMod = 5*log10(CoDistTran*(1+z))+25}else{DistMod=NA}
    return=DistMod
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistAngScale=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    AngScale = (CoDistTran / (1 + z)) * (pi/(180 * 60 * 60)) * 1000
    return=AngScale
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistAngSize=function(z=1, Size=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Dim=1, Dist='Co', outunit='deg', ref){
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
  if (outunit %in% c("deg", "amin", "asec", "rad") == 
        FALSE) {
        stop("outunit must be one of deg, amin, asec or rad")
    }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  temp = function(z, Size, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime, Dim, Dist) {
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    Size=Size/2
    if(Dist=='Co'){Size=Size/(1+z)}
    if(Dim %in% 1:2){Ang = atan(Size/(CoDistTran / (1 + z)))}
    if(Dim == 3){Ang = asin(Size/(CoDistTran / (1 + z)))}
    Ang=Ang*2
    if(outunit=='deg'){Ang=Ang*180/pi}
    if(outunit=='amin'){Ang=60*Ang*180/pi}
    if(outunit=='asec'){Ang=3600*Ang*180/pi}
    return=Ang
  }
  return(Vectorize(temp)(z = z, Size=Size, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime, Dim=Dim, Dist=Dist))
}

cosdistAngArea=function(z=1, Size=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Dim=2, Dist='Co', outunit='deg2', ref){
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
  if (outunit %in% c("deg2", "amin2", "asec2", "rad2", "sr") == 
        FALSE) {
        stop("outunit must be one of deg2, amin2, asec2, rad2 or sr")
    }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  temp = function(z, Size, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime, Dim, Dist) {
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    Size=Size/2
    if(Dist=='Co'){Size=Size/(1+z)}
    if(Dim %in% 1:2){Ang = atan(Size/(CoDistTran / (1 + z)))}
    if(Dim == 3){Ang = asin(Size/(CoDistTran / (1 + z)))}
    if(outunit=='deg2'){Ang=Ang*180/pi}
    if(outunit=='amin2'){Ang=60*Ang*180/pi}
    if(outunit=='asec2'){Ang=3600*Ang*180/pi}
    Area=pi*Ang^2
    return=Area
  }
  return(Vectorize(temp)(z = z, Size=Size, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime, Dim=Dim, Dist=Dist))
}

cosdistCoVol=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    
    HubDist = (299792.458/H0)
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
      CoVol = ((4/3) * pi * CoDist^3)/1e9
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
        CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asinh(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
        CoVol = ((4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asin(sqrt(abs(OmegaK))*(CoDistTran/HubDist))))/1e9
      }
    }
    return=CoVol
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistUniAgeNow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    UniAge = HT * integral(.Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return=UniAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistUniAgeAtz=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    UniAge = HT * integral(.Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    zAge = HT * integral(.Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return=UniAge-zAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistTravelTime=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    zAge = HT * integral(.Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return=zAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistHubTime=function(H0 = 100){
 return((3.08568025e+19/(H0 * 31556926))/1e9)
}

cosdistRelError=function(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
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
  temp = function(z, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    temp = integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
      if(z>0){
        RelError = 1e-8
      }else{
        RelError=0
      }
    return=RelError
  }
  return(Vectorize(temp)(z = z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime))
}

cosdistAngDist12=function(z1=1,z2=2, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  HubDist = (299792.458/H0)
  z1=as.numeric(z1)
  z2=as.numeric(z2)
  if(!all(is.finite(z1))){stop('All z1 must be finite and numeric')}
  if(!all(is.finite(z2))){stop('All z2 must be finite and numeric')}
  if(!all(z1> -1)){stop('All z1 must be > -1')}
  if(!all(z2> -1)){stop('All z2 must be > -1')}
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  OmegaK=1-OmegaM-OmegaL-OmegaR
  if(OmegaK<0){stop('OmegaK must be >=0 to use this function!')}
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    CoDist = HubDist * integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    if(OmegaK==0){
      CoDistTran = CoDist
    }else{
      if(OmegaK>0){
        CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
      }
      if(OmegaK<0){
        CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
      }
    }
    return=CoDistTran
  }
  CoDistTran1=Vectorize(temp)(z = z1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  CoDistTran2=Vectorize(temp)(z = z2, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  AngDist12=(1/(1+z2))*(CoDistTran2*sqrt(1+OmegaK*CoDistTran1^2/HubDist^2) - CoDistTran1*sqrt(1+OmegaK*CoDistTran2^2/HubDist^2))
  return(AngDist12)
}

cosdistCrit=function(z_lens=1, z_source=2, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, ref){
  #if(any(z_lens>z_source)){stop('All z_lens must be less than z_source!')}
  c=299792458
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g/1e6 #Get into Mpc units
  Dl=cosdistAngDist(z=z_lens, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  Ds=cosdistAngDist(z=z_source, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  Dls=cosdistAngDist12(z1=z_lens, z2=z_source, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  SigmaC=(c^2/(4*pi*g))*(Ds/(Dl*Dls))
  SigmaC[z_lens>=z_source]=0
  return(SigmaC)
}

cosdistCoDist12ang=function(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, inunit='deg', ref){
  if(inunit %in% c('deg','amin','asec','rad')==FALSE){stop('inunit must be one of deg, amin, asec, or rad')}
  if(inunit=='amin'){ang=ang/60}
  if(inunit=='asec'){ang=ang/3600}
  if(inunit=='rad'){ang=ang*180/pi}
  HubDist = (299792.458/H0)
  z1 = as.numeric(z1)
  z2 = as.numeric(z2)
  if (!all(is.finite(z1))) {
    stop("All z1 must be finite and numeric")
  }
  if (!all(is.finite(z2))) {
    stop("All z2 must be finite and numeric")
  }
  if (!all(z1 > -1)) {
    stop("All z1 must be > -1")
  }
  if (!all(z2 > -1)) {
    stop("All z2 must be > -1")
  }
  if (!missing(ref)) {
    params = .getcos(ref)
    H0 = as.numeric(params["H0"])
    OmegaM = as.numeric(params["OmegaM"])
    OmegaL = as.numeric(params["OmegaL"])
    if (!is.na(params["OmegaR"])) {
      OmegaR = as.numeric(params["OmegaR"])
    }
  }
  OmegaK = 1 - OmegaM - OmegaL - OmegaR
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    rDist = integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return = rDist
  }
  rDist1 = Vectorize(temp)(z = z1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  rDist2 = Vectorize(temp)(z = z2, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  if(OmegaK==0){
    #k is 0
    Skr12squared= rDist1^2+
      rDist2^2-
      2*rDist1*rDist2*cos(ang*pi/180)
    CoSep12=HubDist*sqrt(Skr12squared)
  }
  if(OmegaK>0){
    #k is -ve so use sinh and cosh
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sinh(rDist1)^2*cosh(rDist2)^2+
      sinh(rDist2)^2*cosh(rDist1)^2-
      sinh(rDist1)^2*sinh(rDist2)^2*sin(ang*pi/180)^2-
      2*sinh(rDist1)*sinh(rDist2)*cosh(rDist1)*cosh(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asinh(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  if(OmegaK<0){
    #k is +ve so use sin and cos
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sin(rDist1)^2*cos(rDist2)^2+
      sin(rDist2)^2*cos(rDist1)^2+
      sin(rDist1)^2*sin(rDist2)^2*sin(ang*pi/180)^2-
      2*sin(rDist1)*sin(rDist2)*cos(rDist1)*cos(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asin(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  return(CoSep12)
}

cosdistAngDist12ang=function(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, inunit='deg', ref){
  if(inunit %in% c('deg','amin','asec','rad')==FALSE){stop('inunit must be one of deg, amin, asec, or rad')}
  if(inunit=='amin'){ang=ang/60}
  if(inunit=='asec'){ang=ang/3600}
  if(inunit=='rad'){ang=ang*180/pi}
  HubDist = (299792.458/H0)
  z1 = as.numeric(z1)
  z2 = as.numeric(z2)
  if (!all(is.finite(z1))) {
    stop("All z1 must be finite and numeric")
  }
  if (!all(is.finite(z2))) {
    stop("All z2 must be finite and numeric")
  }
  if (!all(z1 > -1)) {
    stop("All z1 must be > -1")
  }
  if (!all(z2 > -1)) {
    stop("All z2 must be > -1")
  }
  if (!missing(ref)) {
    params = .getcos(ref)
    H0 = as.numeric(params["H0"])
    OmegaM = as.numeric(params["OmegaM"])
    OmegaL = as.numeric(params["OmegaL"])
    if (!is.na(params["OmegaR"])) {
      OmegaR = as.numeric(params["OmegaR"])
    }
  }
  OmegaK = 1 - OmegaM - OmegaL - OmegaR
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    rDist = integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return = rDist
  }
  rDist1 = Vectorize(temp)(z = z1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  rDist2 = Vectorize(temp)(z = z2, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  if(OmegaK==0){
    #k is 0
    Skr12squared= rDist1^2+
      rDist2^2-
      2*rDist1*rDist2*cos(ang*pi/180)
    CoSep12=HubDist*sqrt(Skr12squared)
  }
  if(OmegaK>0){
    #k is -ve so use sinh and cosh
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sinh(rDist1)^2*cosh(rDist2)^2+
      sinh(rDist2)^2*cosh(rDist1)^2-
      sinh(rDist1)^2*sinh(rDist2)^2*sin(ang*pi/180)^2-
      2*sinh(rDist1)*sinh(rDist2)*cosh(rDist1)*cosh(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asinh(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  if(OmegaK<0){
    #k is +ve so use sin and cos
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sin(rDist1)^2*cos(rDist2)^2+
      sin(rDist2)^2*cos(rDist1)^2+
      sin(rDist1)^2*sin(rDist2)^2*sin(ang*pi/180)^2-
      2*sin(rDist1)*sin(rDist2)*cos(rDist1)*cos(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asin(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  CoDiff=rDist1*HubDist+CoSep12
  MaxCo=cosdistCoDist(1e5, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime)
  if(CoDiff>MaxCo){
    print('Objects not in common light cone!')
    out=NA
  }else{
    zem=cosmapval(CoDiff, 'CoDist', H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime, zrange = c(0,1e5), out='z')
    zeff=cosdistzeff(zref=z1,zem=zem)
    out=CoSep12/(1+zeff)
  }
  return(out)
}

cosdistLumDist12ang=function(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, inunit='deg', ref){
  if(inunit %in% c('deg','amin','asec','rad')==FALSE){stop('inunit must be one of deg, amin, asec, or rad')}
  if(inunit=='amin'){ang=ang/60}
  if(inunit=='asec'){ang=ang/3600}
  if(inunit=='rad'){ang=ang*180/pi}
  HubDist = (299792.458/H0)
  z1 = as.numeric(z1)
  z2 = as.numeric(z2)
  if (!all(is.finite(z1))) {
    stop("All z1 must be finite and numeric")
  }
  if (!all(is.finite(z2))) {
    stop("All z2 must be finite and numeric")
  }
  if (!all(z1 > -1)) {
    stop("All z1 must be > -1")
  }
  if (!all(z2 > -1)) {
    stop("All z2 must be > -1")
  }
  if (!missing(ref)) {
    params = .getcos(ref)
    H0 = as.numeric(params["H0"])
    OmegaM = as.numeric(params["OmegaM"])
    OmegaL = as.numeric(params["OmegaL"])
    if (!is.na(params["OmegaR"])) {
      OmegaR = as.numeric(params["OmegaR"])
    }
  }
  OmegaK = 1 - OmegaM - OmegaL - OmegaR
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    rDist = integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return = rDist
  }
  rDist1 = Vectorize(temp)(z = z1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  rDist2 = Vectorize(temp)(z = z2, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  if(OmegaK==0){
    #k is 0
    Skr12squared= rDist1^2+
      rDist2^2-
      2*rDist1*rDist2*cos(ang*pi/180)
    CoSep12=HubDist*sqrt(Skr12squared)
  }
  if(OmegaK>0){
    #k is -ve so use sinh and cosh
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sinh(rDist1)^2*cosh(rDist2)^2+
      sinh(rDist2)^2*cosh(rDist1)^2-
      sinh(rDist1)^2*sinh(rDist2)^2*sin(ang*pi/180)^2-
      2*sinh(rDist1)*sinh(rDist2)*cosh(rDist1)*cosh(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asinh(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  if(OmegaK<0){
    #k is +ve so use sin and cos
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sin(rDist1)^2*cos(rDist2)^2+
      sin(rDist2)^2*cos(rDist1)^2+
      sin(rDist1)^2*sin(rDist2)^2*sin(ang*pi/180)^2-
      2*sin(rDist1)*sin(rDist2)*cos(rDist1)*cos(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asin(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  CoDiff=rDist1*HubDist+CoSep12
  MaxCo=cosdistCoDist(1e5, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime)
  if(CoDiff>MaxCo){
    print('Objects not in common light cone!')
    out=NA
  }else{
    zem=cosmapval(CoDiff,'CoDist',H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime, zrange = c(0,1e5), out='z')
    zeff=cosdistzeff(zref=z1,zem=zem)
    out=CoSep12*(1+zeff)
  }
  return(out)
}

cosdistzem12ang=function(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, inunit='deg', ref){
  if(inunit %in% c('deg','amin','asec','rad')==FALSE){stop('inunit must be one of deg, amin, asec, or rad')}
  if(inunit=='amin'){ang=ang/60}
  if(inunit=='asec'){ang=ang/3600}
  if(inunit=='rad'){ang=ang*180/pi}
  HubDist = (299792.458/H0)
  z1 = as.numeric(z1)
  z2 = as.numeric(z2)
  if (!all(is.finite(z1))) {
    stop("All z1 must be finite and numeric")
  }
  if (!all(is.finite(z2))) {
    stop("All z2 must be finite and numeric")
  }
  if (!all(z1 > -1)) {
    stop("All z1 must be > -1")
  }
  if (!all(z2 > -1)) {
    stop("All z2 must be > -1")
  }
  if (!missing(ref)) {
    params = .getcos(ref)
    H0 = as.numeric(params["H0"])
    OmegaM = as.numeric(params["OmegaM"])
    OmegaL = as.numeric(params["OmegaL"])
    if (!is.na(params["OmegaR"])) {
      OmegaR = as.numeric(params["OmegaR"])
    }
  }
  OmegaK = 1 - OmegaM - OmegaL - OmegaR
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    rDist = integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return = rDist
  }
  rDist1 = Vectorize(temp)(z = z1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  rDist2 = Vectorize(temp)(z = z2, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  if(OmegaK==0){
    #k is 0
    Skr12squared= rDist1^2+
      rDist2^2-
      2*rDist1*rDist2*cos(ang*pi/180)
    CoSep12=HubDist*sqrt(Skr12squared)
  }
  if(OmegaK>0){
    #k is -ve so use sinh and cosh
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sinh(rDist1)^2*cosh(rDist2)^2+
      sinh(rDist2)^2*cosh(rDist1)^2-
      sinh(rDist1)^2*sinh(rDist2)^2*sin(ang*pi/180)^2-
      2*sinh(rDist1)*sinh(rDist2)*cosh(rDist1)*cosh(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asinh(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  if(OmegaK<0){
    #k is +ve so use sin and cos
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sin(rDist1)^2*cos(rDist2)^2+
      sin(rDist2)^2*cos(rDist1)^2+
      sin(rDist1)^2*sin(rDist2)^2*sin(ang*pi/180)^2-
      2*sin(rDist1)*sin(rDist2)*cos(rDist1)*cos(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asin(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  CoDiff=rDist1*HubDist+CoSep12
  MaxCo=cosdistCoDist(1e5, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime)
  if(CoDiff>MaxCo){
    print('Objects not in common light cone!')
    out=NA
  }else{
    out=cosmapval(CoDiff,'CoDist',H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime, zrange = c(0,1e5), out='z')
  }
  return(out)
}

cosdistzeff12ang=function(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, inunit='deg', ref){
  if(inunit %in% c('deg','amin','asec','rad')==FALSE){stop('inunit must be one of deg, amin, asec, or rad')}
  if(inunit=='amin'){ang=ang/60}
  if(inunit=='asec'){ang=ang/3600}
  if(inunit=='rad'){ang=ang*180/pi}
  HubDist = (299792.458/H0)
  z1 = as.numeric(z1)
  z2 = as.numeric(z2)
  if (!all(is.finite(z1))) {
    stop("All z1 must be finite and numeric")
  }
  if (!all(is.finite(z2))) {
    stop("All z2 must be finite and numeric")
  }
  if (!all(z1 > -1)) {
    stop("All z1 must be > -1")
  }
  if (!all(z2 > -1)) {
    stop("All z2 must be > -1")
  }
  if (!missing(ref)) {
    params = .getcos(ref)
    H0 = as.numeric(params["H0"])
    OmegaM = as.numeric(params["OmegaM"])
    OmegaL = as.numeric(params["OmegaL"])
    if (!is.na(params["OmegaR"])) {
      OmegaR = as.numeric(params["OmegaR"])
    }
  }
  OmegaK = 1 - OmegaM - OmegaL - OmegaR
  temp = function(z, H0, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime) {
    rDist = integral(.Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
    return = rDist
  }
  rDist1 = Vectorize(temp)(z = z1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  rDist2 = Vectorize(temp)(z = z2, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, OmegaK = OmegaK, w0=w0, wprime=wprime)
  if(OmegaK==0){
    #k is 0
    Skr12squared= rDist1^2+
      rDist2^2-
      2*rDist1*rDist2*cos(ang*pi/180)
    CoSep12=HubDist*sqrt(Skr12squared)
  }
  if(OmegaK>0){
    #k is -ve so use sinh and cosh
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sinh(rDist1)^2*cosh(rDist2)^2+
      sinh(rDist2)^2*cosh(rDist1)^2-
      sinh(rDist1)^2*sinh(rDist2)^2*sin(ang*pi/180)^2-
      2*sinh(rDist1)*sinh(rDist2)*cosh(rDist1)*cosh(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asinh(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  if(OmegaK<0){
    #k is +ve so use sin and cos
    rDist1=rDist1*sqrt(abs(OmegaK))
    rDist2=rDist2*sqrt(abs(OmegaK))
    Skr12squared= sin(rDist1)^2*cos(rDist2)^2+
      sin(rDist2)^2*cos(rDist1)^2+
      sin(rDist1)^2*sin(rDist2)^2*sin(ang*pi/180)^2-
      2*sin(rDist1)*sin(rDist2)*cos(rDist1)*cos(rDist2)*cos(ang*pi/180)
    CoSep12=HubDist*asin(sqrt(Skr12squared))/sqrt(abs(OmegaK))
  }
  CoDiff=rDist1*HubDist+CoSep12
  MaxCo=cosdistCoDist(1e5, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime)
  if(CoDiff>MaxCo){
    print('Objects not in common light cone!')
    out=NA
  }else{
    zem=cosmapval(CoDiff,'CoDist',H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0=w0, wprime=wprime, zrange = c(0,1e5), out='z')
    out=cosdistzeff(zref=z1,zem=zem)
  }
  return(out)
}

cosdistzeff=function(zref=1, zem=2){
   if (!all(zref > -1)) {
        stop("All zref must be > -1")
  }
  if (!all(zem > -1)) {
        stop("All zem must be > -1")
  }
  return((1+zem)/(1+zref)-1)
}
