cosdist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, age=FALSE){
  OmegaK=1-OmegaM-OmegaL
  Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
  if(age){Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    
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
      DistMod = 5*log10(LumDist)+25
      AngSize = AngDist*(pi/(180*60*60))*1000
      
      if (age) {
        HT = (3.08568025e+19/(H0*31556926))/1e9
        UniAge = HT*integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
        zAge = HT*integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
      }
      if (age) {
        return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVol = CoVol, HubTime = HT, UniAgeNow = UniAge, UniAgeAtz = UniAge - zAge, TravelTime = zAge)
      }
      else {
        return = c(z = z, a = a, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVol = CoVol)
      }
    }
    return(as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))))
}

cosdistz=function(z = 1){
  return(z)
}

cosdista=function(z = 1){
  return(1/(1 + z))
}

cosdistCoDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    return=CoDist
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistCoDistTran=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
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
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistLumDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
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
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistAngDist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
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
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistDistMod=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
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
    DistMod = 5 * log10(CoDistTran * (1 + z)) + 25
    return=DistMod
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistAngSize=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
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
    AngSize = (CoDistTran / (1 + z)) * (pi/(180 * 60 * 60)) * 1000
    return=AngSize
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistCoVol=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
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
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistUniAgeNow=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    UniAge = HT * integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    return=UniAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistUniAgeAtz=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    UniAge = HT * integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    zAge = HT * integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    return=UniAge-zAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistTravelTime=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    zAge = HT * integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    return=zAge
  }
  return(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK))
}

cosdistHubTime=function(H0 = 100){
 return((3.08568025e+19/(H0 * 31556926))/1e9)
}


