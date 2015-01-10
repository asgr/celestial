cosdist=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, age=FALSE){
  OmegaK=1-OmegaM-OmegaL
  temp = function(z, H0, OmegaM, OmegaL, OmegaK) {
    Einv = function(z, OmegaM, OmegaL, OmegaK) {1/sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL)}
    HubDist = (299792.458/H0)
    CoDist = HubDist * integrate(Einv, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
    
      if(OmegaK==0){
        CoDistTran = CoDist
        CoVol = (4/3) * pi * CoDist^3
      }else{
        if(OmegaK>0){
          CoDistTran = HubDist*(1/sqrt(OmegaK))*sinh(sqrt(OmegaK)*CoDist/HubDist)
          CoVol = (4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asinh(sqrt(abs(OmegaK))*(CoDistTran/HubDist)))
        }
        if(OmegaK<0){
          CoDistTran = HubDist*(1/sqrt(abs(OmegaK)))*sin(sqrt(abs(OmegaK))*CoDist/HubDist)
          CoVol = (4*pi*HubDist^3/(2*OmegaK))*((CoDistTran/HubDist)*sqrt(1+OmegaK*(CoDistTran/HubDist)^2)-(1/sqrt(abs(OmegaK)))*asin(sqrt(abs(OmegaK))*(CoDistTran/HubDist)))
        }
      }
      
      AngDist = CoDistTran/(1 + z)
      LumDist = (1+z) * CoDistTran
      DistMod = 5 * log10(LumDist) + 25
      AngSize = AngDist * (pi/(180 * 60 * 60)) * 1000
      
      if (age) {
        Einvz = function(z, OmegaM, OmegaL, OmegaK){1/(sqrt(OmegaM * (1 + z)^3 + OmegaK * (1 + z)^2 + OmegaL) * (1 + z))}
        HT = 3.08568025e+19/(H0 * 31556926)
        UniAge = HT * integrate(Einvz, 0, Inf, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
        zAge = HT * integrate(Einvz, 0, z, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK, subdivisions = 1000)$value
      }
      if (age) {
        return = c(z = z, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVolGpc3 = CoVol/1e+09, HubTime = HT, UniAgeNow = UniAge, UniAgeAtz = UniAge - zAge, TravelTime = zAge)
      }
      else {
        return = c(z = z, CoDist = CoDist, LumDist = LumDist, AngDist = AngDist, CoDistTran=CoDistTran, DistMod = DistMod, AngSize = AngSize, CoVolGpc3 = CoVol/1e+09)
      }
    }
    return = as.data.frame(t(Vectorize(temp)(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaK = OmegaK)))
}

Hz=function(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  return=H0*sqrt(OmegaM*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL)
}

rhocrit=function(z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM){
  OmegaK=1-OmegaM-OmegaL
  G=6.67384e-11 # m^3 kg^-1 s^-2
  Hub2=Hz(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)^2 # (km/s / Mpc)^2
  km2m=1000
  Mpc2m=3.08567758e22
  Msol2kg=1.9891e30 # kg
  #rhocrit_kgperm3=(3*Hub2)/(8*pi*G) * ((km2m^2)/(Mpc2m^2)) #this is correct, should be ~9.2e-27 for H0=70
  #rhocrit_MsolperMpc3=rhocrit_kgperm3 * (Mpc2m^3)/Msol2kg #this is correct, should be ~2.8e11 for H0=100
  rhocrit_MsolperMpc3=(3*Hub2)/(8*pi*G) * (km2m^2)*Mpc2m/Msol2kg #compact form of the above
  return=rhocrit_MsolperMpc3
}
