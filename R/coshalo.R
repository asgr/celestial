coshaloMvirToSigma=function(Mvir=1e12, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g*Munit/(Lunit*Vunit^2)
  if(Rho=='crit'){RhoVal=(cosgrowRhoCrit(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Rho=='mean'){RhoVal=(cosgrowRhoMean(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  return((Mvir*sqrt(32*pi*g^3*DeltaVir*RhoVal/3))^(1/3))
}

coshaloSigmaToMvir=function(Sigma=230, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g*Munit/(Lunit*Vunit^2)
  if(Rho=='crit'){RhoVal=(cosgrowRhoCrit(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Rho=='mean'){RhoVal=(cosgrowRhoMean(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  return(Sigma^3*sqrt(3/(32*pi*g^3*DeltaVir*RhoVal)))
}

coshaloMvirToRvir=function(Mvir=1e12, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g*Munit/(Lunit*Vunit^2)
  if(Rho=='crit'){RhoVal=(cosgrowRhoCrit(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist=Dist, ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Rho=='mean'){RhoVal=(cosgrowRhoMean(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist=Dist, ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  return(((3*Mvir)/(4*pi*DeltaVir*RhoVal))^(1/3))
}

coshaloRvirToMvir=function(Rvir=162.635, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g*Munit/(Lunit*Vunit^2)
  if(Rho=='crit'){RhoVal=(cosgrowRhoCrit(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist=Dist, ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Rho=='mean'){RhoVal=(cosgrowRhoMean(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist=Dist, ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  return((4*pi/3)*DeltaVir*RhoVal*Rvir^3)
}

coshaloSigmaToRvir=function(Sigma=230, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g*Munit/(Lunit*Vunit^2)
  if(Rho=='crit'){RhoVal=(cosgrowRhoCrit(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Rho=='mean'){RhoVal=(cosgrowRhoMean(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Dist=='Ang'){scale=1}
  if(Dist=='Co'){scale=1+z}
  return(scale*Sigma*(27/(512*pi^3*g^3*DeltaVir^3*RhoVal^3))^(1/6))
}

coshaloRvirToSigma=function(Rvir=162.635, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  g = G*msol_to_kg/(pc_to_m)
  g = g*Munit/(Lunit*Vunit^2)
  if(Rho=='crit'){RhoVal=(cosgrowRhoCrit(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Rho=='mean'){RhoVal=(cosgrowRhoMean(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, Dist='Ang', ref=ref)/1e9)*((Lunit/1e3)^3)/Munit}
  if(Dist=='Ang'){scale=1}
  if(Dist=='Co'){scale=1+z}
  return((Rvir/scale)*((512*pi^3*g^3*DeltaVir^3*RhoVal^3)/27)^(1/6))
}

