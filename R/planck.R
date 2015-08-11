cosplanckCMBTemp=function(z){
  return(2.725*(1+z))
}

cosplanckLawRadFreq=function(nu,Temp=2.725){
  A=2*kNIST2010PlanckConstant*nu^3/kNIST2010speedOfLightInVacuum^2
  return(A/(exp((kNIST2010PlanckConstant*nu)/(kNIST2010BoltzmannConstant*Temp))-1))
}

cosplanckLawRadWave=function(lambda,Temp=2.725){
  A=2*kNIST2010PlanckConstant*kNIST2010speedOfLightInVacuum^2/lambda^5
  return(A/(exp((kNIST2010PlanckConstant*kNIST2010speedOfLightInVacuum)/(lambda*kNIST2010BoltzmannConstant*Temp))-1))
}

cosplanckLawEnFreq=function(nu,Temp=2.725){
  A=2*kNIST2010PlanckConstant*nu^3/kNIST2010speedOfLightInVacuum^2
  return((4*pi/kNIST2010speedOfLightInVacuum)*A/(exp((kNIST2010PlanckConstant*nu)/(kNIST2010BoltzmannConstant*Temp))-1))
}

cosplanckLawEnWave=function(lambda,Temp=2.725){
  A=2*kNIST2010PlanckConstant*kNIST2010speedOfLightInVacuum^2/lambda^5
  return((4*pi/kNIST2010speedOfLightInVacuum)*A/(exp((kNIST2010PlanckConstant*kNIST2010speedOfLightInVacuum)/(lambda*kNIST2010BoltzmannConstant*Temp))-1))
}

cosplanckLawRadFreqN=function(nu,Temp=2.725){
  A=2*kNIST2010PlanckConstant*nu^3/kNIST2010speedOfLightInVacuum^2
  return((1/(kNIST2010PlanckConstant*nu))*A/(exp((kNIST2010PlanckConstant*nu)/(kNIST2010BoltzmannConstant*Temp))-1))
}

cosplanckLawRadWaveN=function(lambda,Temp=2.725){
  A=2*kNIST2010PlanckConstant*kNIST2010speedOfLightInVacuum^2/lambda^5
  nu=kNIST2010speedOfLightInVacuum/lambda
  return((1/(kNIST2010PlanckConstant*nu))*A/(exp((kNIST2010PlanckConstant*kNIST2010speedOfLightInVacuum)/(lambda*kNIST2010BoltzmannConstant*Temp))-1))
}

cosplanckPeakFreq=function(Temp=2.725){
  EnPeakNu=2.821*kNIST2010BoltzmannConstant*Temp
  return(EnPeakNu/kNIST2010PlanckConstant)
}

cosplanckPeakWave=function(Temp=2.725){
  EnPeakLambda=4.965*kNIST2010BoltzmannConstant*Temp
  return(kNIST2010speedOfLightInVacuum*kNIST2010PlanckConstant/EnPeakLambda)
}

cosplanckSBLawRad_sr=function(Temp=2.725){
  return(kNIST2010StefanBoltzmannConstant*Temp^4/pi)
}

cosplanckSBLawRad=function(Temp=2.725){
  return(kNIST2010StefanBoltzmannConstant*Temp^4)
}

cosplanckSBLawEn=function(Temp=2.725){
  return((4/kNIST2010speedOfLightInVacuum)*kNIST2010StefanBoltzmannConstant*Temp^4)
}

cosplanckLawRadPhotEnAv=function(Temp=2.725){
  return(3.729282e-23*Temp)
}

cosplanckLawRadPhotN=function(Temp=2.725){
  return(1.5205e15*Temp^3/pi)
}