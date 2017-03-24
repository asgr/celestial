radec2xy <-
function(RA,Dec,CRVAL1=0,CRVAL2=0,CRPIX1=0,CRPIX2=0,CD1_1=1,CD1_2=0,CD2_1=0,CD2_2=1,header){
# Converts RA/Dec (degrees) to x/y (pixels) position using the Tan Gnomonic projection system
# Translations adapted from: http://mathworld.wolfram.com/GnomonicProjection.html
  if(length(dim(RA))==2){
    Dec=RA[,2]
    RA=RA[,1]
  }
  RA=as.numeric(RA)
  Dec=as.numeric(Dec)
  if(!missing(header)){
    locs=match(c('CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2'),header[,1])
    headerWCS=data.frame(header[locs,1],as.numeric(header[locs,2]))
    if('CRVAL1' %in% headerWCS[,1]){CRVAL1=headerWCS[headerWCS[,1]=='CRVAL1',2]}else{print('Missing CRVAL1')}
    if('CRVAL2' %in% headerWCS[,1]){CRVAL2=headerWCS[headerWCS[,1]=='CRVAL2',2]}else{print('Missing CRVAL1')}
    if('CRPIX1' %in% headerWCS[,1]){CRPIX1=headerWCS[headerWCS[,1]=='CRPIX1',2]}else{print('Missing CRPIX1')}
    if('CRPIX2' %in% headerWCS[,1]){CRPIX2=headerWCS[headerWCS[,1]=='CRPIX2',2]}else{print('Missing CRPIX2')}
    if('CD1_1' %in% headerWCS[,1]){CD1_1=headerWCS[headerWCS[,1]=='CD1_1',2]}else{print('Missing CD1_1')}
    if('CD1_2' %in% headerWCS[,1]){CD1_2=headerWCS[headerWCS[,1]=='CD1_2',2]}else{print('Missing CD1_2')}
    if('CD2_1' %in% headerWCS[,1]){CD2_1=headerWCS[headerWCS[,1]=='CD2_1',2]}else{print('Missing CD2_1')}
    if('CD2_2' %in% headerWCS[,1]){CD2_2=headerWCS[headerWCS[,1]=='CD2_2',2]}else{print('Missing CD2_2')}
  }
  RA0=CRVAL1
  Dec0=CRVAL2
  x0=CRPIX1
  y0=CRPIX2
  x1=CD1_1
  x2=CD1_2
  y1=CD2_1
  y2=CD2_2
  RA0=RA0*(pi/180)
  Dec0=Dec0*(pi/180)
  RA=RA*(pi/180)
  Dec=Dec*(pi/180)
  scalemat=tan(matrix(c(x1,x2,y1,y2),2)*(pi/180))
  cosc=sin(Dec0)*sin(Dec)+(cos(Dec0)*cos(Dec)*cos(RA-RA0))
  xxfunc = function(RA0,Dec0,RA,Dec){
          (cos(Dec)*sin(RA-RA0))/cosc
  }
  yyfunc = function(RA0,Dec0,RA,Dec){
          ((cos(Dec0)*sin(Dec))-(sin(Dec0)*cos(Dec)*cos(RA-RA0)))/cosc
  }
  XX=xxfunc(RA0,Dec0,RA,Dec)
  YY=yyfunc(RA0,Dec0,RA,Dec)
  raw=cbind(XX,YY)
  output=raw %*% solve(scalemat)
  output[,1]=output[,1]+x0
  output[,2]=output[,2]+y0
  colnames(output)=c('x','y')
  return(output)
}

xy2radec <-
function(x,y,CRVAL1=0,CRVAL2=0,CRPIX1=0,CRPIX2=0,CD1_1=1,CD1_2=0,CD2_1=0,CD2_2=1,header) {
  # Converts x/y (pixels) to RA/DEC (degrees) position using the Tan Gnomonic projection system
  # Translations adapted from: http://mathworld.wolfram.com/GnomonicProjection.html
  if(length(dim(x))==2){
    y=x[,2]
    x=x[,1]
  }
  x=as.numeric(x)
  y=as.numeric(y)
  if(!missing(header)){
    locs=match(c('CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2'),header[,1])
    headerWCS=data.frame(header[locs,1],as.numeric(header[locs,2]))
    if('CRVAL1' %in% headerWCS[,1]){CRVAL1=headerWCS[headerWCS[,1]=='CRVAL1',2]}else{print('Missing CRVAL1')}
    if('CRVAL2' %in% headerWCS[,1]){CRVAL2=headerWCS[headerWCS[,1]=='CRVAL2',2]}else{print('Missing CRVAL1')}
    if('CRPIX1' %in% headerWCS[,1]){CRPIX1=headerWCS[headerWCS[,1]=='CRPIX1',2]}else{print('Missing CRPIX1')}
    if('CRPIX2' %in% headerWCS[,1]){CRPIX2=headerWCS[headerWCS[,1]=='CRPIX2',2]}else{print('Missing CRPIX2')}
    if('CD1_1' %in% headerWCS[,1]){CD1_1=headerWCS[headerWCS[,1]=='CD1_1',2]}else{print('Missing CD1_1')}
    if('CD1_2' %in% headerWCS[,1]){CD1_2=headerWCS[headerWCS[,1]=='CD1_2',2]}else{print('Missing CD1_2')}
    if('CD2_1' %in% headerWCS[,1]){CD2_1=headerWCS[headerWCS[,1]=='CD2_1',2]}else{print('Missing CD2_1')}
    if('CD2_2' %in% headerWCS[,1]){CD2_2=headerWCS[headerWCS[,1]=='CD2_2',2]}else{print('Missing CD2_2')}
  }
  RA0=CRVAL1
  Dec0=CRVAL2
  x0=CRPIX1
  y0=CRPIX2
  x1=CD1_1
  x2=CD1_2
  y1=CD2_1
  y2=CD2_2
  RA0=RA0*(pi/180)
  Dec0=Dec0*(pi/180)
  scalemat=tan(matrix(c(x1,x2,y1,y2),2)*(pi/180))
  xytran=cbind(x-x0,y-y0) %*% scalemat
  x = xytran[,1]
  y = xytran[,2]
  rafunc = function(RA0,Dec0,x,y){
      RA0 + atan2(x*sin(atan(sqrt(x^2+y^2))),sqrt(x^2+y^2)*cos(Dec0)*cos(atan(sqrt(x^2+y^2))) - y*sin(Dec0)*sin(atan(sqrt(x^2+y^2))))
  }
  decfunc = function(Dec0,x,y){
      asin(cos(atan(sqrt(x^2+y^2)))*sin(Dec0) + (y*sin(atan(sqrt(x^2+y^2)))*cos(Dec0) / sqrt(x^2+y^2))) 
  }
  RA = rafunc(RA0,Dec0,x,y)*180/pi
  Dec = decfunc(Dec0,x,y)*180/pi
  Dec[which(is.nan(Dec))] = Dec0*180/pi
  output=cbind(as.numeric(RA),as.numeric(Dec))
  colnames(output)=c('RA','Dec')
  return(output)
}