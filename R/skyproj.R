radec2xy <-
function(RA,Dec,header,CRVAL1=0,CRVAL2=0,CRPIX1=0,CRPIX2=0,CD1_1=1,CD1_2=0,CD2_1=0,CD2_2=1,CTYPE1='RA--TAN',CTYPE2='DEC--TAN'){
# Converts RA/Dec (degrees) to x/y (pixels) position using the Tan Gnomonic projection or Sin Orthographic projection systems
# Translations adapted from: http://mathworld.wolfram.com/GnomonicProjection.html and http://mathworld.wolfram.com/OrthographicProjection.html
  if(length(dim(RA))==2){
    Dec=RA[,2]
    RA=RA[,1]
  }
  RA=as.numeric(RA)
  Dec=as.numeric(Dec)
  if(!missing(header)){
    if(is.data.frame(header) | is.matrix(header)){
    locs=match(c('CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','CDELT1','CDELT2'),header[,1])
    locs=locs[is.na(locs)==FALSE]
    headerWCS=data.frame(header[locs,1],as.character(header[locs,2]),stringsAsFactors=FALSE)
      if('CTYPE1' %in% headerWCS[,1]){CTYPE1=headerWCS[headerWCS[,1]=='CTYPE1',2]}else{message('Missing CTYPE1')}
      if('CTYPE2' %in% headerWCS[,1]){CTYPE2=headerWCS[headerWCS[,1]=='CTYPE2',2]}else{message('Missing CTYPE2')}
      if('CRVAL1' %in% headerWCS[,1]){CRVAL1=as.numeric(headerWCS[headerWCS[,1]=='CRVAL1',2])}else{message('Missing CRVAL1')}
      if('CRVAL2' %in% headerWCS[,1]){CRVAL2=as.numeric(headerWCS[headerWCS[,1]=='CRVAL2',2])}else{message('Missing CRVAL2')}
      if('CRPIX1' %in% headerWCS[,1]){CRPIX1=as.numeric(headerWCS[headerWCS[,1]=='CRPIX1',2])}else{message('Missing CRPIX1')}
      if('CRPIX2' %in% headerWCS[,1]){CRPIX2=as.numeric(headerWCS[headerWCS[,1]=='CRPIX2',2])}else{message('Missing CRPIX2')}
      if('CD1_1' %in% headerWCS[,1]){
        CD1_1=as.numeric(headerWCS[headerWCS[,1]=='CD1_1',2])
        if('CD1_2' %in% headerWCS[,1]){CD1_2=as.numeric(headerWCS[headerWCS[,1]=='CD1_2',2])}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% headerWCS[,1]){
          CD1_1=as.numeric(headerWCS[headerWCS[,1]=='CDELT1',2])
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% headerWCS[,1]){
        CD2_2=as.numeric(headerWCS[headerWCS[,1]=='CD2_2',2])
        if('CD2_1' %in% headerWCS[,1]){CD2_1=as.numeric(headerWCS[headerWCS[,1]=='CD2_1',2])}else{message('Missing CD2_1')}
      }else{
        if('CDELT2' %in% headerWCS[,1]){
          CD2_2=as.numeric(headerWCS[headerWCS[,1]=='CDELT2',2])
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }else{
      if('CTYPE1' %in% header){CTYPE1=as.character(header[which(header=='CTYPE1')+1])}else{message('Missing CTYPE1')}
      if('CTYPE2' %in% header){CTYPE2=as.character(header[which(header=='CTYPE2')+1])}else{message('Missing CTYPE2')}
      if('CRVAL1' %in% header){CRVAL1=as.numeric(header[which(header=='CRVAL1')+1])}else{message('Missing CRVAL1')}
      if('CRVAL2' %in% header){CRVAL2=as.numeric(header[which(header=='CRVAL2')+1])}else{message('Missing CRVAL2')}
      if('CRPIX1' %in% header){CRPIX1=as.numeric(header[which(header=='CRPIX1')+1])}else{message('Missing CRPIX1')}
      if('CRPIX2' %in% header){CRPIX2=as.numeric(header[which(header=='CRPIX2')+1])}else{message('Missing CRPIX2')}
      if('CD1_1' %in% header){
        CD1_1=as.numeric(header[which(header=='CD1_1')+1])
        if('CD1_2' %in% header){CD1_2=as.numeric(header[which(header=='CD1_2')+1])}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% header){
          CD1_1=as.numeric(header[which(header=='CDELT1')+1])
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% header){
        CD2_2=as.numeric(header[which(header=='CD2_2')+1])
        if('CD2_1' %in% header){CD2_1=as.numeric(header[which(header=='CD2_1')+1])}else{message('Missing CD2_1')}
      }else{
        if('CDELT1' %in% header){
          CD2_2=as.numeric(header[which(header=='CDELT2')+1])
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
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
  scalemat=matrix(c(x1,x2,y1,y2),2)*(pi/180)
  if(grepl('TAN', CTYPE1)){
    cosc1=sin(Dec0)*sin(Dec)+(cos(Dec0)*cos(Dec)*cos(RA-RA0))
  }else if(grepl('SIN', CTYPE1) | grepl('NCP', CTYPE1)){
    if(grepl('NCP', CTYPE1)){message('Approximating deprecated CTYPE1 NCP with SIN!')}
    cosc1=1
  }else{
    stop('Projection system is not recognised. Must be either TAN, SIN or NCP!')
  }
  if(grepl('TAN', CTYPE2)){
    cosc2=sin(Dec0)*sin(Dec)+(cos(Dec0)*cos(Dec)*cos(RA-RA0))
  }else if(grepl('SIN', CTYPE2) | grepl('NCP', CTYPE2)){
    if(grepl('NCP', CTYPE2)){message('Approximating deprecated CTYPE2 NCP with SIN!')}
    cosc2=1
  }else{
    stop('Projection system is not recognised. Must be either TAN, SIN or NCP!')
  }
  xxfunc = function(RA0,Dec0,RA,Dec){
          (cos(Dec)*sin(RA-RA0))/cosc1
  }
  yyfunc = function(RA0,Dec0,RA,Dec){
          ((cos(Dec0)*sin(Dec))-(sin(Dec0)*cos(Dec)*cos(RA-RA0)))/cosc2
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
function(x,y,header,CRVAL1=0,CRVAL2=0,CRPIX1=0,CRPIX2=0,CD1_1=1,CD1_2=0,CD2_1=0,CD2_2=1,CTYPE1='RA--TAN',CTYPE2='DEC--TAN') {
  # Converts x/y (pixels) to RA/DEC (degrees) position using the Tan Gnomonic or Sin Orthographic projection systems
  # Translations adapted from: http://mathworld.wolfram.com/GnomonicProjection.html and http://mathworld.wolfram.com/OrthographicProjection.html
  if(length(dim(x))==2){
    y=x[,2]
    x=x[,1]
  }
  x=as.numeric(x)
  y=as.numeric(y)
  if(!missing(header)){
    if(is.data.frame(header) | is.matrix(header)){
    locs=match(c('CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','CDELT1','CDELT2'),header[,1])
    locs=locs[is.na(locs)==FALSE]
    headerWCS=data.frame(header[locs,1],as.character(header[locs,2]),stringsAsFactors=FALSE)
      if('CTYPE1' %in% headerWCS[,1]){CTYPE1=headerWCS[headerWCS[,1]=='CTYPE1',2]}else{message('Missing CTYPE1')}
      if('CTYPE2' %in% headerWCS[,1]){CTYPE2=headerWCS[headerWCS[,1]=='CTYPE2',2]}else{message('Missing CTYPE2')}
      if('CRVAL1' %in% headerWCS[,1]){CRVAL1=as.numeric(headerWCS[headerWCS[,1]=='CRVAL1',2])}else{message('Missing CRVAL1')}
      if('CRVAL2' %in% headerWCS[,1]){CRVAL2=as.numeric(headerWCS[headerWCS[,1]=='CRVAL2',2])}else{message('Missing CRVAL2')}
      if('CRPIX1' %in% headerWCS[,1]){CRPIX1=as.numeric(headerWCS[headerWCS[,1]=='CRPIX1',2])}else{message('Missing CRPIX1')}
      if('CRPIX2' %in% headerWCS[,1]){CRPIX2=as.numeric(headerWCS[headerWCS[,1]=='CRPIX2',2])}else{message('Missing CRPIX2')}
      if('CD1_1' %in% headerWCS[,1]){
        CD1_1=as.numeric(headerWCS[headerWCS[,1]=='CD1_1',2])
        if('CD1_2' %in% headerWCS[,1]){CD1_2=as.numeric(headerWCS[headerWCS[,1]=='CD1_2',2])}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% headerWCS[,1]){
          CD1_1=as.numeric(headerWCS[headerWCS[,1]=='CDELT1',2])
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% headerWCS[,1]){
        CD2_2=as.numeric(headerWCS[headerWCS[,1]=='CD2_2',2])
        if('CD2_1' %in% headerWCS[,1]){CD2_1=as.numeric(headerWCS[headerWCS[,1]=='CD2_1',2])}else{message('Missing CD2_1')}
      }else{
        if('CDELT2' %in% headerWCS[,1]){
          CD2_2=as.numeric(headerWCS[headerWCS[,1]=='CDELT2',2])
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }else{
      if('CTYPE1' %in% header){CTYPE1=as.character(header[which(header=='CTYPE1')+1])}else{message('Missing CTYPE1')}
      if('CTYPE2' %in% header){CTYPE2=as.character(header[which(header=='CTYPE2')+1])}else{message('Missing CTYPE2')}
      if('CRVAL1' %in% header){CRVAL1=as.numeric(header[which(header=='CRVAL1')+1])}else{message('Missing CRVAL1')}
      if('CRVAL2' %in% header){CRVAL2=as.numeric(header[which(header=='CRVAL2')+1])}else{message('Missing CRVAL2')}
      if('CRPIX1' %in% header){CRPIX1=as.numeric(header[which(header=='CRPIX1')+1])}else{message('Missing CRPIX1')}
      if('CRPIX2' %in% header){CRPIX2=as.numeric(header[which(header=='CRPIX2')+1])}else{message('Missing CRPIX2')}
      if('CD1_1' %in% header){
        CD1_1=as.numeric(header[which(header=='CD1_1')+1])
        if('CD1_2' %in% header){CD1_2=as.numeric(header[which(header=='CD1_2')+1])}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% header){
          CD1_1=as.numeric(header[which(header=='CDELT1')+1])
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% header){
        CD2_2=as.numeric(header[which(header=='CD2_2')+1])
        if('CD2_1' %in% header){CD2_1=as.numeric(header[which(header=='CD2_1')+1])}else{message('Missing CD2_1')}
      }else{
        if('CDELT1' %in% header){
          CD2_2=as.numeric(header[which(header=='CDELT2')+1])
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
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
  scalemat=matrix(c(x1,x2,y1,y2),2)*(pi/180)
  xytran=cbind(x-x0,y-y0) %*% scalemat
  x = xytran[,1]
  y = xytran[,2]
  rad = sqrt(x^2+y^2)
  if(grepl('TAN', CTYPE1)){
    radproj1=atan(rad)
  }else if(grepl('SIN', CTYPE1) | grepl('NCP', CTYPE1)){
    if(grepl('NCP', CTYPE1)){message('Approximating deprecated CTYPE1 NCP with SIN!')}
    radproj1=asin(rad)
  }else{
    stop('Projection system is not recognised. Must be either TAN, SIN or NCP!')
  }
  if(grepl('TAN', CTYPE2)){
    radproj2=atan(rad)
  }else if(grepl('SIN', CTYPE2) | grepl('NCP', CTYPE2)){
    if(grepl('NCP', CTYPE2)){message('Approximating deprecated CTYPE2 NCP with SIN!')}
    radproj2=asin(rad)
  }else{
    stop('Projection system is not recognised. Must be either TAN, SIN or NCP!')
  }
  rafunc = function(RA0,Dec0,x,y){
      RA0 + atan2(x*sin(radproj1),rad*cos(Dec0)*cos(radproj1) - y*sin(Dec0)*sin(radproj1))
  }
  decfunc = function(Dec0,x,y){
      asin(cos(radproj2)*sin(Dec0) + (y*sin(radproj2)*cos(Dec0) / rad)) 
  }
  RA = rafunc(RA0,Dec0,x,y)*180/pi %% 360
  Dec = decfunc(Dec0,x,y)*180/pi %% 90
  Dec[which(is.nan(Dec))] = Dec0*180/pi
  output=cbind(as.numeric(RA),as.numeric(Dec))
  colnames(output)=c('RA','Dec')
  return(output)
}

getpixscale=function(header, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1){
  if(!missing(header)){
    if(is.data.frame(header) | is.matrix(header)){
    if(length(header)>1e4){
      stop('Header does not seem to be legal- far too long!')
    }
    locs=match(c('CD1_1','CD1_2','CD2_1','CD2_2','CDELT1','CDELT2'),header[,1])
    locs=locs[is.na(locs)==FALSE]
    headerWCS=data.frame(header[locs,1],as.numeric(header[locs,2]),stringsAsFactors=FALSE)
      if('CD1_1' %in% headerWCS[,1]){
        CD1_1=headerWCS[headerWCS[,1]=='CD1_1',2]
        if('CD1_2' %in% headerWCS[,1]){CD1_2=headerWCS[headerWCS[,1]=='CD1_2',2]}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% headerWCS[,1]){
          CD1_1=headerWCS[headerWCS[,1]=='CDELT1',2]
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% headerWCS[,1]){
        CD2_2=headerWCS[headerWCS[,1]=='CD2_2',2]
        if('CD2_1' %in% headerWCS[,1]){CD2_1=headerWCS[headerWCS[,1]=='CD2_1',2]}else{message('Missing CD2_1')}
      }else{
        if('CDELT2' %in% headerWCS[,1]){
          CD2_2=headerWCS[headerWCS[,1]=='CDELT2',2]
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }else{
      if('CD1_1' %in% header){
        CD1_1=as.numeric(header[which(header=='CD1_1')+1])
        if('CD1_2' %in% header){CD1_2=as.numeric(header[which(header=='CD1_2')+1])}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% header){
          CD1_1=as.numeric(header[which(header=='CDELT1')+1])
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% header){
        CD2_2=as.numeric(header[which(header=='CD2_2')+1])
        if('CD2_1' %in% header){CD2_1=as.numeric(header[which(header=='CD2_1')+1])}else{message('Missing CD2_1')}
      }else{
        if('CDELT1' %in% header){
          CD2_2=as.numeric(header[which(header=='CDELT2')+1])
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
  }
  return(3600*(sqrt(CD1_1^2+CD1_2^2)+sqrt(CD2_1^2+CD2_2^2))/2)
}
