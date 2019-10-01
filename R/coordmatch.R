coordmatch=function(coordref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg", radunit='asec', sep = ":", kstart=10, ignoreexact=FALSE, ignoreinternal=FALSE, matchextra=FALSE, smallapprox=FALSE){
  if (inunitref %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitref must be one of deg, rad or sex")
  }
  if (inunitcompare %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitcompare must be one of deg, rad or sex")
  }
  if (radunit %in% c("deg", "amin", "asec", "rad") == FALSE) {
    stop("radunit must be one of deg, amin, asec or rad")
  }
  origrad=rad
  coordref=rbind(coordref)
  if(missing(coordcompare)){
    coordcompare=coordref
    if(missing(ignoreinternal)){
      ignoreinternal=TRUE
    }
  }
  coordcompare=rbind(coordcompare)
  N=length(coordref[,1])
  kmax=length(coordcompare[,1])
  if(length(rad)==1){rad=rep(rad,N)}
  if(length(rad) != N){stop("Length of rad must either be 1 or the same length as rows in coordref")}
  if (inunitref == "sex") {
    coordref[,1] = hms2deg(coordref[,1], sep = sep)
    coordref[,2] = dms2deg(coordref[,2], sep = sep)
    coordref=cbind(as.numeric(coordref[,1]),as.numeric(coordref[,2]))
  }
  if (inunitref == "rad") {
    coordref[,1] = coordref[,1] * 180/pi
    coordref[,2] = coordref[,2] * 180/pi
  }
  if (inunitcompare == "sex") {
    coordcompare[,1] = hms2deg(coordcompare[,1], sep = sep)
    coordcompare[,2] = dms2deg(coordcompare[,2], sep = sep)
    coordcompare=cbind(as.numeric(coordcompare[,1]),as.numeric(coordcompare[,2]))
  }
  if (inunitcompare == "rad") {
    coordcompare[,1] = coordcompare[,1] * 180/pi
    coordcompare[,2] = coordcompare[,2] * 180/pi
  }
  if (radunit == "asec"){
    radmult=(pi/180)/3600
  }
  if (radunit == "amin"){
    radmult=(pi/180)/60
  }
  if (radunit == "deg"){
    radmult=(pi/180)
  }
  rad=rad*radmult
  userad=max(rad,na.rm = TRUE)
  
  coordrefxyz=sph2car(coordref[,1:2,drop=FALSE],deg=TRUE)
  coordcomparexyz=sph2car(coordcompare[,1:2,drop=FALSE],deg=TRUE)
  
  if(matchextra & dim(coordref)[2]>2 & dim(coordcompare)[2]>2){
    if(dim(coordref)[2] == dim(coordcompare)[2]){
      coordrefxyz=cbind(coordrefxyz,coordref[,3:dim(coordref)[2],drop=FALSE]*radmult)
      coordcomparexyz=cbind(coordcomparexyz,coordcompare[,3:dim(coordcompare)[2],drop=FALSE]*radmult)
    }
  }
  
  ksuggest=min(kstart, dim(coordcomparexyz)[1])
  while(is.na(ksuggest)==FALSE){
    tempmatch=nn2(coordcomparexyz,coordrefxyz,searchtype='radius',radius=userad,k=ksuggest)
    ignore=tempmatch[[1]]==0
    tempmatch[[2]][ignore]=NA
    if(smallapprox==FALSE){
      tempmatch[[2]]=2*asin(tempmatch[[2]]/2)
    }
    if(ignoreinternal){
      remove=which(tempmatch[[1]]-1:length(coordcomparexyz[,1])==0)
      tempmatch[[1]][remove]=0
      tempmatch[[2]][remove]=NA
    }
    if(ignoreexact){
      remove=which(!(tempmatch[[2]]<=rad & tempmatch[[2]]>0))
      tempmatch[[1]][remove]=0
      tempmatch[[2]][remove]=NA
    }else{
      remove=which(!tempmatch[[2]]<=rad)
      tempmatch[[1]][remove]=0
      tempmatch[[2]][remove]=NA
    }
    if (all(is.na(tempmatch[[2]][,ksuggest]))){
      kendmin=NA
    }else{
      kendmin=min(tempmatch[[2]][,ksuggest],na.rm = TRUE)
    }
    if(is.na(kendmin)==FALSE & ksuggest<kmax){
      comp=tempmatch[[2]][,ksuggest]
      ksuggest=ceiling(ksuggest*max(rad[comp>0]/comp[comp>0],na.rm=TRUE))
      ksuggest=min(ksuggest,kmax)
    }else{
      ksuggest=NA
    }
  }

  keepcols=which(colSums(is.na(tempmatch[[2]]))<N)
  tempmatch[[1]]=matrix(tempmatch[[1]][,keepcols],nrow=N,byrow=FALSE)
  tempmatch[[2]]=matrix(tempmatch[[2]][,keepcols],nrow=N,byrow=FALSE)
  
  if (radunit == "asec"){
    tempmatch[[2]]=tempmatch[[2]]/((pi/180)/3600)
  }
  if (radunit == "amin"){
    tempmatch[[2]]=tempmatch[[2]]/((pi/180)/60)
  }
  if (radunit == "deg"){
    tempmatch[[2]]=tempmatch[[2]]/((pi/180))
  }
  
  if(length(tempmatch[[1]])>0){
    Nmatch=rowSums(tempmatch[[1]]!=0)
  }else{
    Nmatch=NA
  }
  
  if(is.na(Nmatch[1])==FALSE){
    bestID=tempmatch[[1]][,1]
    bestsep=tempmatch[[2]][,1]
    select=which(bestID>0)
    bestsep=bestsep[select]
    bestID=bestID[select]
    orderID=order(bestsep)
    keep=!duplicated(bestID[orderID])
    bestrefID=select[orderID][keep]
    bestcompareID=bestID[orderID][keep]
    bestsep=bestsep[orderID][keep]
    reorderID=order(bestrefID)
    bestrefID=bestrefID[reorderID]
    bestcompareID=bestcompareID[reorderID]
    bestsep=bestsep[reorderID]
    bestmatch=data.frame(refID=bestrefID, compareID=bestcompareID, sep=bestsep)
    output=list(ID=tempmatch[[1]], sep=tempmatch[[2]], Nmatch=Nmatch, bestmatch=bestmatch)
  }else{
    output=list(ID=tempmatch[[1]], sep=tempmatch[[2]], Nmatch=Nmatch, bestmatch=NA)
  }
  return(invisible(output))
}

coordmatchsing=function(RAref,Decref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg", radunit='asec', sep = ":", ignoreexact=FALSE, smallapprox=FALSE){
  if (inunitref %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitref must be one of deg, rad or sex")
  }
  if (inunitcompare %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitcompare must be one of deg, rad or sex")
  }
  coordref=rbind(c(RAref,Decref))
  coordcompare=rbind(coordcompare)
  N=length(coordref[,1])
  if(N>1){stop("Only 1 RAref and Decref value allowed!")}
  if (inunitref == "sex") {
    coordref[,1] = hms2deg(coordref[,1], sep = sep)
    coordref[,2] = dms2deg(coordref[,2], sep = sep)
    coordref=cbind(as.numeric(coordref[,1]),as.numeric(coordref[,2]))
  }
  if (inunitref == "rad") {
    coordref[,1] = coordref[,1] * 180/pi
    coordref[,2] = coordref[,2] * 180/pi
  }
  if (inunitcompare == "sex") {
    coordcompare[,1] = hms2deg(coordcompare[,1], sep = sep)
    coordcompare[,2] = dms2deg(coordcompare[,2], sep = sep)
    coordcompare=cbind(as.numeric(coordcompare[,1]),as.numeric(coordcompare[,2]))
  }
  if (inunitcompare == "rad") {
    coordcompare[,1] = coordcompare[,1] * 180/pi
    coordcompare[,2] = coordcompare[,2] * 180/pi
  }
  
  coordrefxyz=sph2car(coordref,deg=TRUE)
  coordcomparexyz=sph2car(coordcompare,deg=TRUE)
  
  dotprod=coordcomparexyz[,1]*coordrefxyz[1]+coordcomparexyz[,2]*coordrefxyz[2]+coordcomparexyz[,3]*coordrefxyz[3]
  dotprod[dotprod< -1]=-1
  dotprod[dotprod>1]=1
  if(smallapprox==FALSE){
    ang=acos(dotprod)
  }else{
    ang=pi/2-dotprod
  }
  if (radunit == "asec"){
    ang=ang/((pi/180)/3600)
  }
  if (radunit == "amin"){
    ang=ang/((pi/180)/60)
  }
  if (radunit == "deg"){
    ang=ang/(pi/180)
  }
  if(ignoreexact){select=which(ang<=rad & ang>0)}else{select=which(ang<=rad)}
  if(length(select)==0){
    output=list(ID=NA, sep=NA, Nmatch=NA, bestmatch=NA)
  }else{
    ID=select
    sep=ang[select]
    orderID=order(sep)
    ID=ID[orderID]
    sep=sep[orderID]
    bestmatch=c(compareID=select[which.min(ang[select])],sep=min(ang[select]))
    output=list(ID=ID, sep=sep, Nmatch=length(select), bestmatch=bestmatch)
  }
  return(invisible(output))
}

internalclean=function(RA, Dec, rad=2, tiebreak, decreasing = FALSE, inunit="deg", radunit='asec', sep = ":", Nmatch='all'){
  
  if (length(dim(RA)) == 2){
    RA=as.matrix(RA)
    if(dim(RA)[2]>=2){
      Dec = RA[, 2]
    }
    if(dim(RA)[2]>=3){
      tiebreak = RA[, 3]
    }
    RA = RA[, 1]
  }
  
  RA = as.numeric(RA)
  Dec = as.numeric(Dec)
  
  if(missing(tiebreak)){
    tiebreak=1:length(RA)
  }
  
  bestfunc=function(x){
    if(all(x==0)){
      return(0)
    }else{
      return(min(x[x>0],na.rm = TRUE))
    }
  }
  
  matchorder=order(tiebreak, decreasing=decreasing)
  match=coordmatch(cbind(RA[matchorder],Dec[matchorder]), rad=rad, inunitref = inunit, inunitcompare = inunit, radunit = radunit, sep = sep)
  nearcen=apply(cbind(match$bestmatch[,1],match$ID[match$bestmatch[,1],]),1,bestfunc)
  if(Nmatch[1]=='all'){
    output=matchorder[c(which(match$Nmatch==0),nearcen)]
  }else{
    output=matchorder[nearcen[match$Nmatch(match$bestmatch$ref) %in% Nmatch]]
  }
  return(invisible(sort(unique(output))))
}
