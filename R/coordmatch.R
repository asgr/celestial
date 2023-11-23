coordmatch = function(coordref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg",
                    radunit='asec', sep = ":", kstart=10, ignoreexact=FALSE, ignoreinternal=FALSE,
                    matchextra=FALSE, smallapprox=FALSE, jitter=FALSE, jamount=1e-12, jseed=666){
  if (inunitref %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitref must be one of deg, rad or sex")
  }
  if (inunitcompare %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitcompare must be one of deg, rad or sex")
  }
  if (radunit %in% c("deg", "amin", "asec", "rad") == FALSE) {
    stop("radunit must be one of deg, amin, asec or rad")
  }
  origrad = rad
  coordref = rbind(coordref)
  if(missing(coordcompare)){
    coordcompare = coordref
    if(missing(ignoreinternal)){
      ignoreinternal = TRUE
    }
  }
  coordcompare = rbind(coordcompare)
  N = length(coordref[, 1])
  kmax = length(coordcompare[, 1])
  if(length(rad)==1){rad=rep(rad,N)}
  if(length(rad) != N){stop("Length of rad must either be 1 or the same length as rows in coordref")}
  if (inunitref == "sex") {
    coordref[,1] = hms2deg(coordref[,1], sep = sep)
    coordref[,2] = dms2deg(coordref[,2], sep = sep)
    coordref = cbind(as.numeric(coordref[, 1]), as.numeric(coordref[, 2]))
  }
  if (inunitref == "rad") {
    coordref[,1] = coordref[,1] * 180/pi
    coordref[,2] = coordref[,2] * 180/pi
  }
  if (inunitcompare == "sex") {
    coordcompare[,1] = hms2deg(coordcompare[,1], sep = sep)
    coordcompare[,2] = dms2deg(coordcompare[,2], sep = sep)
    coordcompare = cbind(as.numeric(coordcompare[, 1]), as.numeric(coordcompare[, 2]))
  }
  if (inunitcompare == "rad") {
    coordcompare[,1] = coordcompare[,1] * 180/pi
    coordcompare[,2] = coordcompare[,2] * 180/pi
  }
  if (radunit == "asec"){
    radmult = (pi/180)/3600
  }
  if (radunit == "amin"){
    radmult = (pi/180)/60
  }
  if (radunit == "deg"){
    radmult = (pi/180)
  }
  rad = rad * radmult
  userad = max(rad, na.rm = TRUE)
  
  coordrefxyz = sph2car(coordref[, 1:2, drop = FALSE], deg = TRUE)
  coordcomparexyz = sph2car(coordcompare[, 1:2, drop = FALSE], deg = TRUE)
  
  if(missing(jitter) & ignoreinternal & ignoreexact==FALSE){
    jitter = TRUE
  }
  
  if(jitter){
    set.seed(jseed)
    jit_temp = runif(dim(coordcomparexyz)[1]*3, min=-jamount, max=jamount) #3 because we have 3 columns
    
    if(ignoreinternal){
      coordrefxyz = coordrefxyz + jit_temp
      coordcomparexyz = coordcomparexyz + jit_temp
    }else{
      coordcomparexyz = coordcomparexyz + jit_temp
    }
  }
  
  if(matchextra & dim(coordref)[2]>2 & dim(coordcompare)[2]>2){
    if(dim(coordref)[2] == dim(coordcompare)[2]){
      coordrefxyz = cbind(coordrefxyz, coordref[, 3:dim(coordref)[2], drop = FALSE] * radmult)
      coordcomparexyz = cbind(coordcomparexyz, coordcompare[, 3:dim(coordcompare)[2],
                                                            drop = FALSE] * radmult)
    }
  }
  
  ksuggest = min(kstart, dim(coordcomparexyz)[1])
  while(is.na(ksuggest)==FALSE){
    tempmatch = nn2(
      coordcomparexyz,
      coordrefxyz,
      searchtype = 'radius',
      radius = userad,
      k = ksuggest
    )
    ignore = tempmatch[[1]] == 0
    tempmatch[[2]][ignore] = NA
    if(smallapprox==FALSE){
      tempmatch[[2]] = 2 * asin(tempmatch[[2]] / 2)
    }
    if(ignoreinternal){
      remove = which(tempmatch[[1]] - 1:length(coordcomparexyz[, 1]) == 0)
      tempmatch[[1]][remove] = 0
      tempmatch[[2]][remove] = NA
    }
    if(ignoreexact){
      remove = which(!(tempmatch[[2]] <= rad & tempmatch[[2]] > 0))
      tempmatch[[1]][remove] = 0
      tempmatch[[2]][remove] = NA
    }else{
      remove = which(!tempmatch[[2]] <= rad)
      tempmatch[[1]][remove] = 0
      tempmatch[[2]][remove] = NA
    }
    if (all(is.na(tempmatch[[2]][,ksuggest]))){
      kendmin = NA
    }else{
      kendmin = min(tempmatch[[2]][, ksuggest], na.rm = TRUE)
    }
    if(is.na(kendmin)==FALSE & ksuggest<kmax){
      comp = tempmatch[[2]][, ksuggest]
      ksuggest = ceiling(ksuggest * max(rad[comp > 0] / comp[comp > 0], na.rm =
                                          TRUE))
      ksuggest = min(ksuggest, kmax)
    }else{
      ksuggest = NA
    }
  }

  keepcols = which(colSums(is.na(tempmatch[[2]])) < N)
  tempmatch[[1]] = matrix(tempmatch[[1]][, keepcols], nrow = N, byrow =
                            FALSE)
  tempmatch[[2]] = matrix(tempmatch[[2]][, keepcols], nrow = N, byrow =
                            FALSE)
  
  if (radunit == "asec"){
    tempmatch[[2]] = tempmatch[[2]]/((pi/180)/3600)
  }
  if (radunit == "amin"){
    tempmatch[[2]] = tempmatch[[2]]/((pi/180)/60)
  }
  if (radunit == "deg"){
    tempmatch[[2]] = tempmatch[[2]]/((pi/180))
  }
  
  if(length(tempmatch[[1]])>0){
    Nmatch = rowSums(tempmatch[[1]]!=0)
  }else{
    Nmatch = NA
  }
  
  if(is.na(Nmatch[1])==FALSE){
    bestID = tempmatch[[1]][, 1]
    bestsep = tempmatch[[2]][, 1]
    select = which(bestID > 0)
    bestsep = bestsep[select]
    bestID = bestID[select]
    orderID = order(bestsep)
    keep = !duplicated(bestID[orderID])
    bestrefID = select[orderID][keep]
    bestcompareID = bestID[orderID][keep]
    bestsep = bestsep[orderID][keep]
    reorderID = order(bestrefID)
    bestrefID = bestrefID[reorderID]
    bestcompareID = bestcompareID[reorderID]
    bestsep = bestsep[reorderID]
    bestmatch = data.frame(refID = bestrefID,
                           compareID = bestcompareID,
                           sep = bestsep)
    output = list(
      ID = tempmatch[[1]],
      sep = tempmatch[[2]],
      Nmatch = Nmatch,
      bestmatch = bestmatch
    )
  } else{
    output = list(
      ID = tempmatch[[1]],
      sep = tempmatch[[2]],
      Nmatch = Nmatch,
      bestmatch = NA
    )
  }
  return(invisible(output))
}

coordmatchsing=function(RAref,Decref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg",
                        radunit='asec', sep = ":", ignoreexact=FALSE, smallapprox=FALSE){
  if (inunitref %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitref must be one of deg, rad or sex")
  }
  if (inunitcompare %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitcompare must be one of deg, rad or sex")
  }
  coordref = rbind(c(RAref, Decref))
  coordcompare = rbind(coordcompare)
  N = length(coordref[, 1])
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
  
  coordrefxyz = sph2car(coordref, deg = TRUE)
  coordcomparexyz = sph2car(coordcompare, deg = TRUE)
  
  dotprod = coordcomparexyz[, 1] * coordrefxyz[1] + coordcomparexyz[, 2] *
    coordrefxyz[2] + coordcomparexyz[, 3] * coordrefxyz[3]
  dotprod[dotprod < -1] = -1
  dotprod[dotprod > 1] = 1
  if(smallapprox==FALSE){
    ang = acos(dotprod)
  }else{
    ang = pi / 2 - dotprod
  }
  if (radunit == "asec"){
    ang = ang / ((pi / 180) / 3600)
  }
  if (radunit == "amin"){
    ang = ang / ((pi / 180) / 60)
  }
  if (radunit == "deg"){
    ang = ang / (pi / 180)
  }
  if (ignoreexact) {
    select = which(ang <= rad & ang > 0)
  } else{
    select = which(ang <= rad)
  }
  if(length(select)==0){
    output = list(
      ID = NA,
      sep = NA,
      Nmatch = NA,
      bestmatch = NA
    )
  }else{
    ID = select
    sep = ang[select]
    orderID = order(sep)
    ID = ID[orderID]
    sep = sep[orderID]
    bestmatch = c(compareID = select[which.min(ang[select])], sep = min(ang[select]))
    output = list(
      ID = ID,
      sep = sep,
      Nmatch = length(select),
      bestmatch = bestmatch
    )
  }
  return(invisible(output))
}

internalclean = function(RA, Dec, rad=2, tiebreak, decreasing = FALSE, Nmatch='all', iter=FALSE, group=FALSE, ...){
  
  if (length(dim(RA)) == 2){
    RA = as.matrix(RA)
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
  
  bestfunc = function(x){
    if(all(x==0)){
      return(0)
    }else{
      return(min(x[x>0],na.rm = TRUE))
    }
  }
  
  matchorder = order(tiebreak, decreasing=decreasing)
  match = coordmatch(
    coordref =  cbind(RA[matchorder], Dec[matchorder]),
    rad = rad,
    ...
  )
  
  if(length(match$ID) == 0){
    return(matchorder)
  }
  
  if(group){
    match$ID[match$ID == 0L] = NA
    links = cbind(1:dim(match$ID)[1], as.integer(match$ID))
    groupinfo = group_links(links, selfgroup=TRUE, return_groupinfo=TRUE)$groupinfo
    return(matchorder[groupinfo[,'IDmin']])
  }
  
  nearcen = apply(cbind(match$bestmatch[, 1], match$ID[match$bestmatch[, 1], ]), 1, bestfunc)
  
  if(Nmatch[1]=='all'){
    output = matchorder[c(which(match$Nmatch == 0), nearcen)]
  }else{
    output = matchorder[nearcen[match$Nmatch[match$bestmatch$ref] %in% Nmatch]]
  }
  
  keep = sort(unique(output))
  
  if(iter){
    new_keep = TRUE
    while(length(new_keep) < length(keep)){
      keep = keep[new_keep] #This should keep the correct relative IDs from the first clean
      new_keep = internalclean(RA=RA[keep],
                               Dec=Dec[keep],
                               rad=rad,
                               tiebreak=tiebreak[keep],
                               decreasing=decreasing,
                               Nmatch=Nmatch,
                               iter=FALSE,
                               ...)
    }
  }
  
  return(keep)
}

group_links = function(links, grouptype='list', selfgroup=FALSE, return_groupinfo=FALSE, return_linkslist=FALSE){
  
  links = links[!is.na(links[,1]) & !is.na(links[,2]),, drop=FALSE]
  
  if(selfgroup==FALSE){
    links = links[which(links[,1] != links[,2]),, drop=FALSE]
  }
  
  links = links[order(links[,1]),]
  links_unique = sort(unique(as.vector(links)))
  
  if(links_unique[length(links_unique)] != length(links_unique)){
    remap = TRUE
    mapping = as.matrix(cbind(1:length(links_unique), links_unique))
    colnames(mapping) = NULL
    links = as.matrix(cbind(mapping[match(links[,1], mapping[,2]),1], mapping[match(links[,2], mapping[,2]),1]))
  }else{
    remap = FALSE
  }
  
  names(links) = NULL
  
  table_init = tabulate(links)
  
  grouplist = list()
  if(return_linkslist){
    linkslist = list()
  }
  
  while(length(links) > 0L){
    Ngroup = 0L
    baseref = links[1,1]
    matchref = links[1,2]
    
    if(selfgroup){
      if(table_init[baseref] == 2 & baseref == matchref){
        #This means we are self matching  
        group = baseref
        Ngroup = 1L
      }
    }
    
    if(Ngroup == 0L){
      if(table_init[baseref] == 1L & table_init[matchref] == 1L){
        #This means we have a pure matching pare
        group = c(baseref, matchref)
        Ngroup = 2L
      }else{
        #Otherwise we have a more complicated group
        group = unique(c(baseref, links[links[,1] %in% baseref,2], links[links[,2] %in% baseref,1]))
        Ngroup = length(group)
        Ngroup_old = 0L
        while(Ngroup_old < Ngroup){
          Ngroup_old = Ngroup
          group = unique(c(group, links[which(links[,1] %in% group),2], links[which(links[,2] %in% group),1]))
          Ngroup = length(group)
        }
      }
    }
    
    names(group) = NULL
    grouplist = c(grouplist, list(sort(group)))
    
    if(return_linkslist){
      linkslist = c(linkslist, list(links[links[,1] %in% group,, drop=FALSE]))
    }
    
    if(length(grouplist) > 0){
      links = links[!links[,1] %in% group,, drop=FALSE]
    }
  }
  
  if(remap){
    for(i in 1:length(grouplist)){
      grouplist[[i]] = mapping[grouplist[[i]],2]
      
      if(return_linkslist){
        linkslist[[i]][] = mapping[linkslist[[i]],2]
      }
    }
  }
  
  if(return_groupinfo){
    Ngroup = sapply(grouplist, length)
    IDmin = sapply(grouplist, min)
    IDmax = sapply(grouplist, max)
    
    if(return_linkslist){
      Nlinks = sapply(linkslist, nrow)
      groupinfo = data.frame(Ngroup=Ngroup, Nlinks=Nlinks, IDmin=IDmin, IDmax=IDmax)
    }else{
      groupinfo = data.frame(Ngroup=Ngroup, IDmin=IDmin, IDmax=IDmax)
    }
  }else{
    Ngroup = NULL
    groupinfo = NULL
  }
  
  if(grouptype == 'list'){
    #Do nothing
    group = grouplist
  }else if(grouptype == 'DF'){
    if(is.null(Ngroup)){
      Ngroup = sapply(grouplist, length)
    }
    
    groupID = rep(1:length(grouplist), Ngroup)
    
    group = data.frame(linkID = unlist(grouplist), groupID=groupID)
  }else{
    stop('grouptype must be one of list or DF!')
  }
  
  if(return_groupinfo){
    if(return_linkslist){
      return(list(group=group, groupinfo=groupinfo, linkslist=linkslist))
    }else{
      return(list(group=group, groupinfo=groupinfo))
    }
  }else{
    if(return_linkslist){
      return(list(group=group, linkslist=linkslist))
    }else{
      return(group)
    }
  }
}
