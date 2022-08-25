deg2dms = function(deg, type='mat', sep=':', digits=2){
if(any(deg< -90 | deg>90)){stop('All deg values should be -90<=deg<=90')}
    deg_sign = sign(deg)
    deg = abs(deg)
    DEG = floor(deg)
    MIN = floor((deg - DEG) * 60)
    SEC = (deg - DEG - MIN/60) * 3600
    SEC=round(SEC,digits)
    SEC[SEC<0]=0; SEC[SEC>60]=60
    MIN[SEC == 60] = MIN[SEC == 60] + 1
    SEC[SEC == 60] = 0
    DEG[MIN == 60] = DEG[MIN == 60] + 1
    MIN[MIN == 60] = 0
    if(digits<0){width=9}
    if(digits==0){width=2}
    if(digits>0){width=digits+3}
    SEC = abs(SEC)
    MIN = abs(MIN)
    DEG = abs(DEG)
    SEC = formatC(SEC, format = "f", width = width, flag = 0, digits = digits)
    MIN = formatC(MIN, format = "f", width = 2, flag = 0, digits = 0)
    DEG = formatC(DEG, format = "f", width = 2, flag = 0, digits = 0)
    DEG[deg_sign == -1] = paste0("-", DEG[deg_sign == -1])
    DEG[deg_sign == 1] = paste0("+", DEG[deg_sign == 1])
    DEG[deg_sign == 0] = paste0("+", DEG[deg_sign == 0])
    if(type=='mat'){output = cbind(DEG, MIN, SEC)}
    if(type=='cat' & sep!='DMS' & sep!='dms'){output=apply(cbind(DEG, MIN, SEC),1,paste,collapse=sep)}
    if(type=='cat' & sep=='DMS'){output=paste(DEG,'D',MIN,'M',SEC,'S',sep='')}
    if(type=='cat' & sep=='dms'){output=paste(DEG,'d',MIN,'m',SEC,'s',sep='')}
    if(type=='cat' & sep=='symbol'){
      output={}
      for(i in 1:length(DEG)){
        output=c(output, parse(text=paste(DEG[i],'*degree*',MIN[i],'*m*',SEC[i],'*s', sep='')))
      }
    }
return(output)
}
