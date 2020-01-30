python.setup <- function(){
  
}


python.start <- function(python="~/../AppData/Local/r-miniconda/python", list.libs=list(pygdal ="gdal", pyosr="osr")){
  pystatus <- FALSE
  use_python(python)
  if(length(list.libs)>0){
    for(i in 1:length(list.libs)){
      assign(names(list.libs)[i], import(list.libs[[i]]), envir = .GlobalEnv)
    }
  }
  return(pystatus)
}
