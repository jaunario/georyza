# THIS IS A CONTAINER FOR AUTOMATIC SETUP OF R AND UPDATING OF SCRIPTS

# ENABLE Google Sheet API
# https://developers.google.com/sheets/api/quickstart/python#step_1_turn_on_the

GSHEET_ID = "1FjYVV26sA2JX91N3x5xNucMpmUOqIHIr92meOuBoSig"
TILE_RANGE    = "Tiles!A:G"

DEBUG = length(commandArgs(trailingOnly = TRUE))<=0


library(orysat)
library(jsonlite)

myinfo <- Sys.info()
r.info <- sessionInfo()

if (DEBUG){
  session.args <- c("MODE='DOWNLOAD'", "TILE='h25v07'", "PRODUCT='MOD09A1'", "YEAR=2019", "PROGRESS=1", "ACTION='none'")
} else {
  session.args <- commandArgs(trailingOnly = TRUE)
}

if(sum(grepl("MODE", session.args))==0){
  stop("SESSION Error: MODE not provided.")
}

message(session.args)

eval(parse(text=session.args))

pylibs.start(python="~/../miniconda3", list.libs = 
               list(pickle= "pickle"
                    #,os="os.path"
                    ,gs="googleapiclient.discovery" 
                    #,gs_request="google.auth.transport.requests" 
                    #,gs_flow="google_auth_oauthlib.flow")
                    ))

builtins <- import_builtins()

if (file.exists("token.pickle")){
  token <- builtins$open("token.pickle", "rb")
  cred <- pickle$load(token)
}

if (!exists("cred") || cred$expired){
  if (exists("cred")) rm(token, cred)
  py_run_file("gs_auth.py")
  token <- builtins$open("token.pickle", "rb")
  cred <- pickle$load(token)
} 

service <- gs$build("sheets","v4",credentials=cred)
sheet <- service$spreadsheets()
sv <- sheet$values()

if(MODE=="CHECKIN"){
  SESSION_RANGE = "Checkin!A:H"
  gs_dat <- sv$get(spreadsheetId=GSHEET_ID, range=SESSION_RANGE)$execute()
  SESSION_RANGE <- gs_dat$range
  checkin <- c(myinfo[c("nodename","user","sysname","release", "version")],
                paste0(r.info$R.version$major,".",r.info$R.version$minor), # R version
                format(Sys.Date(), "%m/%d/%Y"), # Date
                format(Sys.Date(), "%m/%d/%Y")) # Time
  
  body_update <- list(values=list(as.vector(checkin)))
  #body_update <- paste0(substr(body_update,1, 14), substr(body_update,17, nchar(body_update)-3), "}")
  sv$append(spreadsheetId=GSHEET_ID, range=SESSION_RANGE, valueInputOption="RAW", body=body_update)$execute()
}

if(MODE=="CHECKOUT"){
  SESSION_RANGE = "Checkout!A:H"
  gs_dat <- sv$get(spreadsheetId=GSHEET_ID, range=SESSION_RANGE)$execute()
  SESSION_RANGE <- gs_dat$range
  checkin <- c(myinfo[c("nodename","user","sysname","release", "version")],
               paste0(r.info$R.version$major,".",r.info$R.version$minor), # R version
               format(Sys.Date(), "%m/%d/%Y"), # Date
               format(Sys.Date(), "%m/%d/%Y")) # Time
  
  body_update <- list(values=list(as.vector(checkin)))
  #body_update <- paste0(substr(body_update,1, 14), substr(body_update,17, nchar(body_update)-3), "}")
  sv$append(spreadsheetId=GSHEET_ID, range=SESSION_RANGE, valueInputOption="RAW", body=body_update)$execute()
}

if(MODE=="DOWNLOAD"){
  SESSION_RANGE = "DOWNLOADS!A:I"
  gs_dat <- sv$get(spreadsheetId=GSHEET_ID, range=SESSION_RANGE)$execute()
  dat.dlstatus <- do.call(rbind, gs_dat$values)
  
  SESSION_RANGE <- gs_dat$range
  
  if(sum(grepl("TILE", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
  if(sum(grepl("PRODUCT", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
  if(sum(grepl("YEAR", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
  if(sum(grepl("PROGRESS", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
  if(sum(grepl("ACTION", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
  if(PROGRESS==0 && ACTION=='wait'){
      body_update <- list(values=list(as.vector(c(TILE, YEAR, PRODUCT, myinfo["nodename"], myinfo["user"], PROGRESS, ACTION, format(Sys.time(), "%b %d, %Y %H:%M %z"), "NA"))))
      sv$append(spreadsheetId=GSHEET_ID, range=SESSION_RANGE, valueInputOption="RAW", body=body_update)$execute()
  } else {
      rr <- which(dat.dlstatus[,1]==TILE & dat.dlstatus[,2]==YEAR & dat.dlstatus[,3]==PRODUCT & dat.dlstatus[,4]==myinfo["nodename"] & dat.dlstatus[,5]==myinfo["user"])
      body_update <- as.vector(dat.dlstatus[rr,])
      body_update[6] <- PROGRESS
      body_update[7] <- ACTION
      body_update[9] <- format(Sys.time(), "%b %d, %Y %H:%M %z")
      body_update <- list(values=list(as.vector(body_update)))  
      SESSION_RANGE <- paste0("Downloads!A",rr,":I",rr)
      sv$update(spreadsheetId=GSHEET_ID, range=SESSION_RANGE, valueInputOption="RAW", body=body_update)$execute()
  }
  
}

if(MODE=="PROGRESS"){
  SESSION_RANGE = "MSF500_palay"
  gs_dat <- sv$get(spreadsheetId=GSHEET_ID, range=SESSION_RANGE)$execute()
  SESSION_RANGE <- gs_dat$range
  
  if(sum(grepl("TILE", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
  if(sum(grepl("YEAR", session.args))==0){
    stop("SESSION Error: MODE not provided.")
  }
  
}
tiles <- sv$get(spreadsheetId=GSHEET_ID, range=TILE_RANGE)
dat.tiles <- tiles$execute()$values
