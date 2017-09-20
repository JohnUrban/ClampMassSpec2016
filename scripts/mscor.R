## Usage:
## Rscript mscor.R combined.ranks.txt

############ GET ARGS ############
args <- commandArgs(trailingOnly=TRUE)
msg <- "Usage: Rscript mscor.R combined.ranks.txt [optional: outprefix]"
if(length(args) == 0 || length(args) > 2 || args[1] == "-h" || args[1] == "--help"){
  cat(msg,"\n")
  quit()
}
## does file exist?
if(!file.exists(args[1])){
  cat("Error: Arg[1] File does not exist. Did you specify correct path to it?\n")
  cat("Your current working dir is:\n")
  cat("\t",getwd(),"\n")
  cat("Are you in correct directory?\n")
  cat(msg,"\n")
  quit()
}
## if outprefix not given, makes args[2] == ""
if (length(args) == 1){args <- c(args, "")}

# print(getwd())
# print(args[1])
############ LIBRARIES ############
library(lattice)


############ FUNCTIONS ############
corplot <- function(cors, titlepre="", labels=c("None", "AA\nLength", "Mol\nWt", "Num\nPeps")){
  levelplot(t(cors), pretty=TRUE, col.regions=colorRampPalette(c("blue", "black","red")), at=seq(-1,1,0.05), 
            main=paste0(titlepre, " Correlation"), cex=2, scales=list(cex=1.5,
                                                   y=list(alternating=1, at=1:5, labels=labels), 
                                                   x=list(alternating=1, at=1:5, labels=labels)), 
            xlab="", ylab="",
            panel = function(...){ panel.levelplot(...); panel.abline(h = seq(1.5,4.5,1), col="white"); panel.abline(v = seq(1.5,4.5,1), col="white")})
}

get_df <- function(ranks, getcol="score"){
  none=c(ranks[[getcol]][ranks$file == "Kc_CLAMP_IP"], ranks[[getcol]][ranks$file == "Kc_IgG_IP"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP"], ranks[[getcol]][ranks$file == "S2_IgG_IP"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP"], ranks[[getcol]][ranks$file == "S2_IgG_IP.XLIP"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted"])
  aalength=c(ranks[[getcol]][ranks$file == "Kc_CLAMP_IP.length"], ranks[[getcol]][ranks$file == "Kc_IgG_IP.length"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.length"], ranks[[getcol]][ranks$file == "S2_IgG_IP.length"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP.length"], ranks[[getcol]][ranks$file == "S2_IgG_IP.XLIP.length"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted.length"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted.length"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted.length"])
  molwt=c(ranks[[getcol]][ranks$file == "Kc_CLAMP_IP.molweight"], ranks[[getcol]][ranks$file == "Kc_IgG_IP.molweight"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.molweight"], ranks[[getcol]][ranks$file == "S2_IgG_IP.molweight"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP.molweight"], ranks[[getcol]][ranks$file == "S2_IgG_IP.XLIP.molweight"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted.molweight"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted.molweight"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted.molweight"])
  npep=c(ranks[[getcol]][ranks$file == "Kc_CLAMP_IP.numpeps"], ranks[[getcol]][ranks$file == "Kc_IgG_IP.numpeps"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.numpeps"], ranks[[getcol]][ranks$file == "S2_IgG_IP.numpeps"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP.numpeps"], ranks[[getcol]][ranks$file == "S2_IgG_IP.XLIP.numpeps"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted.numpeps"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted.numpeps"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted.numpeps"])
  cordf <- data.frame(none=none, aalength=aalength, molwt=molwt, npep=npep)
}

get_df_igg <- function(ranks, getcol="score"){
  subtract <- c(ranks[[getcol]][ranks$file == "S2.IgGsubtracted"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted.length"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted.length"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted.length"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted.molweight"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted.molweight"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted.molweight"], ranks[[getcol]][ranks$file == "S2.IgGsubtracted.numpeps"], ranks[[getcol]][ranks$file == "Kc.IgGsubtracted.numpeps"], ranks[[getcol]][ranks$file == "S2.XLIP.IgGsubtracted.numpeps"])
  clamp <- c(ranks[[getcol]][ranks$file == "S2_CLAMP_IP"], ranks[[getcol]][ranks$file == "Kc_CLAMP_IP"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.length"], ranks[[getcol]][ranks$file == "Kc_CLAMP_IP.length"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP.length"])
  clamp <- c(clamp, ranks[[getcol]][ranks$file == "S2_CLAMP_IP.molweight"], ranks[[getcol]][ranks$file == "Kc_CLAMP_IP.molweight"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP.molweight"])
  clamp <- c(clamp, ranks[[getcol]][ranks$file == "S2_CLAMP_IP.numpeps"], ranks[[getcol]][ranks$file == "Kc_CLAMP_IP.numpeps"], ranks[[getcol]][ranks$file == "S2_CLAMP_IP.XLIP.numpeps"])
  cordf <- data.frame(subtract=subtract, clamp=clamp)
}



############ EXECUTE ############
ranks <- read.table(args[1], col.names = c('rank', 'acc', 'score', 'name', 'file'), colClasses = c("numeric", "character", "numeric", "character", "factor"))

cortype <- "scores"
cordf <- get_df(ranks, "score") 
cors <- cor(cordf)
fname <-  paste0(args[2], cortype, ".correlation")
cat("Correlation Matrix for Scores\n")
print(cors)
cat("\n")
write.table(cors, file = paste0(fname,".txt"), sep=",",quote = FALSE, col.names=NA)
pdf(file = paste0(fname,".pdf"))
corplot(cors, "Scores")
garbage <- dev.off()

cortype <- "rank"
cordf <- get_df(ranks, "rank")
cors <- cor(cordf)
fname <- paste0(args[2], cortype, ".correlation")
cat("Correlation Matrix for Rankings\n")
print(cors)
cat("\n")
write.table(cors, file = paste0(fname,".txt"), sep=",",quote = FALSE, col.names=NA)
pdf(file = paste0(fname,".pdf"))
corplot(cors, "Rank")
garbage <- dev.off()


#### Looking at subtracting out IgG vs not doing so #####
cortype <- "score"
cordf <- get_df_igg(ranks, "score") 
cors <- cor(cordf)
fname <-  paste0(args[2], cortype, ".clampVsClamp-IgG.correlation")
cat("Correlation Matrix for Scores: CLAMP vs. CLAMP-IgG\n")
print(cors)
cat("\n")
write.table(cors, file = paste0(fname,".txt"), sep=",",quote = FALSE, col.names=NA)
pdf(file = paste0(fname,".pdf"))
corplot(cors, titlepre="CLAMP vs. CLAMP-IgG: Scores", labels=c("CLAMP-IgG","CLAMP"))
garbage <- dev.off()


cortype <- "rank"
cordf <- get_df_igg(ranks, "rank") 
cors <- cor(cordf)
fname <-  paste0(args[2], cortype, ".clampVsClamp-IgG.correlation")
cat("Correlation Matrix for Rankings: CLAMP vs. CLAMP-IgG\n")
print(cors)
cat("\n")
write.table(cors, file = paste0(fname,".txt"), sep=",",quote = FALSE, col.names=NA)
pdf(file = paste0(fname,".pdf"))
corplot(cors, titlepre="CLAMP vs. CLAMP-IgG: Rank", labels=c("CLAMP-IgG","CLAMP"))
garbage <- dev.off()


