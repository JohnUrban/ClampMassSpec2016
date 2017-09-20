## Usage:
## Rscript mscor.R combined.ranks.txt

############ GET ARGS ############
args <- commandArgs(trailingOnly=TRUE)
msg <- "Usage: Rscript glu.R data prefix"
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
## if outprefix not given, makes args[2] == args[1]
if (length(args) == 1){args <- c(args, args[1])}

# print(getwd())
# print(args[1])
############ LIBRARIES ############


############ FUNCTIONS ############



############ EXECUTE ############
counts <- read.table(args[1], col.names = c('count'), colClasses = c("numeric"))

pdf(file = paste0(args[2],".pdf"))

plot(1:length(counts$count), counts$count, xlab="Position (AA)", ylab="Density", las=1, type="h", col="blue", ylim=c(0,1))

garbage <- dev.off()

