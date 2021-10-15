options(stringsAsFactors=F)
if (interactive()) {
	setwd('~/d/sci/src/cns_prp_quant/')
}
library(foreign)

# for a while the Optima Data Analysis program on our Fluostar Optima platereader's computer was broken
# therefore we couldn't export results as csv, but had to grab the raw dbf files from BMG's Program Files
# directory and process them through this script

platedate = '2021-10-14'
plateno = '00001'
a450_dbf = '001'
abck_dbf = '002'

dat = read.dbf(paste('data/salvage/',plateno,'/',a450_dbf,'.dbf',sep=''))
dat = dat[3:98,]
dat$well_row = substr(dat$WELLNUM,1,1)
dat$well_col = as.integer(substr(dat$WELLNUM,2,3))
master = dat[,c('well_row','well_col','CONTENT','M1')]
colnames(master)[1:3] = c("Well Row","Well Col","CONTENT")
colnames(master)[4] = 'Raw Data (450BP) 1'

dat = read.dbf(paste('data/salvage/',plateno,'/',abck_dbf,'.dbf',sep=''))
dat = dat[3:98,]
dat$well_row = substr(dat$WELLNUM,1,1)
dat$well_col = as.integer(substr(dat$WELLNUM,2,3))
master$V5 = dat$M1[match(master$CONTENT,dat$CONTENT)]
colnames(master)[5] = 'Raw Data (620BP) 1'

master$V6 = '' # optima output has a blank column at end

write.table(master,paste('data/raw/',ifelse(is.null(platedate),Sys.Date(),platedate),'-',plateno,'-salvage.csv',sep=''),sep=',',row.names=F,col.names=T,quote=F)

