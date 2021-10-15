options(stringsAsFactors=FALSE)
suppressMessages(library(reshape2))
suppressMessages(library(sqldf))
suppressMessages(library(plyr))
suppressMessages(library(minpack.lm))
if(interactive()) {
  setwd('~/d/sci/src/cns_prp_quant/')
}

source('src/shared_functions.R')

# set up for 4-parameter logistic curve fitting
# nice explanation here: https://www.myassays.com/four-parameter-logistic-regression.html
# the following are initial values; the nls() function will iterate until it finds the best fit values for these 4 parameters
a_init = 0 # a = the minimum value that can be obtained (i.e. what happens at 0 dose)
b_init = 1 # b = Hill’s slope of the curve (i.e. this is related to the steepness of the curve at point c).
c_init = 500 # c = the point of inflection (i.e. the point on the S shaped curve halfway between a and d)
d_init = 3 # d = the maximum value that can be obtained (i.e. what happens at infinite dose)

calculate_conc = function(abs, fit) {
  a = coefficients(fit)["a"]
  b = coefficients(fit)["b"]
  c = coefficients(fit)["c"]
  d = coefficients(fit)["d"]
  conc = c * ((a-d)/(abs-d) - 1)^(1/b)
  return(conc)
}

stdcol = '#9D1309'
concol = '#444456'
csfcol = '#3333FF'
spkcol = '#4CBB17'
lincol = '#777777'

hucol = '#239ACD'
ocol = '#28170D'

espresso = '#28170D'

scurve = data.frame(nominal=c(5/(2.5^(0:6)),0),name=paste0('Std0',1:8))

plates = read.table('data/meta/plates.tsv',sep='\t',header=TRUE,comment.char='#')

cat(paste0("\n"))

for (plate in plates$plateno) {
  
  cat(paste0("\rNow processing plate ",plate,"..."))
  flush.console()
  
  datafile = plates$datafile[plates$plateno==plate]
  metafile = plates$metafile[plates$plateno==plate]
  
  datafile_fullpath = paste('data/raw/',datafile,sep='')
  
  if(grepl('\\.txt',datafile,ignore.case=TRUE)) {
    rawdata = process_spectramax_data(datafile_fullpath)
  } else if (grepl('\\.csv',datafile,ignore.case=TRUE) & any(grepl('EnSpire',readLines(datafile_fullpath)))) {
    rawdata = process_enspire_data(datafile_fullpath)
  } else if (grepl('\\.csv',datafile,ignore.case=TRUE)) {
    rawdata = process_optima_data(datafile_fullpath)
  } else {
    stop('unknown file extension: ',datafile)
  }
  
  meta = read.table(paste('data/meta/',metafile,sep=''),sep='\t',header=T,quote='',comment.char='')
  meta = meta[,!grepl('^X',colnames(meta))] # remove extra blank cols
  
  rawdata$plate = plate
  meta$plate = plate
  
  if (castable_to_integer(plate)) {
    plate_prefix = formatC(as.integer(plate),width=3,flag='0')
  } else {
    plate_prefix = plate
  }
  
  data = merge(rawdata,meta,by=c("plate","row","col"),sort=FALSE)
  # merge's sort option sorts as if they were characters, so 1, 10, 11, 12, 2, etc.
  # so I set sort=FALSE in the merge above and then sort separately:
  data = data[order(data$row, data$col),]
  
  if (any(is.na(data$dilution))) {
    data$dilution[is.na(data$dilution)] = 1
  }
  
  #### QC plots
  
  subtitle = paste('cross-species ELISA plate ',plate_prefix,'',sep='')
  
  # plate heatmap
  mat = daply(rawdata, .(row, col), function(x) x$a450_bck)
  png(paste('data/qc/',plate_prefix,'_heatmap.png',sep=''),width=600,height=400,res=100)
  image(x=1:12, y=1:8, z=t(apply(mat, 2, rev)), col=colorRampPalette(c("black", "red"))(10),axes=FALSE, ann=FALSE)
  axis(side=1, at=(1:12), labels=1:12, lwd=0, lwd.ticks=0)
  axis(side=2, at=(1:8), labels=LETTERS[8:1], lwd=0, lwd.ticks=0, las=2)
  title(main='Heatmap')
  mtext(side=3,line=0.5,text=subtitle,cex=.7)
  dev.off()
  
  data$above_baseline = TRUE # set presumption that all data points are above negative controls - will set to FALSE as needed in below loop
  
  data = data[data$stype != '',] # remove unused wells
  
  if (!('standard_curve' %in% colnames(data))) { 
    data$standard_curve = 'only'
  }
  
  for (stdcurvename in unique(data$standard_curve)) {
    this_scurve_rows = data$standard_curve == stdcurvename
    
    # flag wells with absorbance less than blank wells (presumptive pipetting or other errors)
    negative_control_rows = data$stype=='standard' & (data$detail == 'Std08' | data$detail == '0')
    data$above_baseline[this_scurve_rows & data$a450_bck < mean(data$a450_bck[negative_control_rows])] = FALSE
    
    # fill in std curve nominal concs
    data$nominal = as.numeric(NA)
    if ('Std01' %in% data$detail[this_scurve_rows]) {
      data$nominal[data$stype=='standard' & this_scurve_rows] = scurve$nominal[match(data$detail[this_scurve_rows & data$stype=='standard'], scurve$name)]
    } else {
      data$nominal[data$stype=='standard' & this_scurve_rows] = as.numeric(data$detail[this_scurve_rows & data$stype=='standard'])
    }
    
    # per FDA 2018 guidance, anchor points should not be included in curve fit
    fitdata = data[this_scurve_rows & data$stype=='standard' & data$nominal > 0.02048,]
    
    # fit standard curve - use standards from this standard curve, except the zero and any readings below baseline:
    fit = suppressWarnings(nlsLM(a450_bck ~ d + (a - d) / (1 + (nominal/c)^b), 
                start=list(a=a_init,b=b_init,c=c_init,d=d_init), 
                data=fitdata))
    
    # apply the model
    x = (1:400)/100
    f_of_x = calculate_conc(x, fit)
    
    # plot standard curve and its fit
    png(paste('data/qc/',plate_prefix,'_standards_',gsub('[^a-z0-9_]','_',tolower(stdcurvename)),'.png',sep=''),width=600,height=400,res=100)
    par(mar=c(4,4,4,8))
    plot(NA, NA, xlim=c(0,5.1), ylim=c(0,3.5), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
    axis(side=1, at=scurve$nominal, labels=formatC(scurve$nominal,format='fg',digits=2), lwd.ticks=1, lwd=1, cex.axis=.8)
    axis(side=2, at=0:4, lwd.ticks=1, lwd=1, las=2)
    abline(h=0:3,lwd=0.25)
    mtext(side=1, line=2.5, text='nominal concentration (ng/mL)')
    mtext(side=2, line=2.5, text='A450 - A620')
    title(main=paste0(subtitle,' standard curve: ',stdcurvename),cex=0.7)
    # plot the fit model
    points(x = f_of_x, y = x, type='l', lwd=.5)
    # plot the actual points
    points(as.numeric(data$nominal[this_scurve_rows]),data$a450_bck[this_scurve_rows], pch=20, col=stdcol)
    dev.off()
    
    # calculate the inferred concentration for each sample
    data$ngml[this_scurve_rows] = calculate_conc(data$a450_bck[this_scurve_rows], fit)*data$dilution[this_scurve_rows]
    
    # occasionally because of Abck normalization the A450_bck may be <0, in which case ngml is NaN.
    # set these ngml values to 0.
    data$ngml[this_scurve_rows & is.na(data$ngml) & data$a450_bck < 0] = 0.0
    # flag the below/above LLQ cases
    
    # handle the Cambridge Biomedical plates with the earlier std curve configuration:
    if (!all(scurve$nominal %in% data$detail)) {
      llq = min(data$nominal[data$nominal > 0], na.rm=T)
      ulq = max(data$nominal, na.rm=T)
    } else {
      llq = scurve$nominal[scurve$name=='Std06']
      ulq = scurve$nominal[scurve$name=='Std01']
    }
    llq_fluor = mean(data$a450_bck[this_scurve_rows & data$stype=='standard' & data$nominal > 0 & data$nominal <= llq]) 
    ulq_fluor = mean(data$a450_bck[this_scurve_rows & data$stype=='standard' & data$nominal >= ulq]) 
    data$ngml_trunc[this_scurve_rows] = data$ngml[this_scurve_rows] # make a copy of the ngml column, but with range only in [LLQ, ULQ]
    data$flag[this_scurve_rows] = '' # default flag is none
    llq_to_fix = this_scurve_rows & (data$a450_bck < llq_fluor)
    ulq_to_fix = this_scurve_rows & (data$a450_bck > ulq_fluor)
    data$ngml_trunc[llq_to_fix] = llq*data$dilution[llq_to_fix]
    data$flag[llq_to_fix] = 'LLQ'
    data$ngml_trunc[ulq_to_fix] = ulq*data$dilution[ulq_to_fix]
    data$flag[ulq_to_fix] = 'ULQ'
  }
  
  data$ngml = round(data$ngml, 3)
  data$ngml_trunc = round(data$ngml_trunc, 3)
  
  output_cols = c('plate','row','col','a450_bck','stype','detail','dilution','standard_curve','above_baseline','nominal','ngml','ngml_trunc','flag')
  
  # table with well-level detail
  write.table(data[,output_cols],paste('data/processed/',plate_prefix,'.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  
  # summary table with sample-level detail. average for each sample
  # - throw out observations outside dynamic range of assay
  # - if all observations outside dynamic range, set value to limit but also set flag to indicate LLQ/ULQ violated
  smry = sqldf("select detail sample, standard_curve from data where stype = 'sample' group by 1 order by 1;")
  colnames(smry) = c('sample','standard_curve')
  if (nrow(smry) > 0) {
    smry$ngml_av = as.numeric(NA)
    smry$se_mean = as.numeric(NA)
    smry$flag = ''
    for (i in 1:nrow(smry)) {
      relevant_rows = data$stype == 'sample' & data$detail == smry$sample[i] & data$standard_curve == smry$standard_curve[i]
      in_range = relevant_rows & data$flag == ''
      if (sum(in_range) > 0) {
        # handle a special case first: if <100% of replicates are in dynamic range even at the one valid dilution,
        # then base the estimate on averaging the measured value and the LLQ or ULQ as appropriate.
        dilutions_in_range = unique(data$dilution[in_range])
        if (length(dilutions_in_range)==1 & sum(in_range) < sum(relevant_rows)) {
          at_relevant_dilution = relevant_rows & data$dilution %in% dilutions_in_range
          smry$ngml_av[i] = mean(data$ngml_trunc[at_relevant_dilution])
          smry$se_mean[i] = sd(data$ngml_trunc[at_relevant_dilution]) / sqrt(sum(at_relevant_dilution))
        } else { # most of the time, just average all the observations that do fall within dynamic range
          smry$ngml_av[i] = mean(data$ngml[in_range])
          smry$se_mean[i] = sd(data$ngml[in_range]) / sqrt(sum(in_range))
        }
      } else if (all(data$a450_bck[relevant_rows] > ulq_fluor)) {
        smry$ngml_av[i] = ulq * max(data$dilution[relevant_rows])
        smry$se_mean[i] = 0
        smry$flag[i] = 'ULQ'
      } else if (all(data$a450_bck[relevant_rows] < llq_fluor)) {
        smry$ngml_av[i] = llq * min(data$dilution[relevant_rows])
        smry$se_mean[i] = 0
        smry$flag[i] = 'LLQ'
      }
    }
    smry$ngml_av = round(smry$ngml_av, 3)
    smry$se_mean = round(smry$se_mean, 3)
  }
  
  write.table(smry,paste('data/processed/',plate_prefix,'_summary.tsv',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  
  flush.console()
  
}

cat(paste0("\n"))