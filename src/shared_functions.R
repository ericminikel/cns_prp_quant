options(stringsAsFactors=FALSE)
library(reshape2)
library(sqldf)
library(plyr)
if(interactive()) {
  setwd('~/d/sci/src/cns_prp_quant/')
}

expand_range = function(x, by=.5) {
  return ( c(min(x)-by,max(x)+by) )
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}


# take a path to a raw datafile (.txt) from the spectramax platereader, and return a data frame
process_spectramax_data = function(path) {
  file_lines = readLines(path,skipNul=T) # after 2019 May software update, Spectramax now adds nulls on txt export so skipNul=T is needed
  lines_to_skip = suppressWarnings(grep('Plate:',file_lines) + 1)
  rawdata = read.table(textConnection(file_lines),sep='\t',header=F,skip=lines_to_skip,fill=T,skipNul=T) # again, skipNul is now needed
  if (ncol(rawdata) == 16) {
    rawdata = as.matrix(rawdata[1:8,3:14])
    colnames(rawdata) = 1:12
    rownames(rawdata) = LETTERS[1:8]
    rawdata = data.frame(melt(rawdata))
    wavelength = strsplit(suppressWarnings(readLines(path))[lines_to_skip-1],"\t")[[1]][16]
    a_wavelength = paste('a',wavelength,sep='')
    colnames(rawdata) = c('row','col',a_wavelength)
    rawdata$row = as.character(rawdata$row)
  } else if (ncol(rawdata) == 29) {
    rawdata = rawdata[1:8,c(3:14,16:27)]
    abs450 = as.matrix(rawdata[1:8,1:12])
    absbck = as.matrix(rawdata[1:8,13:24])
    abs450_bck = abs450 - absbck
    colnames(abs450_bck) = 1:12
    rownames(abs450_bck) = LETTERS[1:8]
    rawdata = data.frame(melt(abs450_bck))
    colnames(rawdata) = c('row','col','a450_bck')
    rawdata$row = as.character(rawdata$row)
  }
  return (rawdata)
}

# same but for the Fluostar Optima platereader
process_optima_data = function(path) {
  lines_to_skip = grep("Well Row,Well Col",readLines(path))-1 # figure out where the data actually starts
  rawdata = read.table(path,skip=lines_to_skip,header=TRUE,sep=',') # read it in
  rawdata = rawdata[,c(-3,-(ncol(rawdata)))] # remove "Content" column and empty column at end
  colnames(rawdata)[c(1,2)] = c('row','col')
  for (i in 3:ncol(rawdata)) {
    colnames(rawdata)[i] = paste('a',gsub('[A-Za-z\\.]+','',gsub('BP.*','',colnames(rawdata)[i])),sep='')
  }
  if ('a450' %in% colnames(rawdata) & any(grepl('^a6',colnames(rawdata)))) {
    bck_col = grep('^a6',colnames(rawdata))
    bck_colname = colnames(rawdata)[bck_col]
    rawdata$a450_bck = rawdata$a450 - rawdata[,bck_col]
    rawdata$a450_bck[rawdata$a450_bck < 0] = 0
  }
  return (rawdata)
}

process_enspire_data = function(path) {
  # for now, assume A is always 450 and B is 630. In fact, they are specified much further below in the file.
  
  # measurement A
  lines_to_skip_a = grep('Results for Meas A',readLines(path)) # figure out where the data actually starts
  rawdata_a = read.table(path,skip=lines_to_skip_a,header=TRUE,sep=',',nrows = 8) # read it in
  meltdata_a = melt(rawdata_a[,1:13], id.vars = "X") # remove empty col at end
  colnames(meltdata_a) = c('row','col','a450')
  meltdata_a$col = as.integer(gsub('X','',meltdata_a$col))
  
  # measurement B
  lines_to_skip_b = grep('Results for Meas B',readLines(path)) # figure out where the data actually starts
  rawdata_b = read.table(path,skip=lines_to_skip_b,header=TRUE,sep=',',nrows = 8) # read it in
  meltdata_b = melt(rawdata_b[,1:13], id.vars = "X") # remove empty col at end
  colnames(meltdata_b) = c('row','col','a630')
  meltdata_b$col = as.integer(gsub('X','',meltdata_b$col))
  
  stopifnot(all(meltdata_a$row == meltdata_b$row))
  stopifnot(all(meltdata_a$col == meltdata_b$col))
  meltdata_a$a630 = meltdata_b$a630
  meltdata_a$a450_bck = meltdata_a$a450 - meltdata_a$a630

  return (meltdata_a[,c('row','col','a450_bck')])  
}

castable_to_integer = function(x) {
  tryCatch ({
    as.integer(x)
    return (TRUE)
  }, warning = function(w) {
    return (FALSE)
  })
}

# 
# # test
# filename = 'totprot4-2016-11-17.CSV'
# filename = 'SV-2016-09-30-123319-CSF-tp.txt'
# full_path = paste('data/dc/raw/',filename,sep='')
# path = full_path
# 
# filename = 'csf-elisa-08-30-16-164815.txt'
# full_path = paste('data/elisa/raw/',filename,sep='')
# path = full_path
# 


stdcol = '#9D1309'
concol = '#444456'
csfcol = '#3333FF'
spkcol = '#4CBB17'
lincol = '#777777'

hucol = '#239ACD'
ocol = '#28170D'

espresso = '#28170D'
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

process_one_cselisa_plate = function(rawdatatext, metadatatext, plate) {
  
  plate = as.integer(plate)
  
  scurve = data.frame(nominal=c(5/(2.5^(0:6)),0),name=paste0('Std0',1:8))
  
  # to do: make this able to handle Optima data in addition to Spectramax
  rawdata = process_spectramax_data(textConnection(rawdatatext))
  
  meta = read.table(textConnection(metadatatext),sep='\t',header=T,quote='',comment.char='')
  
  if (castable_to_integer(plate)) {
    plate_prefix = formatC(as.integer(plate),width=3,flag='0')
  } else {
    plate_prefix = plate
  }
  
  rawdata$plate = plate
  meta$plate = plate
  
  data = merge(rawdata,meta,by=c("plate","row","col"),sort=FALSE)
  # merge's sort option sorts as if they were characters, so 1, 10, 11, 12, 2, etc.
  # so I set sort=FALSE in the merge above and then sort separately:
  data = data[order(data$row, data$col),]
  
  ## error check the metadata file
  sample_stdcurvenames = unique(data$standard_curve[data$stype=='sample'])
  standard_stdcurvenames = unique(data$standard_curve[data$stype=='standard'])
  orphan_stdcurvenames = setdiff(sample_stdcurvenames, standard_stdcurvenames)
  if (length(orphan_stdcurvenames) > 0) {
    stop(paste0('the following standard curve names appear for samples but not standards, in your metadata file: ',paste(orphan_stdcurvenames,collapse=' ')))
  }
  
  if (any(is.na(data$dilution))) {
    stop(paste0('some wells do not have dilution filled in: ',paste(data$detail[is.na(data$dilution)],collapse=' ')))
  }
  
  if (!('standard_curve') %in% colnames(data)) {
    data$standard_curve = ''
  }
  
  
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
  
  data$above_baseline = TRUE # set presumption that all data points are above negative controls - will set to FALSE if needed later
  
  data = data[data$stype != '',] # remove unused wells
  
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
    
    # for nlsLM 4-param curve fit, we do not need to throw out the zero
    fitdata = data[this_scurve_rows & data$stype=='standard',]
    
    # fit standard curve - use standards from this standard curve, except the zero and any readings below baseline:
    fit = nlsLM(a450_bck ~ d + (a - d) / (1 + (nominal/c)^b), 
                start=list(a=a_init,b=b_init,c=c_init,d=d_init), 
                data=fitdata)
    
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
    llq = scurve$nominal[scurve$name=='Std06']
    ulq = scurve$nominal[scurve$name=='Std01']
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
  }
  
  return(list(data, smry))
}

