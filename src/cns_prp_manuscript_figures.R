options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/cns_prp_quant')
}
options(warn=1)

start_time = Sys.time()

cat(file=stderr(), 'Loading packages...')

### PACKAGES

suppressMessages(library(sqldf))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(minpack.lm))
suppressMessages(library(magick))

cat(file=stderr(), 'Done.\nLoading functions and datasets...')

### FUNCTIONS

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

# format percents like a percent instead of proportion. two significant digits by default
percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

signed_percent = function(proportion,digits=2) {
  result = gsub(' ','',paste(ifelse(proportion > 0, '+', ''),round(proportion*100,0),"%",sep="") )
  result[is.na(proportion)] = '-'
  return ( result )
}


format_p = function(p) {
  formatted_p = paste0(' = ',formatC(signif(p,2), format='g', digits=2))
  formatted_p[signif(p,2) < 0.1] = paste0(' = ',formatC(signif(p[signif(p,2) < 0.1],2), format='f', digits=3))
  formatted_p[signif(p,2) < 0.01] = paste0(' = ',formatC(signif(p[signif(p,2) < 0.1],2), format='e', digits=1))
  formatted_p[signif(p,2) == 1] = ' = 1.0'
  formatted_p[signif(p,2) < 1e-10] = ' < 1.0e-10'
  return (formatted_p)
}

prpcol = '#0001CD'

# stuff for 4-point curve fitting:
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

clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

summarize_ci = function(df, bycols, valcol) {
  bycol_commadelim = paste0(bycols, collapse=', ')
  bycol_colnums = paste0(1:length(bycols), collapse=', ')
  query = paste0("
                 select   ",bycol_commadelim,", avg(",valcol,") mean, stdev(",valcol,") sd, count(",valcol,") n
                 from     df
                 where    ",valcol," is not null
                 group by ",bycol_colnums,"
                 order by ",bycol_colnums,"
                 ;")
  smry = sqldf(query)
  smry$l95 = smry$mean - 1.96 * smry$sd / sqrt(smry$n)
  smry$u95 = smry$mean + 1.96 * smry$sd / sqrt(smry$n)
  return (smry)
}


startplot = function(xlims, ylims) {
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
}


species_meta = data.frame(std = c('MoPrP16','RaPrP50','CyPrP51','HuPrP5','BvPrP37','SHaPrP71'),
                          name = c('mouse','rat','cynomolgus','human','bank vole','hamster'),
                          color = c('#d95f02','#1b9e77','#7570b3','#c69b02','#66a61e','#e7298a')) #c('#7fc97f','#beaed4','#675423','#978213','#386cb0','#f0027f')) 



### START AN OUTPUT FILE FOR STATS TO QUOTE IN MANUSCRIPT TEXT

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T

### READ IN DATA

peptides = read.table('data/mrm/ref/peptides.tsv',sep='\t',header=T,quote='',comment.char='')

validation_plates = 12:18
plates = read.table('data/meta/plates.tsv',sep='\t',header=T)
plates$validation = FALSE
plates$validation[plates$plateno %in% validation_plates] = TRUE

for (p in 1:nrow(plates)) {
  plateno = plates$plateno[p]
  plate_prefix = formatC(plates$plateno[p],width=3,flag='0')
  this_elisa_raw = read.table(paste0('data/processed/',plate_prefix,'.tsv'),sep='\t',header=T)
  this_elisa = read.table(paste0('data/processed/',plate_prefix,'_summary.tsv'),sep='\t',header=T)
  this_elisa_raw$plate = plateno
  this_elisa$plate = plateno
  if (p == 1) {
    elisa = this_elisa
    elisa_raw = this_elisa_raw
  } else {
    elisa = rbind(elisa, this_elisa)
    elisa_raw = rbind(elisa_raw, this_elisa_raw)
  }
}

cat(file=stderr(), 'Done.\nCreating Table S1...')

### TABLE S1

calibration_intra = sqldf("
select   detail, plate, count(*) n, avg(a450_bck) mean_abs, stdev(a450_bck)/avg(a450_bck) cv_abs, 
         avg(ngml) mean_backfit, stdev(ngml)/avg(ngml) cv_backfit 
from     elisa_raw
where    stype = 'standard'
and      standard_curve = 'only'
and      plate in (select plateno from plates where validation)
group by 1, 2
having   count(*) = 6
order by 1
;")
calibration_inter = sqldf("
select   detail, count(*) n_plate, avg(mean_abs) mean_abs1, avg(mean_backfit) mean_backfit1, stdev(mean_backfit)/avg(mean_backfit) cv_backfit
from (
select   detail, plate, count(*) n_repl, avg(a450_bck) mean_abs, avg(ngml) mean_backfit
from     elisa_raw
where    stype = 'standard'
and      (standard_curve is null or standard_curve = 'MoPrP16')
and      plate in (select plateno from plates where validation)
group by 1, 2
)
group by 1
order by 1 desc
;")

calibration = calibration_inter
calibration$nominal = formatC(as.numeric(calibration$detail), format='f', digits=2) # round(as.numeric(calibration_intra$detail),2)
calibration$mean_backfit = formatC(calibration$mean_backfit1, format='f', digits=2)
calibration$mean_abs = formatC(calibration$mean_abs1, format='f', digits=3)
calibration$cv_backfit_interplate = percent(calibration$cv_backfit, digits=1)
calibration$cv_backfit_intraplate = percent(calibration_intra$cv_backfit[match(calibration$detail, calibration_intra$detail)], digits=1)
calibration$cv_abs_intraplate = paste0(formatC(calibration_intra$cv_abs[match(calibration$detail, calibration_intra$detail)]*100, format='f', digits=1),'%')
calibration$fold_blank = formatC(calibration$mean_abs1 / calibration$mean_abs1[calibration$detail=='0'],format='f',digits=1)

calibration_disp_cols = c('nominal','mean_abs','cv_abs_intraplate','mean_backfit','cv_backfit_intraplate','cv_backfit_interplate','fold_blank')
# calibration[,calibration_disp_cols]

qc_identities = data.frame(y=4:1,
                           fullname=c('Mo Pos Hi QC','Mo Pos Mid QC','Mo Pos Lo QC','Mo Neg QC'),
                           shortname=c('High QC','Mid QC','Low QC','Neg QC'),
                           identity=c('WT','het KO','90% hom KO / 10% WT','hom KO'),
                           color=c('#4a1486', '#6a51a3', '#807dba', '#bdbdbd'))

qc_intra = sqldf("
                 select   detail, plate, count(*) n, avg(a450_bck) mean_abs, stdev(a450_bck)/avg(a450_bck) cv_abs, 
                 avg(ngml) mean_backfit, stdev(ngml)/avg(ngml) cv_backfit 
                 from     elisa_raw
                 where    stype = 'sample' and detail like '%QC%' and dilution = 200
                 and      plate = 15 -- the validation plate with 6 replicates of each QC
                 group by 1, 2
                 -- having   count(*) >= 5
                 order by 1
                 ;")
qc_inter = sqldf("
                 select   detail, count(*) n_plate, avg(mean_abs) mean_abs1, avg(mean_backfit) mean_backfit1, stdev(mean_backfit)/avg(mean_backfit) cv_backfit
                 from (
                 select   detail, plate, count(*) n_repl, avg(a450_bck) mean_abs, avg(ngml) mean_backfit
                 from     elisa_raw
                 where    stype = 'sample' and detail like '%QC%' and dilution = 200
                 and      standard_curve in ('only','MoPrP16')
                 and      plate in (select plateno from plates where validation)
                 group by 1, 2
                 )
                 group by 1
                 order by 1 desc
                 ;")

qc = qc_inter
qc$qc = qc$detail
qc$mean_backfit = formatC(qc$mean_backfit1, format='f', digits=2)
qc$mean_abs = formatC(qc$mean_abs1, format='f', digits=3)
qc$cv_backfit_interplate = percent(qc$cv_backfit, digits=1)
qc$cv_backfit_intraplate = percent(qc_intra$cv_backfit[match(qc$detail, qc_intra$detail)], digits=1)
qc$cv_abs_intraplate = percent(qc_intra$cv_abs[match(qc$detail, qc_intra$detail)], digits=1)
qc$percent_high = percent(qc$mean_backfit1 / qc$mean_backfit1[qc$detail=='Mo Pos Hi QC'],digits=1)

qc$y = qc_identities$y[match(qc$detail, qc_identities$fullname)]
qc$shortname = qc_identities$shortname[match(qc$detail, qc_identities$fullname)]
qc$identity = qc_identities$identity[match(qc$detail, qc_identities$fullname)]

qc_disp_cols = c('shortname','identity','mean_abs','cv_abs_intraplate','mean_backfit','cv_backfit_intraplate','cv_backfit_interplate','percent_high')

# TABLE S1 is in two halves. suppressWarnings on 2nd half b/c appending colnames to file
write.table(qc[with(qc, order(-y)),qc_disp_cols], 'display_items/table-s1.tsv', sep='\t', row.names=F, quote=F, col.names=T)
suppressWarnings(write.table(calibration[,calibration_disp_cols], 'display_items/table-s1.tsv', sep='\t', row.names=F, quote=F, col.names=T, append=T))








cat(file=stderr(), 'Done.\nCreating Figure S1...')

### FIGURE S1

resx=300
png('display_items/figure-s1.png',width=6.5*resx,height=6*resx,res=resx)

layout_matrix = matrix(c(1,1,6,6,
                         1,1,7,7,
                         2,3,8,8,
                         4,5,8,8), nrow=4, byrow=T)
layout(layout_matrix)

panel = 1

par(mar=c(6,6,3,0.5))
ab_screen = read.table('data/assaydev/opt1_summary.tsv',sep='\t',header=T)
ab_screen_meta = data.frame(ab=unique(ab_screen$capture_ab), x=1:4)
ab_screen$capno = ab_screen_meta$x[match(ab_screen$capture_ab, ab_screen_meta$ab)]
ab_screen$detno = ab_screen_meta$x[match(ab_screen$detect_ab, ab_screen_meta$ab)]
ab_screen$sn = ab_screen$huprp200ngml / ab_screen$huprp_0ngml
ab_screen_img = as.matrix(acast(data=ab_screen, capno ~ detno, value.var='sn'))
image(ab_screen_img, x=ab_screen_meta$x, y=ab_screen_meta$x, axes=F, xlab='', ylab='', col=hcl.colors(16, "YlOrRd", rev = TRUE)[1:12])
mtext(side=1, at=ab_screen_meta$x, text=ab_screen_meta$ab, las=2)
mtext(side=2, at=ab_screen_meta$x, text=ab_screen_meta$ab, las=2)
text(x=ab_screen$capno, y=ab_screen$detno, labels=formatC(ab_screen$sn, format='f', digits=1), font=2, cex=1)
mtext(side=1, line=4, text='capture Ab')
mtext(side=2, line=4, text='detection Ab')
rect(xleft=0.5,xright=1.5,ybottom=2.5,ytop=3.5,border='#000000',lwd = 3)
rect(xleft=0.5,xright=1.5,ybottom=1.5,ytop=2.5,border='#000000',lwd = 3)
rect(xleft=1.5,xright=2.5,ybottom=0.5,ytop=1.5,border='#000000',lwd = 3)
rect(xleft=2.5,xright=3.5,ybottom=0.5,ytop=1.5,border='#000000',lwd = 3)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1


ab_dr = read.table('data/assaydev/opt3_4_summary.tsv',sep='\t',header=T)
ab_dr_meta = sqldf("select capture_ab, detect_ab from ab_dr group by 1, 2 order by 1, 2;")
ab_dr_meta$pair = 1:nrow(ab_dr_meta)
ab_dr_meta$disp = paste0(ab_dr_meta$capture_ab, ' capture\n', ab_dr_meta$detect_ab, ' detect')
xlims = c(0.1, 100)
xbigs = c(.1, 1, 10, 100)
xats = rep(1:9, 4) * rep(c(.1, 1, 10, 100), each=9)
ylims = c(0, 4)

par(mar=c(3,3,2.5,0.5))
for (i in 1:nrow(ab_dr_meta)) {
  subs = ab_dr[ab_dr$capture_ab == ab_dr_meta$capture_ab[i] & ab_dr$detect_ab == ab_dr_meta$detect_ab[i],]
  plot(NA, NA, xlim=xlims, ylim=ylims, log='x', axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
  axis(side=1, at=xats, labels=NA, tck=-0.025)
  axis(side=2, at=0:4, labels=NA, tck=-0.05, las=2)
  axis(side=2, at=0:4, lwd=0, las=2, line=-0.5)
  axis(side=2, at=0:8/2, labels=NA, tck=-0.025)
  mtext(side=1, line=1.25, text='nominal [PrP] (ng/mL)', cex=0.7)
  mtext(side=2, line=1.25, text='absorbance', cex=0.7)
  points(subs$prp_ngml, subs$rat, pch=20, col=species_meta$color[species_meta$name=='rat'])
  fit = suppressWarnings(nlsLM(rat ~ d + (a - d) / (1 + (prp_ngml/c)^b), 
              start=list(a=a_init,b=b_init,c=c_init,d=d_init), 
              data=subs[subs$prp_ngml > 0,]))
  x = 1:400/100
  f_of_x = calculate_conc(x, fit)
  points(x=f_of_x, y=x, type='l', lwd=.5, col=species_meta$color[species_meta$name=='rat'])
  points(subs$prp_ngml, subs$human, pch=20, col=species_meta$color[species_meta$name=='human'])
  fit = suppressWarnings(nlsLM(human ~ d + (a - d) / (1 + (prp_ngml/c)^b), 
              start=list(a=a_init,b=b_init,c=c_init,d=d_init), 
              data=subs))
  x = 1:400/100
  f_of_x = calculate_conc(x, fit)
  points(x=f_of_x, y=x, type='l', lwd=.5, col=species_meta$color[species_meta$name=='human'])
  text(x=subs$prp_ngml[5], y=  subs$rat[5], labels='rat',   col=species_meta$color[species_meta$name=='rat'], pos=2)
  text(x=subs$prp_ngml[5], y=subs$human[5], labels='human', col=species_meta$color[species_meta$name=='human'], pos=4)
  mtext(side=3, line=0, text=ab_dr_meta$disp[i], cex=0.7)
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.2, line = 0.5)
  panel = panel + 1
}

# detergents experiment
data = read.table('data/assaydev/006_summary.tsv',sep='\t',header=T)
data$animal = substr(data$sample,1,3)
data$hemi = substr(data$sample,5,5)
data$spin = substr(data$sample, 7,10)

detergents = data.frame(animal=formatC(1:9,width=3,flag='0'), detergent=rep(c('CHAPS 0.2%','CHAPS 0.03%','RIPA'),each=3), x=rep(1:3, each=3))

deterglabs = sqldf("select detergent, x from detergents group by 1, 2 order by 2, 1;")
deterglabs$disp = paste0(gsub('[A-Z ]*','',deterglabs$detergent),'\n',gsub(' .*','',deterglabs$detergent))

data$deterg = detergents$detergent[match(data$animal, detergents$animal)]
data$x_det = detergents$x[match(data$animal, detergents$animal)]
data$x = data$x_det
data$x[data$spin=='post'] = data$x[data$spin=='post'] + 3
x_offset = 0.25
data$x[data$hemi=='L'] = data$x[data$hemi=='L'] - x_offset
data$x[data$hemi=='R'] = data$x[data$hemi=='R'] + x_offset

step1 = sqldf("
              select   deterg, spin, animal, stdev(ngml_av)/avg(ngml_av) cv, avg(ngml_av) mean
              from     data
              group by 1, 2, 3
              order by 2,1,3
              ;")

step2 = sqldf("
              select   x x_det, deterg, spin, avg(cv) mean_inter_hemisphere_cv, avg(mean) mean_ngml, stdev(mean)/avg(mean) inter_animal_cv
              from     step1 s, detergents d
              where    s.animal = d.animal
              group by 1, 2, 3
              order by 2, 1
              ;")
step2$x = step2$x_det + ifelse(step2$spin=='post',3,0)

norm_value = mean(step2$mean_ngml[step2$deterg=='CHAPS 0.2%' & step2$spin=='pre'])
step2$rel = step2$mean_ngml / norm_value
data$rel = data$ngml_av / norm_value


llq_plate6 = .148*100
llq_rel = llq_plate6 / norm_value

par(mar=c(0.5,6,3,4))
xlims = c(0.5, 6.5)
plot(NA, NA, xlim=xlims, ylim=c(0,0.25), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, at=0:5/20, labels=percent(0:5/20, digits=0), las=2)
barwidth = 0.4
rect(xleft=step2$x - barwidth, xright=step2$x + barwidth, ybottom=rep(0,nrow(step2)), ytop=step2$mean_inter_hemisphere_cv, col='#7C7C7C', border=NA)
mtext(side=2, line=3.0, text='inter-hemisphere\nmean CV', cex=0.8)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

par(mar=c(0.5,6,3,4))
plot(NA, NA, xlim=xlims, ylim=c(0,0.25), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, at=0:5/20, labels=percent(0:5/20, digits=0), las=2)
barwidth = 0.4
rect(xleft=step2$x - barwidth, xright=step2$x + barwidth, ybottom=rep(0,nrow(step2)), ytop=step2$inter_animal_cv, col='#7C7C7C', border=NA)
mtext(side=2, line=3.0, text='inter-animal CV', cex=0.8)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

par(mar=c(6,6,3,4))
plot(NA, NA, xlim=xlims, ylim=c(0,1.25), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
#axis(side=1, line=-1, at=rep(1:6,each=2)+rep(c(-x_offset,x_offset),6), labels=rep(c('L','R'),6), lwd=0, cex.axis=0.5)
mtext(side=1, line=-0.25, at=rep(1:6,each=2)+rep(c(-x_offset,x_offset),6), text=rep(c('L','R'),6), cex=0.4, font=2)
mtext(side=1, line=1.25, at=deterglabs$x,   text=deterglabs$disp, cex=0.5)
mtext(side=1, line=1.25, at=deterglabs$x+3, text=deterglabs$disp, cex=0.5)
axis(side=1, at=c(0.6, 3.4), line=2.75, tck=0.025, labels=NA)
axis(side=1, at=c(3.6, 6.4), line=2.75, tck=0.025, labels=NA)
mtext(side=1, line=2.75, at=2, text='pre-spin', cex=0.8)
mtext(side=1, line=2.75, at=5, text='post-spin', cex=0.8)
axis(side=2, at=0:5/4, labels=percent(0:5/4), las=2)
abline(h=1, lty=3, lwd=0.75)
abline(h=llq_rel, lty=3, lwd=0.75, col='red')
mtext(side=4, at=llq_rel, text='LLQ', col='red', las=2, line=0.25, cex=0.9)
mtext(side=2, line=3.0, text='brain [PrP]\n(normalized)', cex=0.8)
points(data$x, data$rel, pch=20)
for (animal in data$animal) {
  for (spin in c('pre','post')) {
    subs = data[data$animal==animal & data$spin==spin,]
    points(subs$x, subs$rel, type='l', lwd=0.5)
  }
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

silence_message = dev.off()




cat(file=stderr(), 'Done.\nCreating Figure S2...')

### FIGURE S2

resx=300
png('display_items/figure-s2.png',width=6.5*resx,height=5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,
                         4,5,6
                         ), byrow=T, nrow=2)
layout(layout_matrix)

panel = 1

par(mar=c(4,4,3,2))

bh_parallelism = elisa_raw[elisa_raw$plate==12 & elisa_raw$stype=='sample',]
bh_plism_smry = sqldf("select detail, dilution, avg(ngml) mean from bh_parallelism group by 1,2 order by 1;")
bh_parallelism$reference = bh_plism_smry$mean[bh_plism_smry$detail=='Mo Pos Hi QC'][match(bh_parallelism$dilution, bh_plism_smry$dilution[bh_plism_smry$detail=='Mo Pos Hi QC'])]
bh_parallelism$rel = bh_parallelism$ngml / bh_parallelism$reference
bh_parallelism$color = qc_identities$color[match(bh_parallelism$detail, qc_identities$fullname)]
bh_plism_smry = sqldf("select detail, dilution, color, avg(ngml) mean, avg(rel) mean_rel from bh_parallelism group by 1,2,3 order by 1;")

xlims = c(100, 10000)
xats = rep(1:9,3) * rep(c(100,1000,10000), each=9)
xbigs = c(100, 1000, 10000)
ylims = c(0, 400)

plot(NA, NA, xlim=xlims, ylim=ylims, log='x', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=paste0('1:',xbigs), lwd=0, line=-0.5, las=2)
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=2, at=0:4*100, labels=NA, tck=-0.05, las=2)
axis(side=2, at=0:4*100, lwd=0, las=2, line=-0.5)
axis(side=2, at=0:8/2*100, labels=NA, tck=-0.025)
mtext(side=1, line=3, text='dilution', cex=0.8)
mtext(side=2, line=2.5, text='adjusted [PrP] (ng/g)', cex=0.8)
points(x=bh_parallelism$dilution, y=bh_parallelism$ngml, col=bh_parallelism$color, pch=20)
for (smpl in unique(bh_plism_smry$detail)) {
  subs = bh_plism_smry[bh_plism_smry$detail==smpl,]
  points(x=subs$dilution, y=subs$mean, col=subs$color, type='l', lwd=.75)
}
llq = 0.05
ulq = 5
adjusted_llqs = llq * xats 
adjusted_ulqs = ulq * xats
points(x=xats, y=adjusted_llqs, col='red', type='l', lty=3, lwd=1.25)
points(x=xats, y=adjusted_ulqs, col='red', type='l', lty=3, lwd=1.25)
par(xpd=T)
text(x=xats[max(which(xats < xlims[2]))], y=ylims[2], pos=3, col='red', labels='LLQ')
par(xpd=F)
legend('topleft',qc_identities$fullname,lwd=.75,pch=20,col=qc_identities$color,text.col=qc_identities$color,cex=0.8,text.font=2,bty='n')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

ylims = c(0,1.25)
plot(NA, NA, xlim=xlims, ylim=ylims, log='x', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=paste0('1:',xbigs), lwd=0, line=-0.5, las=2)
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=2, at=0:5/4, labels=NA, tck=-0.05, las=2)
axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, las=2, line=-0.5)
abline(h=c(1,.5,.1), lty=3, lwd=.75)
mtext(side=1, line=3, text='dilution', cex=0.8)
mtext(side=2, line=2.5, text='[PrP] (% Hi QC)', cex=0.8)
points(x=bh_parallelism$dilution, y=bh_parallelism$rel, col=bh_parallelism$color, pch=20)
for (smpl in unique(bh_plism_smry$detail)) {
  subs = bh_plism_smry[bh_plism_smry$detail==smpl,]
  points(x=subs$dilution, y=subs$mean_rel, col=subs$color, type='l', lwd=.75)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

### WT to KO dilution series
wtko = sqldf("
select   detail, ngml, ngml_trunc
from     elisa_raw
where    plate = 13
and      dilution = 200
and      detail not like '%QC%'
;")
wtko$nominal = as.numeric(gsub('%.*','',wtko$detail))/100
wtko$nominal_plot = wtko$nominal + rep(c(-0.0125,0.0125),7) # jitter the points just a bit so you can see both
multiplier = mean(wtko$ngml_trunc[wtko$nominal==1])
wtko$rel = wtko$ngml_trunc / multiplier
xlims = c(-0.05,1.05)
ylims = c(0,1.05)
par(mar=c(4,4,3,5))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=unique(wtko$nominal), labels=NA)
axis(side=1, at=sort(unique(wtko$nominal)), labels=c('100% KO','10/90','25/75','50/50','75/25','90/10','100% WT'),lwd=0,line=-0.5,las=2,cex.axis=0.9)
axis(side=2, at=unique(wtko$nominal), labels=NA)
axis(side=2, at=unique(wtko$nominal), labels=percent(unique(wtko$nominal)),lwd=0,line=-0.5,las=2,cex.axis=0.9)
axis(side=4, at=unique(wtko$nominal), labels=formatC(multiplier*unique(wtko$nominal),format='f',digits=0), las=2)
abline(v=unique(wtko$nominal),lwd=0.125)
abline(h=unique(wtko$nominal),lwd=0.125)
llq = 0.05 * 200 / multiplier
abline(h=llq, col='red', lwd=1.5, lty=2)
text(x=.5, y=llq, pos=3, labels='LLQ', col='red', lty=3)
mtext(side=1, line=2.5, text='WT/KO mix', cex=0.8)
mtext(side=2, line=3.5, text='PrP concentration (%WT)', cex=0.8)
par(xpd=T)
text(x=1.5, y=.5, srt=270, labels='PrP concentration (ng/g)', cex=1)
par(xpd=F)
points(wtko$nominal_plot, wtko$rel, pch=19, col=alpha(prpcol,ci_alpha))
abline(lm(rel ~ nominal, data=wtko), col=prpcol)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

csf_parallelism = elisa_raw[elisa_raw$plate==28 & elisa_raw$stype=='sample' & !grepl('QC',elisa_raw$detail),]
csf_plism_meta = data.frame(sample=c('9-IPC','CSF 1922','CSF 1847'), disp=c('interplate control','high CSF','low CSF'), color=c('#7fc97f','#beaed4','#fdc086'))
csf_parallelism$color = csf_plism_meta$color[match(csf_parallelism$detail, csf_plism_meta$sample)]
csf_plism_smry = sqldf("select detail, dilution, color, avg(ngml) mean from csf_parallelism group by 1,2,3 order by 1;")
csf_parallelism$reference = csf_plism_smry$mean[csf_plism_smry$detail=='9-IPC'][match(csf_parallelism$dilution, csf_plism_smry$dilution[csf_plism_smry$detail=='9-IPC'])]
csf_parallelism$rel = csf_parallelism$ngml / csf_parallelism$reference
csf_plism_smry = sqldf("select detail, dilution, color, avg(ngml) mean, avg(rel) mean_rel from csf_parallelism group by 1,2,3 order by 1;")

xlims = c(1, 100)
xats = rep(1:9,3) * rep(c(1,10,100), each=9)
xbigs = c(1,10,100)
ylims = c(0, 50)
par(mar=c(4,4,3,2))
plot(NA, NA, xlim=xlims, ylim=ylims, log='x', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=paste0('1:',xbigs), lwd=0, line=-0.5, las=2)
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=2, at=0:5*10, labels=NA, tck=-0.05, las=2)
axis(side=2, at=0:5*10, lwd=0, las=2, line=-0.5)
axis(side=2, at=0:10*5, labels=NA, tck=-0.025)
mtext(side=1, line=3, text='dilution', cex=0.8)
mtext(side=2, line=2.5, text='adjusted [PrP] (ng/mL)', cex=0.8)
points(x=csf_parallelism$dilution, y=csf_parallelism$ngml, col=csf_parallelism$color, pch=20)
for (smpl in unique(csf_plism_smry$detail)) {
  subs = csf_plism_smry[csf_plism_smry$detail==smpl,]
  points(x=subs$dilution, y=subs$mean, col=subs$color, type='l', lwd=.75)
}
llq = 0.05
ulq = 5
adjusted_llqs = llq * xats 
adjusted_ulqs = ulq * xats
points(x=xats, y=adjusted_llqs, col='red', type='l', lty=3, lwd=1.25)
points(x=xats, y=adjusted_ulqs, col='red', type='l', lty=3, lwd=1.25)
mtext(side=3, at=10, text='ULQ', col='red', cex=0.7)
mtext(side=4, las=2, at=5, text='LLQ', col='red', cex=0.7)
par(xpd=T)
legend('left',csf_plism_meta$disp,lwd=.75,pch=20,col=csf_plism_meta$color,text.col=csf_plism_meta$color,cex=0.8,text.font=2)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

### species cross-reactivity
species = sqldf("
                select   *
                from     elisa_raw
                where    stype = 'standard'
                and      (  (plate == 14 and standard_curve in ('MoPrP16','RaPrP50','CyPrP51'))
                or (plate == 16 and standard_curve in ('MoPrP16','HuPrP5'))
                or (plate == 18 and standard_curve in ('MoPrP16','BvPrP37','SHaPrP71')) )
                ;")

species$nominal = as.numeric(species$detail)
# make the adjustment for bank vole having been run based on the 1.15 ng/mL NanoDrop conc, later AAA came back at 1.11 ng/mL
species$nominal[species$standard_curve=='BvPrP37'] = species$nominal[species$standard_curve=='BvPrP37'] * 1.11/1.15
# do throw out the lowest std curve point, but don't throw out the zero
species = species[species$nominal != 0.02048,]


# re-calculate back conc using the Mo standard curve to compare
for (plateno in c(14,16,18)) {
  
  this_plate_rows = species$plate==plateno
  
  mocurve = species[this_plate_rows & species$standard_curve=='MoPrP16',]
  mocurve
  
  # fit standard curve - use standards from this standard curve, except the zero and any readings below baseline:
  fit = nlsLM(a450_bck ~ d + (a - d) / (1 + (nominal/c)^b), 
              start=list(a=a_init,b=b_init,c=c_init,d=d_init), 
              data=mocurve)
  
  species$ngml_by_mo[this_plate_rows] = calculate_conc(species$a450_bck[this_plate_rows], fit)
  
}

species_smry = sqldf("
                     select   standard_curve,
                     nominal,
                     avg(ngml_by_mo) mean,
                     avg(ngml_by_mo) - 1.96*stdev(ngml_by_mo)/count(*) l95,
                     avg(ngml_by_mo) + 1.96*stdev(ngml_by_mo)/count(*) u95
                     from     species
                     where    nominal > 0
                     group by 1, 2
                     order by 1, 2
                     ;")


par(mar=c(4,5,3,1))
limits = c(0.02, 7)
plot(NA, NA, xlim=limits, ylim=limits, log='xy', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=limits, labels=NA, lwd.ticks=0)
axis(side=2, at=limits, labels=NA, lwd.ticks=0)
axis(side=1, at=unique(species_smry$nominal), labels=NA, tck=-0.025)
axis(side=1, at=unique(species_smry$nominal), labels=formatC(unique(species_smry$nominal),format='f',digits=2), las=2, lwd=0, line=-0.5, cex.axis=0.8)
axis(side=2, at=unique(species_smry$nominal), labels=NA, tck=-0.025)
axis(side=2, at=unique(species_smry$nominal), labels=formatC(unique(species_smry$nominal),format='f',digits=2), las=2, lwd=0, line=-0.5, cex.axis=0.8)
mtext(side=1, line=2.0, text='nominal (ng/mL)', cex=0.8)
mtext(side=2, line=2.0, text='concentration back-calculated\nfrom mouse standard (ng/mL)', cex=0.8)
abline(v=unique(species_smry$nominal),lwd=0.125)
abline(h=unique(species_smry$nominal),lwd=0.125)

for (i in 1:nrow(species_meta)) {
  subs = subset(species_smry, standard_curve == species_meta$std[i])
  points(subs$nominal, subs$mean, pch=20, cex=0.75, col=species_meta$color[i])
  points(subs$nominal, subs$mean, type='l', lwd=0.75, col=species_meta$color[i])
}
legend('topleft',species_meta$name,pch=20,lwd=2,col=species_meta$color,text.col=species_meta$color, text.font=2, cex=0.8, bty='n')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


cs_vs_betaprion = read.table('data/assaydev/cs_vs_betaprion.tsv',sep='\t',header=T,quote='',comment.char='')
par(mar=c(4,4,3,1))
plot(NA, NA, xlim=c(0,150), ylim=c(0,700), xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=0:3*50)
axis(side=2, at=0:7*100, las=2)
mtext(side=1, line=2.5, text='cross-species PrP ELISA', cex=0.8)
mtext(side=2, line=3.0, text='BetaPrion ELISA', cex=0.8)
points(cs_vs_betaprion$cs_prp, cs_vs_betaprion$betaprion_prp, pch=20, col=alpha(prpcol, ci_alpha))
m = lm(betaprion_prp ~ cs_prp, data=cs_vs_betaprion)
abline(m, col=prpcol)
spearman_obj = suppressWarnings(cor.test(cs_vs_betaprion$betaprion_prp, cs_vs_betaprion$cs_prp, method='spearman'))
msg = paste0('rho = ',formatC(spearman_obj$estimate,format='f',digits=2),', P',format_p(spearman_obj$p.value))
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write(paste('Figure S2F results: N=',nrow(cs_vs_betaprion),' CSF samples analyzed by both our ELISA and BetaPrion, correlation was ',msg, ' and mean ratio was ',formatC(mean(cs_vs_betaprion$betaprion_prp / cs_vs_betaprion$cs_prp), digits=2),'\n',sep=''),text_stats_path,append=T)

silence_message = dev.off()






cat(file=stderr(), 'Done.\nCreating Figure S4E-F...')

### FIGURE S4E-F

pdf('display_items/vector/figure-s4-e-f.pdf',width=6.5, height=3.25)
xlims = c(0, 5.5)
ylims = c(0, 3)
par(mfrow=c(1,2), mar=c(4,4,3,1))

sha_reps = sqldf("
select   *
                 from     elisa_raw
                 where    stype = 'standard'
                 and      plate in (18,45) 
                 and      standard_curve in ('MoPrP16','SHaPrP71')
                 ;")
sha_reps$color = species_meta$color[match(sha_reps$standard_curve, species_meta$std)]
sha_reps$color[sha_reps$standard_curve=='MoPrP16'] = '#000000'
sha_reps$nominal = as.numeric(sha_reps$detail)

for (plate in c(18, 45)) {
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=unique(sha_reps$nominal), labels=formatC(unique(sha_reps$nominal), format='fg', digits=2))
  mtext(side=1, line=2.5, text='nominal')
  axis(side=2, at=0:3, las=2)
  mtext(side=2, line=2.5, text='absorbance')
  for (std in c('MoPrP16', 'SHaPrP71')) {
    subs = sha_reps[sha_reps$plate==plate & sha_reps$standard_curve==std,]
    points(x=subs$nominal, y=subs$a450_bck, pch=20, col=subs$color)
    fit = nlsLM(a450_bck ~ d + (a - d) / (1 + (nominal/c)^b),
                start=list(a=a_init,b=b_init,c=c_init,d=d_init),
                data=subs)
    y = 0:300/100
    x = calculate_conc(abs=y, fit)
    points(x, y, type='l', lwd=0.75, col=subs$color)
  }
  legend('topleft',
         c('mouse','hamster'),
         col=c('#000000',species_meta$color[species_meta$name=='hamster']),
         text.col=c('#000000',species_meta$color[species_meta$name=='hamster']),
         pch=19,lty=1,bty='n')
}
silence_message = dev.off()






cat(file=stderr(), 'Done.\nCreating Table S2...')

### TABLE S2

stab = elisa_raw[elisa_raw$plate==17,]
stab$stabsample = gsub(' - .*','',stab$detail)
stab$condition = gsub('.+ - ','',stab$detail)
stab$forceorder = 1
stab$forceorder[stab$condition == 'Freshly Thawed'] = 0

table_s2 = sqldf("
select   stabsample, condition, avg(ngml) mean, stdev(ngml)/avg(ngml) cv, count(*) n
from     stab
where    stype = 'sample' and detail not like '%QC%'
group by 1, 2
order by stabsample, forceorder
;")
table_s2$reference = table_s2$mean[table_s2$condition=='Freshly Thawed'][match(table_s2$stabsample, table_s2$stabsample[table_s2$condition=='Freshly Thawed'])]
table_s2$re = abs(table_s2$reference - table_s2$mean)/table_s2$reference
table_s2$re = percent(table_s2$re)
table_s2$mean = formatC(table_s2$mean, format='f', digits=1)
table_s2$cv = percent(table_s2$cv)

write.table(table_s2[,c('stabsample','condition','n','mean','cv','re')], 'display_items/table-s2.tsv', sep='\t', col.names=T, row.names=F, quote=F)







cat(file=stderr(), 'Done.\nCreating Figure 1...')

### FIGURE 1

brain_meta = read.table(sep='|',comment.char='',col.names=c('color','region'),textConnection("
#820BBB|hippocampus
#9715AF|parietal cortex BA7
#AB1EA4|prefrontal cortex
#C02898|frontal cortex
#D4318C|visual cortex
#458B00|thalamus
#3F741A|striatum
#395D33|putamen
#EE7600|cerebellum
#777777|white matter
#584E56|pons
#3B2E3A|medulla
#352734|olivary nucleus
"))

p67 = elisa[elisa$plate==67,]

hubrn = p67[!grepl('QC',p67$sample),]
hubrn$indiv = substr(hubrn$sample,1,4)
hubrn$region = trimws(substr(hubrn$sample,6,10))
hubrn_meta = data.frame(region=c('HP',"BA7", "CB", "ON", "Pu", "Thal", "WM"),
                        x=c(0,1,4,6,2,3,5),
                        disp=c('hippocampus','parietal cortex BA7','cerebellum','olivary nucleus','putamen','thalamus','white matter'))
hubrn_meta$color = brain_meta$color[match(hubrn_meta$disp, brain_meta$region)]
hubrn$x = hubrn_meta$x[match(hubrn$region, hubrn_meta$region)]
hubrn$disp = hubrn_meta$disp[match(hubrn$region, hubrn_meta$region)]
hubrn$color = hubrn_meta$color[match(hubrn$region, hubrn_meta$region)]
hubrn$y = max(hubrn$x) + 1 - hubrn$x

# hubrn_hiqc = p67$ngml_av[p67$sample=='Mo Pos Hi QC']
# hubrn$ngml_adj = hubrn$ngml_av / hubrn_hiqc

hubrn_smry = summarize_ci(hubrn, bycols=c('region','disp','y','color'), valcol='ngml_av')

human_region_aov = aov(ngml_av ~ region, data=hubrn)

write(paste('Figure 1 results: human brain regions Type I ANOVA P = ',formatC(summary(human_region_aov)[[1]]['region','Pr(>F)'],format='g',digits=2),'\n',sep=''),text_stats_path,append=T)

mobrn = elisa[elisa$plate %in% c(65,66) & grepl(' - ',elisa$sample) & !grepl('old',elisa$sample),]
mobrn$indiv = substr(mobrn$sample,1,2)
mobrn$region = substr(mobrn$sample,6,9)

mobrn_meta = data.frame(region=c('CB','HP','MC','Thal','PFC','Str','Pons','Medu'),
                        y=c(3,8,6,4,7,5,2,1),# c(1,6,4,2,5,3),
                        disp=c('cerebellum','hippocampus','visual cortex','thalamus','prefrontal cortex','striatum','pons','medulla'))
mobrn_meta$color = brain_meta$color[match(mobrn_meta$disp, brain_meta$region)]
mobrn$y = mobrn_meta$y[match(mobrn$region, mobrn_meta$region)]
mobrn$disp = mobrn_meta$disp[match(mobrn$region, mobrn_meta$region)]
mobrn$color = mobrn_meta$color[match(mobrn$region, mobrn_meta$region)]

# mobrn_hiqc = elisa[elisa$plate %in% c(65,66) & elisa$sample=='Mo Pos Hi QC',c('plate','ngml_av')]
# mobrn$ngml_adj = mobrn$ngml_av/mobrn_hiqc$ngml_av[match(mobrn$plate,mobrn_hiqc$plate)]

mobrn_smry = summarize_ci(mobrn, bycols=c('region','disp','y','color'), valcol='ngml_av')

mobrn$sex = substr(mobrn$indiv, 1,1)

mouse_region_aov = aov(ngml_av ~ region, data=mobrn)
write(paste('Figure 1 results: mouse brain regions Type I ANOVA P = ',formatC(summary(mouse_region_aov)[[1]]['region','Pr(>F)'],format='g',digits=2),'\n',sep=''),text_stats_path,append=T)


cybrn = elisa[elisa$plate==54 & grepl('-',elisa$sample),]
cybrn$indiv = substr(cybrn$sample,1,4)
cybrn$sex = substr(cybrn$sample,5,5)
cybrn$region = substr(cybrn$sample,7,11)



cybrn_meta = data.frame(region=c('CB','HP','PFC','Thal','Pons'),
                        y=c(2,5,4,3,1),
                        disp=c('cerebellum','hippocampus','frontal cortex','thalamus','pons'))
cybrn_meta$color = brain_meta$color[match(cybrn_meta$disp, brain_meta$region)]
cybrn$y = cybrn_meta$y[match(cybrn$region, cybrn_meta$region)]
cybrn$disp = cybrn_meta$disp[match(cybrn$region, cybrn_meta$region)]
cybrn$color = cybrn_meta$color[match(cybrn$region, cybrn_meta$region)]

# cybrn_hiqc = elisa$ngml_av[elisa$plate==54 & elisa$sample=='Mo Pos Hi QC']
# cybrn$ngml_adj = cybrn$ngml_av/cybrn_hiqc

cybrn_smry = summarize_ci(cybrn, bycols=c('region','disp','y','color'), valcol='ngml_av')

cyno_region_aov = aov(ngml_av ~ region, data=cybrn)
write(paste('Figure 1 results: cynomologus brain regions Type I ANOVA P = ',formatC(summary(cyno_region_aov)[[1]]['region','Pr(>F)'],format='g',digits=2),'\n',sep=''),text_stats_path,append=T)


# diagram
# see instructions: 
# http://help.brain-map.org/display/api/Downloading+and+Displaying+SVG
# this should work but it won't render in Chrome: http://api.brain-map.org/api/v2/svg_download/100960033?groups=28

resx=300
png('display_items/figure-1.png',width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(1,2,
                         3,4,
                         5,6), nrow=3, byrow=T)
layout(layout_matrix, heights=c(1, 1))

panel = 1


llq = 200*0.05


img = image_read('data/diagrams/human-regions-cropped.png')
img_png = image_convert(img, 'png')
par(mar=c(0,0,0.5,0))
plot(as.raster(img_png))
mtext(LETTERS[panel], side=3, cex=2, adj = 0.2, line = -1.75)
panel = panel + 1

par(mar=c(4,5,3,1))
ylims = range(hubrn_smry$y) + c(-0.5, 0.5)
xlims = c(0, 375)

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, at=hubrn_smry$y, text=hubrn_smry$disp, col=hubrn_smry$color, las=2, cex=0.8)
axis(side=1, at=0:4*100)
axis(side=1, at=0:40*10, lwd=0, lwd.ticks=1, tck=-0.025, labels=NA)
mtext(side=1, line=2.5, text='PrP (ng/g)')
abline(v=llq, col='red', lty=3, lwd=1.25)
mtext(side=1, line=0.25, at=llq, text='LLQ', col='red', cex=.6)
points(hubrn$ngml_av, hubrn$y, pch=19, col=hubrn$color)
for (each_indiv in unique(hubrn$indiv)) {
  subs = subset(hubrn, indiv==each_indiv)
  subs = subs[with(subs, order(y)),]
  points(subs$ngml_av, subs$y, type='l', lwd=0.125)  
}
barwidth=0.3
segments(y0=hubrn_smry$y-barwidth, y1=hubrn_smry$y+barwidth, x0=hubrn_smry$mean, col=hubrn_smry$color)
arrows(y0=hubrn_smry$y, x0=hubrn_smry$l95, x1=hubrn_smry$u95, col=hubrn_smry$color, code=3, angle=90, length=0.05)
mtext(side=3, line=0, text='human')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


img = image_read('data/diagrams/cynomolgus-regions-cropped.png')
img_png = image_convert(img, 'png')
par(mar=c(0,0,3,0))
plot(as.raster(img_png))
mtext(LETTERS[panel], side=3, cex=2, adj = 0.2, line = 0.5)
panel = panel + 1

par(mar=c(4,5,3,1))
ylims = range(cybrn_smry$y) + c(-0.5, 0.5)
xlims = c(0, 150)

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, at=cybrn_smry$y, text=cybrn_smry$disp, col=cybrn_smry$color, las=2, cex=0.8)
axis(side=1, at=0:3*100)
axis(side=1, at=0:30*10, lwd=0, lwd.ticks=1, tck=-0.025, labels=NA)
abline(v=llq, col='red', lty=3, lwd=1.25)
mtext(side=1, line=0.25, at=llq, text='LLQ', col='red', cex=.6)
mtext(side=1, line=2.5, text='PrP (ng/g)')
points(cybrn$ngml_av, cybrn$y, pch=19, col=cybrn$color)
for (each_indiv in unique(cybrn$indiv)) {
  subs = subset(cybrn, indiv==each_indiv)
  subs = subs[with(subs, order(y)),]
  points(subs$ngml_av, subs$y, type='l', lwd=0.125)  
}
barwidth=0.3
segments(y0=cybrn_smry$y-barwidth, y1=cybrn_smry$y+barwidth, x0=cybrn_smry$mean, col=cybrn_smry$color)
arrows(y0=cybrn_smry$y, x0=cybrn_smry$l95, x1=cybrn_smry$u95, col=cybrn_smry$color, code=3, angle=90, length=0.05)
mtext(side=3, line=0, text='cynomolgus')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



img = image_read('data/diagrams/mouse-regions-cropped.png')
img_png = image_convert(img, 'png')
par(mar=c(0,0,3,0))
plot(as.raster(img_png))
mtext(LETTERS[panel], side=3, cex=2, adj = 0.2, line = 0.5)
panel = panel + 1


par(mar=c(4,5,3,1))
ylims = range(mobrn_smry$y) + c(-0.5, 0.5)
xlims = c(0, 200)

plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, at=mobrn_smry$y, text=mobrn_smry$disp, col=mobrn_smry$color, las=2, cex=0.8)
axis(side=1, at=0:3*100)
axis(side=1, at=0:30*10, lwd=0, lwd.ticks=1, tck=-0.025, labels=NA)
abline(v=llq, col='red', lty=3, lwd=1.25)
mtext(side=1, line=0.25, at=llq, text='LLQ', col='red', cex=.6)
mtext(side=1, line=2.5, text='PrP (ng/g)')
points(mobrn$ngml_av, mobrn$y, pch=19, col=mobrn$color)
for (each_indiv in unique(mobrn$indiv)) {
  subs = subset(mobrn, indiv==each_indiv)
  subs = subs[with(subs, order(y)),]
  points(subs$ngml_av, subs$y, type='l', lwd=0.125)  
}
barwidth=0.3
segments(y0=mobrn_smry$y-barwidth, y1=mobrn_smry$y+barwidth, x0=mobrn_smry$mean, col=mobrn_smry$color)
arrows(y0=mobrn_smry$y, x0=mobrn_smry$l95, x1=mobrn_smry$u95, col=mobrn_smry$color, code=3, angle=90, length=0.05)
mtext(side=3, line=0, text='mouse')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


silence_message = dev.off()




cat(file=stderr(), 'Done.\nCreating Figure 2...')

tpm = read.table('data/gtex/prnp_tpm_t.txt',header=F,skip=2)
colnames(tpm) = c('sampid','tpm')

samp = read.table('data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t',header=T,quote='',comment.char='')
colnames(samp) = gsub('[^a-z0-9_]','_',tolower(colnames(samp)))

samp$subjid = substr(samp$sampid,1,10)

indiv = read.table('data/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt',sep='\t',header=T,quote='',comment.char='')
colnames(indiv) = gsub('[^a-z0-9_]','_',tolower(colnames(indiv)))

indiv$medage = as.integer(substr(indiv$age,1,1))*10+5
indiv$dthhrdy = as.factor(indiv$dthhrdy)

dat = sqldf("
            select   i.subjid, i.medage, i.dthhrdy, s.smts, s.smtsd, t.tpm
            from     tpm t, samp s, indiv i
            where    t.sampid = s.sampid
            and      s.subjid = i.subjid
            order by 1, 2
            ;")

all_dat = sqldf("
                select   i.subjid, i.medage, i.dthhrdy, s.smtsd, i.sex, count(*) n, avg(t.tpm) mean_tpm
                from     tpm t, samp s, indiv i
                where    t.sampid = s.sampid
                and      s.subjid = i.subjid
                group by 1, 2, 3, 4, 5
                order by 1, 2, 3, 4, 5
                ;")

# definition of hardy scale: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/variable.cgi?study_id=phs000424.v4.p1&phv=169092
all_dat$dthhrdy = as.factor(all_dat$dthhrdy)
all_dat$sex = ifelse(all_dat$sex==1,'M','F')
all_dat$medage = as.integer(all_dat$medage)
# p values for confounding:
sex_pval = 1-pchisq(q=as.numeric(chisq.test(table(all_dat[,c('dthhrdy','sex')]))$statistic),df=4)
age_pval = 1-pchisq(q=as.numeric(chisq.test(table(all_dat[,c('dthhrdy','medage')]))$statistic),df=20)

write(paste('Figure 2 results: in GTEx, sex is confounded with Hardy scale, Chi-square test P',format_p(sex_pval),'\n',sep=''),text_stats_path,append=T)
write(paste('Figure 2 results: in GTEx, age is confounded with Hardy scale, Chi-square test P',format_p(age_pval),'\n',sep=''),text_stats_path,append=T)

cors = sqldf("
             select   smtsd, count(*) n
             from     all_dat
             group by 1
             ;")

cors$spearman_r = as.numeric(NA)
cors$spearman_p = as.numeric(NA)

cors$lm_medage_annual_change = as.numeric(NA)
cors$lm_medage_annual_change_l95 = as.numeric(NA)
cors$lm_medage_annual_change_u95 = as.numeric(NA)
cors$lm_medage_p = as.numeric(NA)

cors$lm_sex_beta = as.numeric(NA)
cors$lm_sex_diff_l95 = as.numeric(NA)
cors$lm_sex_diff_u95 = as.numeric(NA)
cors$lm_sex_p = as.numeric(NA)

for (i in 1:nrow(cors)) {
  subs = subset(all_dat, smtsd==cors$smtsd[i])
  spearman_obj = suppressWarnings(cor.test(subs$medage, subs$mean_tpm, method='spearman', alternative='two.sided'))
  cors$spearman_r[i] = spearman_obj$estimate
  cors$spearman_p[i] = spearman_obj$p.value
  
  subs$log_tpm = log(subs$mean_tpm)
  
  if (length(unique(subs$sex))==2) {
    m = lm(log_tpm ~ medage + dthhrdy + sex, data=subs)
  } else {
    m = lm(log_tpm ~ medage + dthhrdy, data=subs)
  }
  
  cors$lm_medage_annual_change[i] = exp(summary(m)$coefficients['medage','Estimate'])-1
  cors$lm_medage_annual_change_l95[i] = exp(summary(m)$coefficients['medage','Estimate'] - 1.96*summary(m)$coefficients['medage','Std. Error'])-1
  cors$lm_medage_annual_change_u95[i] = exp(summary(m)$coefficients['medage','Estimate'] + 1.96*summary(m)$coefficients['medage','Std. Error'])-1
  cors$lm_medage_p[i] = summary(m)$coefficients['medage','Pr(>|t|)']
  
  if (length(unique(subs$sex))==2) {
    cors$lm_sex_diff[i] = exp(summary(m)$coefficients['sexM','Estimate'])-1
    cors$lm_sex_diff_l95[i] = exp(summary(m)$coefficients['sexM','Estimate'] - 1.96*summary(m)$coefficients['sexM','Std. Error'])-1
    cors$lm_sex_diff_u95[i] = exp(summary(m)$coefficients['sexM','Estimate'] + 1.96*summary(m)$coefficients['sexM','Std. Error'])-1
    cors$lm_sex_p[i] = summary(m)$coefficients['sexM','Pr(>|t|)']
  } else {
    cors$lm_sex_beta[i] = NA
    cors$lm_sex_diff[i] = NA
    cors$lm_sex_diff_l95[i] = NA
    cors$lm_sex_diff_u95[i] = NA
    cors$lm_sex_p[i] = NA
  }
}

cors$p_highlight = ''
cors$p_highlight[cors$lm_medage_p < 0.05/nrow(cors)] = '*'
cors$p_highlight[cors$lm_medage_p < 0.01/nrow(cors)] = '**'
cors$p_highlight[cors$lm_medage_p < 0.001/nrow(cors)] = '***'

cors$sex_p_highlight = ''
cors$sex_p_highlight[cors$lm_sex_p < 0.05/sum(!is.na(cors$lm_sex_p))] = '*'
cors$sex_p_highlight[cors$lm_sex_p < 0.01/sum(!is.na(cors$lm_sex_p))] = '**'
cors$sex_p_highlight[cors$lm_sex_p < 0.001/sum(!is.na(cors$lm_sex_p))] = '***'

write.table(cors, 'data/gtex/prnp_tpm_age_correlations.tsv',sep='\t',row.names=F,col.names=T,quote=F)

tissue_meta = read.table('data/gtex/tissue_metatissue.tsv',sep='\t',header=T,quote='',comment.char='')

cors$color = tissue_meta$color[match(cors$smtsd, tissue_meta$dispname)]
cors$y = nrow(cors):1



resx=300
png('display_items/figure-2.png',width=6.5*resx,height=8*resx,res=resx)

layout_matrix = matrix(c(rep(1,4), rep(2,4), rep(3,4),
                         rep(4,6), rep(6,6),
                         rep(5,6), rep(7,6)), nrow=3, byrow=T)
layout(layout_matrix, widths=c(1,1,1), heights=c(3,1,1))
#layout_matrix = matrix(1:3, nrow=1, byrow=T)
#layout(layout_matrix, widths=c(1.25,1,1))

panel = 1

p_offset = .25

ylims = range(cors$y) + c(-0.5, 0.5)
xlims = c(0,1)
par(mar=c(4,1,3,0))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
mtext(side=2, at=cors$y, text=cors$smtsd, line=-15, las=2, cex=.6)

par(mar=c(4,0,3,3))
#xlims = c(floor(min(cors$lm_medage_annual_change)*100-1)/100, ceiling(max(cors$lm_medage_annual_change)*100+1)/100)
#xats = seq(min(xlims), max(xlims), .01)
xlims = c(-.025, .025)
xats = seq(-.02,.02,.01)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, labels=signed_percent(xats,digits=0))
abline(v=xats, lwd=.125)
abline(v=0, lwd=1)
mtext(side=1, line=2.25, text='annual change', cex=0.8)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
#mtext(side=2, at=cors$y, text=cors$smtsd, line=0.25, las=2)
points(x=cors$lm_medage_annual_change, y=cors$y, pch=19, col=cors$color)
segments(x0=cors$lm_medage_annual_change_l95, x1=cors$lm_medage_annual_change_u95, y0=cors$y, lwd=3, col=cors$color)
mtext(side=4, at=cors$y-p_offset, text=cors$p_highlight, las=2, line=0.25)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

par(mar=c(4,0,3,3))
xlims = c(-.4, .4) #c(floor(min(cors$lm_sex_diff_l95)*100-1)/100, ceiling(max(cors$lm_sex_diff_u95)*100+1)/100)
xats = seq(min(xlims), max(xlims), .1)
ylims = range(cors$y) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, labels=signed_percent(abs(xats),digits=0))
abline(v=xats, lwd=.125)
abline(v=0, lwd=1)
mtext(side=1, at=min(xlims)/2, line=2.25, text='female bias', cex=.8)
mtext(side=1, at=max(xlims)/2, line=2.25, text='male bias', cex=.8)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
#mtext(side=2, at=cors$y, text=cors$smtsd, line=0.25, las=2)
points(x=cors$lm_sex_diff, y=cors$y, pch=19, col=cors$color)
segments(x0=cors$lm_sex_diff_l95, x1=cors$lm_sex_diff_u95, y0=cors$y, lwd=3, col=cors$color)
mtext(side=4, at=cors$y-p_offset, text=cors$sex_p_highlight, las=2, line=0.25)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

human_csf_llq = 0.05 * 40

if (file.exists('ignore/')) {
  mgh_by_demo = read.table('ignore/mgh_by_demo.tsv', sep='\t', header=T)
} else {
  cat(file=stderr(), 'sensitive patient data unavailable, leaving panels 2C-2D empty...')
  mgh_by_demo = data.frame(sex=character(0), decade=integer(0), prp_mean=numeric(0))
}

sex_params = data.frame(sex=c('M','F'),disp=c('male','female'),x=c(1,2),color=c('#226eb2','#b13024')) # the colors are from Aguet 2019 GTEx v8 paper
mgh_by_demo$color = sex_params$color[match(mgh_by_demo$sex, sex_params$sex)]
mgh_by_demo$sex_x = sex_params$x[match(mgh_by_demo$sex, sex_params$sex)]
mgh_by_demo$sex_x_plot = jitter(mgh_by_demo$sex_x, amount=.125)
mgh_by_demo_smry = summarize_ci(mgh_by_demo, bycols=c('sex','sex_x','color'), valcol=c('prp_mean'))

barwidth=0.25  

par(mar=c(2,5,3,2))
xlims = c(0.5, 2.5)
ylims = c(0,120)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=0.25, at=sex_params$x, text=sex_params$disp)
abline(h=human_csf_llq, col='red', lty=3, lwd=1.25)
mtext(side=4, line=0.25, at=human_csf_llq, text='LLQ', las=2, col='red', cex=.6)
axis(side=2, at=0:20*10, tck=-0.025, labels=NA, las=2)
axis(side=2, at=0:2*50, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3.0, text='CSF [PrP] (ng/mL)')
points(mgh_by_demo$sex_x_plot, mgh_by_demo$prp_mean, pch=20, col=alpha(mgh_by_demo$color, ci_alpha))
segments(x0=mgh_by_demo_smry$sex_x-barwidth, x1=mgh_by_demo_smry$sex_x+barwidth, y0=mgh_by_demo_smry$mean, lwd=1, col=mgh_by_demo_smry$color)
arrows(x0=mgh_by_demo_smry$sex_x, y0=mgh_by_demo_smry$l95, y1=mgh_by_demo_smry$u95, code=3, angle=90, length=0.1, col=mgh_by_demo_smry$color)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

mgh_by_decade_smry = summarize_ci(mgh_by_demo, bycols=c('decade'), valcol=c('prp_mean'))

decade_params = data.frame(decade=c(20,30,40,50,60), disp=c('20-29','30-39','40-49','50-59','60+'))
set.seed(1)
mgh_by_demo$decade_plot = jitter(mgh_by_demo$decade, amount=1)

xlims = c(15,65)
ylims = c(0,120)
par(mar=c(4,5,3,2))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=0.25, at=decade_params$decade, text=decade_params$disp)
mtext(side=1, line=2.5, text='age')
abline(h=human_csf_llq, col='red', lty=3, lwd=1.25)
mtext(side=4, line=0.25, at=human_csf_llq, text='LLQ', las=2, col='red', cex=.6)
axis(side=2, at=0:20*10, tck=-0.025, labels=NA, las=2)
axis(side=2, at=0:2*50, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3.0, text='CSF [PrP] (ng/mL)')
points(mgh_by_demo$decade_plot, mgh_by_demo$prp_mean, pch=20, col=alpha('#000000', ci_alpha))
barwidth=2.5
segments(x0=mgh_by_decade_smry$decade-barwidth, x1=mgh_by_decade_smry$decade+barwidth, y0=mgh_by_decade_smry$mean, lwd=1, col='#000000')
arrows(x0=mgh_by_decade_smry$decade, y0=mgh_by_decade_smry$l95, y1=mgh_by_decade_smry$u95, code=3, angle=90, length=0.05, col='#000000')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

rat_brain_llq = 0.05 * 200
rat_csf_llq = 0.05 * 20

agerabrn = elisa[elisa$plate==34 & !grepl('QC',elisa$sample),]
agerat_meta = read.table('data/samples/PRP210114_meta.tsv',sep='\t',header=T)
agerabrn$prefix = substr(agerabrn$sample,1,1)
agerabrn$age = agerat_meta$age[match(agerabrn$prefix, agerat_meta$prefix)]
agerabrn_smry = summarize_ci(agerabrn, bycols='age', valcol='ngml_av')     

m = lm(ngml_av ~ age, data=agerabrn)
write(paste('Figure 2 results: aged rat brain, linear regression of [PrP] vs age P',format_p(summary(m)$coefficients['age','Pr(>|t|)']),'\n'),text_stats_path,append=T)

ageracsf = elisa[elisa$plate==35 & !grepl('QC',elisa$sample),]
agerat_meta = read.table('data/samples/PRP210114_meta.tsv',sep='\t',header=T)
ageracsf$prefix = substr(ageracsf$sample,1,1)
ageracsf$age = agerat_meta$age[match(ageracsf$prefix, agerat_meta$prefix)]
ageracsf_smry = summarize_ci(ageracsf, bycols='age', valcol='ngml_av')     

m = lm(ngml_av ~ age, data=ageracsf)
write(paste('Figure 2 results: aged rat CSF, linear regression of [PrP] vs age P',format_p(summary(m)$coefficients['age','Pr(>|t|)']),'\n'),text_stats_path,append=T)

barwidth = 5

brncol = '#5A6351'
csfcol = '#8A2BE2'

par(mar=c(1,4,3,2))

plot(NA, NA, xlim=c(0,366), ylim=c(0,150), axes=F, ann=F, xaxs='i', yaxs='i')
#axis(side=1, at=seq(0,366,30.44), labels=seq(0,12,1))
#mtext(side=1, line=2.5, text='age (months)')
axis(side=1, at=c(0,366), labels=NA, lwd.ticks=0)
axis(side=2, at=seq(0,150,50), las=2)
mtext(side=2, line=2.5, text='brain [PrP] (ng/g)')
abline(h=rat_brain_llq, col='red', lty=3, lwd=1.25)
mtext(side=4, line=0.0, at=rat_brain_llq, text='LLQ', las=2, col='red', cex=.6)
points(agerabrn$age, agerabrn$ngml_av, pch=20, col=alpha(brncol,ci_alpha))
segments(x0=agerabrn_smry$age, y0=agerabrn_smry$l95, y1=agerabrn_smry$u95, lwd=1.5)
segments(x0=agerabrn_smry$age-barwidth, x1=agerabrn_smry$age+barwidth, y0=agerabrn_smry$mean, lwd=1.5)
m = lm(ngml_av ~ age, data=agerabrn)
abline(m, col=brncol, lwd=.75)
mtext(side=3, line=0, text=paste0('P = ',formatC(summary(m)$coefficients['age','Pr(>|t|)'],format='f',digits=2)))
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

par(mar=c(4,4,3,2))
plot(NA, NA, xlim=c(0,366), ylim=c(0,15), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,366,30.44), labels=seq(0,12,1), cex.axis=0.8)
mtext(side=1, line=2.5, text='age (months)')
axis(side=2, at=seq(0,15,5), las=2)
mtext(side=2, line=2.5, text='CSF [PrP] (ng/mL)')
abline(h=rat_csf_llq, col='red', lty=3, lwd=1.25)
mtext(side=4, line=0.0, at=rat_csf_llq, text='LLQ', las=2, col='red', cex=.6)
points(ageracsf$age, ageracsf$ngml_av, pch=20, col=alpha(csfcol,ci_alpha))
segments(x0=ageracsf_smry$age, y0=ageracsf_smry$l95, y1=ageracsf_smry$u95, lwd=1.5)
segments(x0=ageracsf_smry$age-barwidth, x1=ageracsf_smry$age+barwidth, y0=ageracsf_smry$mean, lwd=1.5)
m = lm(ngml_av ~ age, data=ageracsf)
abline(m, col=csfcol, lwd=.75)
mtext(side=3, line=0, text=paste0('P = ',formatC(summary(m)$coefficients['age','Pr(>|t|)'],format='f',digits=2)))
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

silence_message = dev.off()


















if (file.exists('ignore/')) {
  
  cat(file=stderr(), 'Done.\nCreating Figure 3...')
  
### FIGURE 3

resx=300
png('display_items/figure-3.png',width=6.5*resx,height=5.5*resx,res=resx)

layout_matrix = matrix(c(1,1,1,1,1,1,1,
                         2,3,4,5,6,7,8,
                         9,10,11,12,13,14,15), nrow=3, byrow=T)
layout(layout_matrix,  heights=c(1.25,1,1))

panel = 1

indiv_means = read.table('ignore/mgh_indiv_means.tsv',sep='\t',header=T,quote='',comment.char='')

write(paste('Figure 3 results: in N=',sum(indiv_means$n > 1),' subjects with >1 lumbar puncture, CSF [PrP] mean test-retest CV = ',percent(mean(indiv_means$cv[indiv_means$n > 1]),digits=3),'\n',sep=''),text_stats_path,append=T)

indiv_means %>%
  group_by(gtdisp, x, color) %>%
  summarize(rel_mean = mean(prp_rel), rel_sd=sd(prp_rel), n=n(),
            rel_l95=mean(prp_rel) - 1.96*sd(prp_rel)/sqrt(n()),
            rel_u95=mean(prp_rel) + 1.96*sd(prp_rel)/sqrt(n()),
            abs_mean = mean(prp_mean), abs_sd=sd(prp_mean), 
            abs_l95=mean(prp_mean) - 1.96*sd(prp_mean)/sqrt(n()),
            abs_u95=mean(prp_mean) + 1.96*sd(prp_mean)/sqrt(n()),
            .groups='keep') -> gt_means

set.seed(1)
indiv_means$x_plot = jitter(indiv_means$x, amount=.25)

xlims = range(indiv_means$x) + c(-0.5, 0.5)
ylims = c(0, 1.65)
par(mar=c(1,8,3,6))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=0.5, at=gt_means$x, text=gt_means$gtdisp, cex=1.0, col=gt_means$color, font=2)
axis(side=2, at=0:20/10, labels=NA, tck=-0.025)
axis(side=2, at=0:4/2, labels=percent(0:4/2), tck=-0.05, las=2)
par(xpd=T)
text(x=-0.25, y=mean(ylims), srt=90, labels='relative CSF [PrP]\n(% control average)', cex=1.1)
par(xpd=F)
barwidth=.4
#rect(xleft=gt_means$x-barwidth, xright=gt_means$x+barwidth, ybottom=rep(0,nrow(gt_means)), ytop=gt_means$rel_mean, col=alpha(gt_means$color, ci_alpha), border=NA)
segments(x0=gt_means$x-barwidth, x1=gt_means$x+barwidth, y0=gt_means$rel_mean, col=gt_means$color, lwd=2)
arrows(x0=gt_means$x, y0=gt_means$rel_l95, y1=gt_means$rel_u95, col=gt_means$color, lwd=2, code=3, angle=90, length=0.1)
points(x=indiv_means$x_plot, y=indiv_means$prp_rel, col=alpha(indiv_means$color, ci_alpha), pch=19)
abline(h=c(1,.5), lty=3, lwd=.5)
multiplier = gt_means$abs_mean[gt_means$gtdisp=='no mutation']
abs_ylims = c(0, ceiling(max(indiv_means$prp_mean)/10)*10)
axis(side=4, at=seq(0,max(abs_ylims),10)/multiplier, tck=-0.025, labels=NA)
axis(side=4, at=seq(0,max(abs_ylims),50)/multiplier, tck=-0.05, labels=seq(0,max(abs_ylims),50), las=2)
#mtext(side=4, line=3, text='absolute CSF [PrP]\n(ng/mL)')
abline(h=human_csf_llq/multiplier, col='red', lty=3, lwd=1.25)
text(x=0.65, y=human_csf_llq/multiplier, label='LLQ', pos=3, col='red', cex=.9)
#mtext(side=4, line=0.25, at=human_csf_llq/multiplier, text='LLQ', las=2, col='red', cex=.6)
par(xpd=T)
text(x=6.0, y=mean(ylims), srt=270, labels='absolute CSF [PrP]\n(ng/mL)', cex=1.1)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

### MGH study samples MRM data
peptides = read.table('data/mrm/ref/peptides.tsv',sep='\t',header=T,quote='',comment.char='')
expts = read.table('data/mrm/processed/expt28.tsv',sep='\t',header=T,quote='',comment.char='')[0,]
metas = read.table('ignore/expt28.tsv',sep='\t',header=T,quote='',comment.char='')[0,]
for (exptno in 28:32) {
  expt = read.table(paste0('data/mrm/processed/expt',exptno,'.tsv'),sep='\t',header=T,quote='',comment.char='')
  expt$run = exptno
  meta = read.table(paste0('ignore/expt',exptno,'.tsv'),sep='\t',header=T,quote='',comment.char='')
  meta$run = exptno
  expts = rbind(expts, expt)
  metas = rbind(metas, meta)
}
expts$ln = expts$light_area/expts$n15_area
expts$lh = expts$light_area/expts$heavy_area
expts$indiv_visit = metas$sample[match(expts$id, metas$label)]
expts$indiv = substr(expts$indiv_visit,1,4)

cv1 = sqldf("
            select   run, peptide, indiv, indiv_visit, avg(ln) mean_ln, stdev(ln)/avg(ln) cv_ln, avg(lh) mean_lh, stdev(lh)/avg(lh) cv_lh
            from     expts
            group by 1, 2, 3, 4
            ;")

pep_tech_perf = sqldf("
                      select   peptide,
                      avg(mean_ln) mean_ln1, avg(cv_ln) mean_techrep_cv_ln,
                      avg(mean_lh) mean_lh1, avg(cv_lh) mean_techrep_cv_lh
                      from     cv1
                      group by 1
                      order by 1
                      ;")

test_retest = sqldf("
                    select   peptide, indiv, 
                    avg(mean_ln) mean_ln1, stdev(mean_ln)/avg(mean_ln) testretest_cv_ln,
                    avg(mean_lh) mean_lh1, stdev(mean_lh)/avg(mean_lh) testretest_cv_lh
                    from     cv1
                    where    run = 28 -- the run where we did test-retest
                    group by 1, 2
                    ;")

pep_retest = sqldf("
                   select   p.ncorder, t.peptide, 
                   avg(mean_ln1) mean_ln, avg(testretest_cv_ln) mean_testretest_cv_ln,
                   avg(mean_lh1) mean_lh, avg(testretest_cv_lh) mean_testretest_cv_lh
                   from     test_retest t, peptides p
                   where    t.peptide = p.peptide
                   and      t.peptide in (select peptide from peptides where human = 'yes')
                   group by 1, 2
                   order by 1
                   ;")

pep_perf_and_retest = sqldf("
                            select   rt.ncorder, rt.peptide,
                            rt.mean_ln, p.mean_techrep_cv_ln, rt.mean_testretest_cv_ln,
                            rt.mean_lh, p.mean_techrep_cv_lh, rt.mean_testretest_cv_lh
                            from     pep_retest rt, pep_tech_perf p
                            where    rt.peptide = p.peptide
                            and      rt.peptide in (select peptide from peptides where human ='yes')
                            order by 1
                            ;")

write.table(format(pep_perf_and_retest[,c('ncorder','peptide','mean_ln','mean_techrep_cv_ln','mean_testretest_cv_ln')], digits=2),
            'display_items/table_s4.tsv', sep='\t', col.names=T, quote=F, row.names=F)

pep_tech_perf$x = peptides$ncorder[match(pep_tech_perf$peptide, peptides$peptide)]
pep_tech_perf = pep_tech_perf[pep_tech_perf$peptide %in% peptides$peptide[peptides$human=='yes'],]

all_mgh_mrm = expts

mgh_samples = read.table('ignore/gt.tsv',sep='\t',header=T)
mgh_samples$gt = mgh_samples$mut
mgh_samples$gt[mgh_samples$mut=='other/unknown'] = 'other'
mgh_samples$gt[mgh_samples$mut=='none'] = 'no mutation'

all_mgh_mrm$visit = substr(all_mgh_mrm$indiv_visit,6,6)
all_mgh_mrm$gt = mgh_samples$gt[match(all_mgh_mrm$indiv, mgh_samples$id)]
all_mgh_mrm$gt_exact = mgh_samples$empirical_mutation[match(all_mgh_mrm$indiv, mgh_samples$id)]

all_mgh_mrm$codon = suppressWarnings(as.integer(gsub('[A-Z]*','',all_mgh_mrm$gt_exact)))
all_mgh_mrm$mismatch = all_mgh_mrm$codon %in% 195:204 & all_mgh_mrm$peptide =='GENFTETDVK' |
  all_mgh_mrm$codon %in% 209:220 & all_mgh_mrm$peptide == 'VVEQMCITQYER'

all_mgh_mrm$flag[!is.na(all_mgh_mrm$flag) & all_mgh_mrm$flag==''] = NA

all_mgh_mrm = all_mgh_mrm[all_mgh_mrm$peptide %in% peptides$peptide[peptides$human=='yes'],]

mrm_indiv_means = sqldf("
                    select   peptide, indiv, gt, gt_exact, mismatch, avg(ln) mean_ln, stdev(ln)/avg(ln) cv_ln,
                    avg(lh) mean_lh, stdev(lh)/avg(lh) cv_lh, count(*) n
                    from     all_mgh_mrm
                    where    gt is not null 
                    and      flag is null
                    and      peptide in (select peptide from peptides where human='yes')
                    group by 1, 2, 3, 4, 5
                    ;")

bygt = sqldf("
             select   peptide, gt, count(*) n,
             avg(mean_ln) mean_ln1, stdev(mean_ln)/avg(mean_ln) cv_ln, 
             avg(mean_lh) mean_lh1, stdev(mean_lh)/avg(mean_lh) cv_lh
             from     mrm_indiv_means
             group by 1, 2
             order by 1, 2
             ;")

indiv_rel = sqldf("
                  select   a.peptide, a.indiv, a.gt, a.gt_exact, a.mismatch,
                  avg(a.ln) mean_ln, stdev(a.ln)/avg(a.ln) cv_ln,
                  avg(a.lh) mean_lh, stdev(a.lh)/avg(a.lh) cv_lh, count(*) n
                  from     all_mgh_mrm a
                  where    (a.gt is not null)
                  and      (a.peptide in (select peptide from peptides where human='yes'))
                  and      a.flag is null
                  group by 1, 2, 3, 4, 5
                  having   stdev(a.ln)/avg(a.ln) < .15 -- filter on L/15N CV but not L/H CV
                  ;")

mrm_nomut_mean = sqldf("
                   select   peptide, avg(mean_ln) ln_nomut, avg(mean_lh) lh_nomut, count(*) n
                   from     indiv_rel
                   where    gt = 'no mutation'
                   group by 1
                   order by 1
                   ;")

indiv_rel$ln_nomut = mrm_nomut_mean$ln_nomut[match(indiv_rel$peptide, mrm_nomut_mean$peptide)]
indiv_rel$lh_nomut = mrm_nomut_mean$lh_nomut[match(indiv_rel$peptide, mrm_nomut_mean$peptide)]

indiv_rel$ln_rel = indiv_rel$mean_ln / indiv_rel$ln_nomut
indiv_rel$lh_rel = indiv_rel$mean_lh / indiv_rel$lh_nomut

rel_smry = sqldf("
                 select     peptide, gt, count(*) n,
                 avg(ln_rel) mean_ln_rel, stdev(ln_rel) sd_ln_rel,
                 avg(lh_rel) mean_lh_rel, stdev(lh_rel) sd_lh_rel
                 from     indiv_rel
                 group by 1, 2
                 order by 1, 2
                 ;")

rel_smry$ln_rel_l95 = rel_smry$mean_ln_rel - 1.96 * rel_smry$sd_ln_rel/sqrt(rel_smry$n)
rel_smry$ln_rel_u95 = rel_smry$mean_ln_rel + 1.96 * rel_smry$sd_ln_rel/sqrt(rel_smry$n)

rel_smry$lh_rel_l95 = rel_smry$mean_lh_rel - 1.96 * rel_smry$sd_lh_rel/sqrt(rel_smry$n)
rel_smry$lh_rel_u95 = rel_smry$mean_lh_rel + 1.96 * rel_smry$sd_lh_rel/sqrt(rel_smry$n)

params = read.table('data/samples/mgh_study_params.tsv',sep='\t',header=T,quote='',comment.char='')
params = params[1:5,]

indiv_rel$color = params$color[match(indiv_rel$gt, params$disp)]
indiv_rel$x = params$x[match(indiv_rel$gt, params$disp)]
rel_smry$color = params$color[match(rel_smry$gt, params$disp)]
rel_smry$x = params$x[match(rel_smry$gt, params$disp)]

xlims = range(params$x) + c(-0.5, 0.5)
ylims = c(0, 1.9)

par(mar=c(0.1,0.5,2,0.5))
plot(NA, NA, xlim=c(0,1), ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, line=-7, at=0:4/2, labels=percent(0:4/2), las=2, lwd=0, lwd.ticks=1, tck=-0.1)
axis(side=2, line=-7, at=0:20/10, labels=NA,lwd=0, lwd.ticks=1, tck=-0.025)
mtext(side=2, line=-4, text=expression(atop('relative L/'^'15'*'N',
                                       '(% control average)')), cex=0.8)
mtext(LETTERS[panel], side=3, cex=2, adj = 0.2, line = 0.5)
panel = panel + 1

barwidth = .33
for (this_peptide in peptides$peptide[peptides$human=='yes']) {
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
  #axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  #axis(side=2, at=0:3/2, labels=percent(0:3/2,digits=0), las=2, lwd=0, lwd.ticks=1)
  abline(h=1, lwd=.75, lty=3)
  #mtext(side=1, line=0.25, text=this_peptide, cex=0.6)
  smry_subs = rel_smry[rel_smry$peptide==this_peptide,]
  # rect(xleft=smry_subs$x - barwidth, xright=smry_subs$x + barwidth, ybottom=rep(0, nrow(smry_subs)), ytop=smry_subs$mean_ln_rel, col=alpha(smry_subs$color, ci_alpha), border=NA)
  # arrows(x0=smry_subs$x, y0=smry_subs$ln_rel_l95, y1=smry_subs$ln_rel_u95, angle=90, length=0.025, code=3, col='black')
  segments(x0=smry_subs$x - barwidth, x1=smry_subs$x + barwidth, y0=smry_subs$mean_ln_rel, col=smry_subs$color)
  arrows(x0=smry_subs$x, y0=smry_subs$ln_rel_l95, y1=smry_subs$ln_rel_u95, angle=90, length=0.025, code=3, col=smry_subs$color)
  indiv_subs = indiv_rel[indiv_rel$peptide==this_peptide,]
  points(x=indiv_subs$x, y=indiv_subs$ln_rel, pch=19, col=alpha(indiv_subs$color, ci_alpha))
  points(x=indiv_subs$x[indiv_subs$mismatch], y=indiv_subs$ln_rel[indiv_subs$mismatch], pch=0, cex=1.25, col='red')
  # mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  # panel = panel + 1
  
}

indiv_rel$elisa_rel = indiv_means$prp_rel[match(indiv_rel$indiv, indiv_means$indiv)]

xlims = ylims = c(0, 1.9)

par(mar=c(3,0.5,1.5,0.5))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, line=-7, at=0:4/2, labels=percent(0:4/2), las=2, lwd=0, lwd.ticks=1, tck=-0.1)
axis(side=2, line=-7, at=0:20/10, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
# line breaks within an expression: https://stackoverflow.com/a/20549830/3806692
mtext(side=2, line=-4, text=expression(atop('relative L/'^'15'*'N',
                                            '(% control average)')), cex=0.8)
mtext(LETTERS[panel], side=3, cex=2, adj = 0.2, line = 0.5)
panel = panel + 1

barwidth = .33
for (this_peptide in peptides$peptide[peptides$human=='yes']) {
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=2, at=ylims, lwd=1,lwd.ticks=0, labels=NA)
  axis(side=1, at=xlims, lwd=1,lwd.ticks=0, labels=NA)
  axis(side=1, at=0:4/2, labels=percent(0:4/2), lwd=0, lwd.ticks=0, line=-0.75)
  axis(side=1, at=0:4/2, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=1, at=0:20/10, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
  mtext(side=1, line=1.5, text='ELISA PrP', cex=0.6)
  abline(h=1, lwd=.75, lty=3)
  abline(v=1, lwd=.75, lty=3)
  abline(a=0, b=1, lwd=.25, col='#777777')
  if (this_peptide=='GENFTETDVK') {
    abline(a=0, b=0.5, lwd=.25, col='#FF7777')
  }
  mtext(side=3, line=.5, text=this_peptide, cex=0.6)
  mtext(side=3, line=-.25, text=peptides$codons_hu[peptides$peptide==this_peptide], cex=0.6)
  indiv_subs = indiv_rel[indiv_rel$peptide==this_peptide,]
  points(x=indiv_subs$elisa_rel, y=indiv_subs$ln_rel, pch=19, col=alpha(indiv_subs$color, ci_alpha))
  points(x=indiv_subs$elisa_rel[indiv_subs$mismatch], y=indiv_subs$ln_rel[indiv_subs$mismatch], pch=0, cex=1.25, col='red')
  # mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  # panel = panel + 1
  
}

silence_message = dev.off()

p102l_ks = ks.test(indiv_means$prp_mean[indiv_means$gtdisp=='no mutation'], indiv_means$prp_mean[indiv_means$gtdisp=='P102L'])
d178n_ks = ks.test(indiv_means$prp_mean[indiv_means$gtdisp=='no mutation'], indiv_means$prp_mean[indiv_means$gtdisp=='D178N'])
e200k_ks = ks.test(indiv_means$prp_mean[indiv_means$gtdisp=='no mutation'], indiv_means$prp_mean[indiv_means$gtdisp=='E200K'])

write(paste('Figure 3 results: P102L vs. no mutation CSF [PrP] KS test P',format_p(p102l_ks$p.value),'\n',sep=''),text_stats_path,append=T)
write(paste('Figure 3 results: D178N vs. no mutation CSF [PrP] KS test P',format_p(d178n_ks$p.value),'\n',sep=''),text_stats_path,append=T)
write(paste('Figure 3 results: E200K vs. no mutation CSF [PrP] KS test P',format_p(e200k_ks$p.value),'\n',sep=''),text_stats_path,append=T)

cat(file=stderr(), 'Done.\nCreating Figure S5...')

eqtl_c129 = read.table('data/gtex/chr20_4699605_A_G_multitissue_eqtl_data.tsv',sep='\t',header=T)
eqtl_ue = read.table('data/gtex/chr20_4613886_C_T_multitissue_eqtl_data.tsv',sep='\t',header=T)
eqtl = sqldf("select c.tissue, c.nes c_nes, c.l95 c_l95, c.u95 c_u95, u.nes u_nes, u.l95 u_l95, u.u95 u_u95 from eqtl_c129 c, eqtl_ue u where c.tissue = u.tissue order by c.tissue;")
eqtl$tmatch = gsub('[^a-z]','',tolower(eqtl$tissue))
cors$tmatch = gsub('[^a-z]','',tolower(cors$smtsd))
eqtl$disp = cors$smtsd[match(eqtl$tmatch, cors$tmatch)]
eqtl$y = cors$y[match(eqtl$tmatch, cors$tmatch)]
eqtl$color = cors$color[match(eqtl$tmatch, cors$tmatch)]

eqtl_c129_stats = read.table('data/gtex/gtex_multitissue_rs1799990.csv',sep=',',header=T)
colnames(eqtl_c129_stats) = gsub('[^a-z0-9_]','_',tolower(colnames(eqtl_c129_stats)))
eqtl$c_p = eqtl_c129_stats$p_value[match(eqtl$disp, eqtl_c129_stats$tissue)]
eqtl$c_p_symb = ''
eqtl$c_p_symb[eqtl$c_p < 1e-5] = '*'
eqtl$c_p_symb[eqtl$c_p < 1e-6] = '**'
eqtl$c_p_symb[eqtl$c_p < 1e-7] = '***'

eqtl_ue_stats = read.table('data/gtex/gtex_multitissue_rs17327121.csv',sep=',',header=T)
colnames(eqtl_ue_stats) = gsub('[^a-z0-9_]','_',tolower(colnames(eqtl_ue_stats)))
eqtl$u_p = eqtl_ue_stats$p_value[match(eqtl$disp, eqtl_ue_stats$tissue)]
eqtl$u_p_symb = ''
eqtl$u_p_symb[eqtl$u_p < 1e-5] = '*'
eqtl$u_p_symb[eqtl$u_p < 1e-6] = '**'
eqtl$u_p_symb[eqtl$u_p < 1e-7] = '***'

resx=300
png('display_items/figure-s5.png',width=6.5*resx,height=8.5*resx,res=resx)

layout_matrix = matrix(c(rep(1,4), rep(2,4), rep(3,4),
                         rep(4,4), rep(5,4), rep(6,4)
                         ), nrow=2, byrow=T)
layout(layout_matrix, heights=c(3,1))

panel = 1

p_offset = .25

ylims = range(eqtl$y) + c(-0.5, 0.5)
xlims = c(0,1)
par(mar=c(4,1,3,0))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
mtext(side=2, at=cors$y, text=cors$smtsd, line=-15, las=2, cex=.6)

xlims = c(-.7, .7)
xats = seq(-.7, .7, .1)
xbigs = seq(-.5, .5, .5)

par(mar=c(4,0,3,3))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=ifelse(xbigs > 0,paste0('+',formatC(xbigs,format='f',digits=1)),formatC(xbigs,format='f',digits=1)), lwd.ticks=0, line=-0.25)
abline(v=xats, lwd=.125)
abline(v=0, lwd=1)
mtext(side=1, line=3.25, at=0, text='NES', cex=0.8)
mtext(side=1, line=2.25, at=.5, text='V higher', cex=0.7)
mtext(side=1, line=2.25, at=-.5, text='M higher', cex=0.7)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
points(x=eqtl$c_nes, y=eqtl$y, pch=19, col=eqtl$color)
segments(x0=eqtl$c_l95, x1=eqtl$c_u95, y0=eqtl$y, lwd=3, col=eqtl$color)
mtext(side=4, at=eqtl$y-p_offset, text=eqtl$c_p_symb, las=2, line=0.25)
panel = 1
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)

par(mar=c(4,0,3,3))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, labels=ifelse(xbigs > 0,paste0('+',formatC(xbigs,format='f',digits=1)),formatC(xbigs,format='f',digits=1)), lwd.ticks=0, line=-0.25)
abline(v=xats, lwd=.125)
abline(v=0, lwd=1)
mtext(side=1, line=3.25, at=0, text='NES', cex=0.8)
mtext(side=1, line=2.25, at=.5, text='T higher', cex=0.7)
mtext(side=1, line=2.25, at=-.5, text='C higher', cex=0.7)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
points(x=eqtl$u_nes, y=eqtl$y, pch=19, col=eqtl$color)
segments(x0=eqtl$u_l95, x1=eqtl$u_u95, y0=eqtl$y, lwd=3, col=eqtl$color)
mtext(side=4, at=eqtl$y-p_offset, text=eqtl$u_p_symb, las=2, line=0.25)
panel = 3
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)

plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F)

par(mar = c(3,4,3,1))

indiv_means %>%
  group_by(codon129, cx, ccolor) %>%
  summarize(abs_mean = mean(prp_mean), abs_sd=sd(prp_mean), n=n(),
            abs_l95=mean(prp_mean) - 1.96*sd(prp_mean)/sqrt(n()),
            abs_u95=mean(prp_mean) + 1.96*sd(prp_mean)/sqrt(n()),
            .groups='keep') -> c129_means
c129_means$abs_mean[c129_means$n==1] = NA
set.seed(1)
indiv_means$cx_plot = jitter(indiv_means$cx, amount=.25)

xlims = range(c129_means$cx) + c(-0.5, 0.5)
ylims = c(0,ceiling(max(indiv_means$prp_mean)/10)*10)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=1, at=c129_means$cx, text=c129_means$codon129, col=c129_means$ccolor, font=2)
axis(side=2, at=0:20*10, labels=NA, tck=-0.025)
axis(side=2, at=0:4*50, tck=-0.05, las=2)
mtext(side=2, line=2.5, text='CSF [PrP] (ng/mL)')
barwidth = 0.4
segments(x0=c129_means$cx-barwidth, x1=c129_means$cx+barwidth, y0=c129_means$abs_mean, col=c129_means$ccolor, lwd=2)
arrows(x0=c129_means$cx, y0=c129_means$abs_l95, y1=c129_means$abs_u95, angle=90, length=0.05, code=3, col=c129_means$ccolor, lwd=2)
set.seed(1)
points(x=indiv_means$cx_plot, y=indiv_means$prp_mean, pch=20, cex=.8, col=alpha(indiv_means$ccolor, ci_alpha))
panel = 2
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)

c129_ks = ks.test(indiv_means$prp_mean[indiv_means$codon129=='MM'], indiv_means$prp_mean[indiv_means$codon129=='MV'])
c129_ks_nomut = ks.test(indiv_means$prp_mean[indiv_means$codon129=='MM' & indiv_means$gtdisp=='no mutation'], indiv_means$prp_mean[indiv_means$codon129=='MV' & indiv_means$gtdisp=='no mutation'])

write(paste('Figure S5 results: codon 129 MM vs. MV CSF [PrP] KS test P',format_p(c129_ks$p.value),' in all individuals or P',format_p(c129_ks_nomut$p.value),' restricted to individuals without mutations\n',sep=''),text_stats_path,append=T)


indiv_means %>%
  group_by(ue, ux, ucolor) %>%
  summarize(abs_mean = mean(prp_mean), abs_sd=sd(prp_mean), n=n(),
            abs_l95=mean(prp_mean) - 1.96*sd(prp_mean)/sqrt(n()),
            abs_u95=mean(prp_mean) + 1.96*sd(prp_mean)/sqrt(n()),
            .groups='keep') -> ue_means
ue_means$abs_mean[ue_means$n==1] = NA
set.seed(1)
indiv_means$ux_plot = jitter(indiv_means$ux, amount=.25)


xlims = range(ue_means$ux, na.rm=T) + c(-0.5, 0.5)
ylims = c(0,ceiling(max(indiv_means$prp_mean)/10)*10)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=1, at=ue_means$ux, text=ue_means$ue, col=ue_means$ucolor, font=2)
axis(side=2, at=0:20*10, labels=NA, tck=-0.025)
axis(side=2, at=0:4*50, tck=-0.05, las=2)
mtext(side=2, line=2.5, text='CSF [PrP] (ng/mL)')
barwidth = 0.4
segments(x0=ue_means$ux-barwidth, x1=ue_means$ux+barwidth, y0=ue_means$abs_mean, col=ue_means$ucolor, lwd=2)
arrows(x0=ue_means$ux, y0=ue_means$abs_l95, y1=ue_means$abs_u95, angle=90, length=0.05, code=3, col=ue_means$ucolor, lwd=2)
set.seed(1)
points(x=indiv_means$ux_plot[!is.na(indiv_means$ue)], y=indiv_means$prp_mean[!is.na(indiv_means$ue)], pch=20, cex=.8, col=alpha(indiv_means$ucolor[!is.na(indiv_means$ue)], ci_alpha))
panel = 4
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)

eqtl_cc_ct_ks = ks.test(indiv_means$prp_mean[indiv_means$ue=='CC'], indiv_means$prp_mean[indiv_means$ue=='CT'])
eqtl_ct_tt_ks = ks.test(indiv_means$prp_mean[indiv_means$ue=='CT'], indiv_means$prp_mean[indiv_means$ue=='TT'])
eqtl_cc_tt_ks = ks.test(indiv_means$prp_mean[indiv_means$ue=='CC'], indiv_means$prp_mean[indiv_means$ue=='TT'])

write(paste('Figure S5 results: rs17327121 eQTL CC vs. CT CSF [PrP] KS test P',format_p(eqtl_cc_ct_ks$p.value),'\n',sep=''),text_stats_path,append=T)
write(paste('Figure S5 results: rs17327121 eQTL CT vs. TT CSF [PrP] KS test P',format_p(eqtl_ct_tt_ks$p.value),'\n',sep=''),text_stats_path,append=T)
write(paste('Figure S5 results: rs17327121 eQTL CC vs. TT CSF [PrP] KS test P',format_p(eqtl_cc_tt_ks$p.value),'\n',sep=''),text_stats_path,append=T)

silence_message = dev.off()

cat(file=stderr(), 'Done.\n')

} else {
  cat(file=stderr(), 'Done.\nSensitive patient data unavailable, skipping Figures 3 and S5.\n')
}






# are there proxies for ue?
gwas = read.table('data/gwas/cjd_gwas_rs17327121_50kb.tsv',sep='\t',header=T)
ld = read.table('data/gwas/ldlink_proxy20576.txt',header=T)
colnames(ld) = gsub('[^a-z0-9_]','_',tolower(colnames(ld)))
inner_join(gwas, ld, by=c('rs_id' = 'rs_number')) %>%
  arrange(-r2) ->
  gwas_ld

write(paste('Figure S5 results: highest LD with rs17327121 for any SNP in the CJD GWAS data: r^2=',max(gwas_ld$r2),'\n',sep=''),text_stats_path,append=T)


cat(file=stderr(), 'Creating Figure 4...')

### FIGURE 4

modr_2wk = read.table('data/samples/PRP210122.tsv',sep='\t',header=T)
modr_2wk$rna = modr_2wk$mrna/100
modr_2wk$animal = gsub('-','',modr_2wk$animal)
#elisa[elisa$plate==38 & grepl('QC',elisa$sample),]
cs38 = elisa[elisa$plate==38 & !grepl('QC|-',elisa$sample),]
modr_2wk$prp = cs38$ngml_av[match(modr_2wk$animal, cs38$sample)]
saline_mean_prp = mean(modr_2wk$prp[modr_2wk$dose==0])
modr_2wk$prp_rel = modr_2wk$prp / saline_mean_prp
modr_2wk_params = read.table('data/samples/PRP210122_meta.tsv',sep='\t',header=T,quote='',comment.char='')
modr_2wk$color = modr_2wk_params$color[match(modr_2wk$dose, modr_2wk_params$dose)]

smry_mrna = summarize_ci(modr_2wk, valcol='rna',bycols='dose')
smry_prot = summarize_ci(modr_2wk, valcol='prp_rel',bycols='dose')

modr_2wk_smry = sqldf("
select   r.dose, r.mean rna_mean, r.l95 rna_l95, r.u95 rna_u95, p.mean prot_mean, p.l95 prot_l95, p.u95 prot_u95, pm.color 
from     smry_mrna r, smry_prot p, modr_2wk_params pm
where    r.dose = p.dose and r.dose = pm.dose
order by 1;")

### RATS
samples = read.table('data/mrm/samples/rats-4wk-dr.tsv',sep='\t',header=T)
samparams = read.table('data/mrm/samples/dr_cohorts.tsv',sep='\t',header=T,quote='',comment.char='')

cs22 = elisa[elisa$plate==22,]
cs22_raw = elisa_raw[elisa_raw$plate==22,]
cs22$dose = samples$dose[match(cs22$sample, samples$animal)]
cs22$rna = samples$ipsi_hemisphere_prnp_mrna[match(cs22$sample, samples$animal)]/100
cs22$color = samparams$color[match(cs22$dose, samparams$dose)]

cs24 = elisa[elisa$plate==24,]
cs24_raw = elisa_raw[elisa_raw$plate==24,]
samples = rbind(samples, data.frame(animal='D4/5',
                                    tx=samples$tx[samples$animal=='D4'],
                                    dose=samples$dose[samples$animal=='D4'],
                                    ipsi_hemisphere_prnp_mrna=mean(samples$ipsi_hemisphere_prnp_mrna[samples$animal %in% c('D4','D5')]))
)
cs24$dose = samples$dose[match(cs24$sample, samples$animal)]
cs24$rna = samples$ipsi_hemisphere_prnp_mrna[match(cs24$sample, samples$animal)]/100
cs24$color = samparams$color[match(cs24$dose, samparams$dose)]


brain_csf = sqldf("
                  select   b.sample, b.dose, b.color, b.ngml_av b_ngml, c.ngml_av c_ngml
                  from     cs22 b, cs24 c
                  where    (b.sample = c.sample or (b.sample in ('D4','D5') and c.sample = 'D4/5'))
                  and      b.dose is not null
                  order by 2, 1
                  ")

brain_csf$b_rel = brain_csf$b_ngml / mean(brain_csf$b_ngml[brain_csf$dose==0])
brain_csf$c_rel = brain_csf$c_ngml / mean(brain_csf$c_ngml[brain_csf$dose==0])


samparams = read.table('data/mrm/samples/dr_cohorts.tsv',sep='\t',header=T,quote='',comment.char='')
pepparams = read.table('data/mrm/ref/peptides.tsv',sep='\t',header=T,quote='',comment.char='')
samples = read.table('data/mrm/samples/rats-4wk-dr.tsv',sep='\t',header=T)
samples = rbind(samples, data.frame(animal='D4/5',tx=samples$tx[samples$animal=='D4'], dose=samples$dose[samples$animal=='D4'],
                                    ipsi_hemisphere_prnp_mrna=mean(samples$ipsi_hemisphere_prnp_mrna[samples$animal %in% c('D4','D5')])))
samples = samples[with(samples, order(animal)),]


expt26 = read.table('data/mrm/processed/expt26.tsv',sep='\t',header=T,quote='',comment.char='')
meta26 = read.table('data/mrm/meta/expt26.tsv',sep='\t',header=T,quote='',comment.char='')

expt27 = read.table('data/mrm/processed/expt27.tsv',sep='\t',header=T,quote='',comment.char='')
meta27 = read.table('data/mrm/meta/expt27.tsv',sep='\t',header=T,quote='',comment.char='')

meta26$dose = samples$dose[match(meta26$sample, samples$animal)]
meta26$txcolor = samparams$color[match(meta26$dose, samparams$dose)]
expt26$pepcolor = pepparams$color[match(expt26$peptide, pepparams$peptide)]
expt26$rat_expected = pepparams$rat[match(expt26$peptide, pepparams$peptide)]
expt26$lh = expt26$light_area / expt26$heavy_area
expt26$txcolor = meta26$txcolor[match(expt26$sample,meta26$skyline_id)]
expt26$dose = meta26$dose[match(expt26$sample,meta26$skyline_id)]
expt26$animal = meta26$sample[match(expt26$sample, meta26$skyline_id)]
expt26$pepx = pepparams$allorder[match(expt26$peptide, pepparams$peptide)]

meta27$dose = samples$dose[match(meta27$sample, samples$animal)]
meta27$txcolor = samparams$color[match(meta27$dose, samparams$dose)]
expt27$pepcolor = pepparams$color[match(expt27$peptide, pepparams$peptide)]
expt27$rat_expected = pepparams$rat[match(expt27$peptide, pepparams$peptide)]
expt27$lh = expt27$light_area / expt27$heavy_area
expt27$txcolor = meta27$txcolor[match(expt27$sample,meta27$skyline_id)]
expt27$dose = meta27$dose[match(expt27$sample,meta27$skyline_id)]
expt27$animal = meta27$sample[match(expt27$sample, meta27$skyline_id)]
expt27$pepx = pepparams$allorder[match(expt27$peptide, pepparams$peptide)]

norm26 = sqldf("
               select   e.peptide, e.pepx, e.dose, e.animal, e.txcolor, avg(e.lh) mean_lh, avg(e.lh / control.mean_saline_lh) norm_lh
               from     expt26 e, (select peptide, avg(lh) mean_saline_lh from expt26 where dose = 0 and (flag = '' or flag is null) group by 1 order by 1) control
               where    e.peptide = control.peptide
               and      (e.flag = '' or e.flag is null)
               and      e.animal in (select animal from samples)
               group by 1, 2, 3, 4, 5
               ;")
norm27 = sqldf("
               select   e.peptide, e.pepx, e.dose, e.animal, e.txcolor, avg(e.lh) mean_lh, avg(e.lh / control.mean_saline_lh) norm_lh
               from     expt27 e, (select peptide, avg(lh) mean_saline_lh from expt27 where dose = 0 and (flag = '' or flag is null) group by 1 order by 1) control
               where    e.peptide = control.peptide
               and      (e.flag = '' or e.flag is null)
               and      e.animal in (select animal from samples)
               group by 1, 2, 3, 4, 5
               ;")

# need to make this a full outer join because different samples in A (the reference 0 dose group) are missing for brain vs. csf
brain_csf_mrm = sqldf("
select   b.peptide, b.pepx, b.dose, b.animal, b.txcolor, b.norm_lh brain_norm_lh, c.norm_lh csf_norm_lh
from     norm27 b, norm26 c
where    b.peptide = c.peptide
and      (b.animal = c.animal or (b.animal in ('D4','D5') and c.animal = 'D4/5'))
;")

norm26 %>%
  group_by(peptide, pepx, dose, txcolor) %>%
  summarize(.groups='keep', n=n(), 
            csf_mean=mean(norm_lh, na.rm=T), 
            csf_sd=sd(norm_lh, na.rm=T), 
            csf_l95=mean(norm_lh, na.rm=T)-1.96*sd(norm_lh, na.rm=T)/sqrt(n()), 
            csf_u95=mean(norm_lh, na.rm=T)+1.96*sd(norm_lh, na.rm=T)/sqrt(n())) -> csf_mrm_smry

norm27 %>%
  group_by(peptide, pepx, dose, txcolor) %>%
  summarize(.groups='keep', n=n(), 
            brain_mean=mean(norm_lh, na.rm=T), 
            brain_sd=sd(norm_lh, na.rm=T), 
            brain_l95=mean(norm_lh, na.rm=T)-1.96*sd(norm_lh, na.rm=T)/sqrt(n()), 
            brain_u95=mean(norm_lh, na.rm=T)+1.96*sd(norm_lh, na.rm=T)/sqrt(n())) -> brain_mrm_smry

inner_join(csf_mrm_smry, brain_mrm_smry, 
           by=c('peptide', 'pepx', 'dose', 'txcolor')) -> brain_csf_mrm_smry

modr_4wk_params = read.table('data/samples/PRP210608_meta.tsv',sep='\t',header=T,quote='',comment.char='')
modr_4wk_rna = read.table('data/samples/PRP210608.tsv',sep='\t',header=T)
elisa %>%
  filter(plate==47) %>%
  inner_join(modr_4wk_rna, by=c('sample'='animal')) %>%
  inner_join(modr_4wk_params, by=c('dose'='dose')) %>%
  select(animal=sample, dose, rna, ngml_av, color) -> modr_4wk

saline_mean = mean(modr_4wk$ngml_av[modr_4wk$dose==0])
modr_4wk$rel = modr_4wk$ngml_av / saline_mean

modr_4wk %>%
  group_by(dose, color) %>%
  summarize(.groups='keep', 
            mean=mean(ngml_av), 
            l95=mean(ngml_av) - 1.96*sd(ngml_av)/sqrt(n()),
            u95=mean(ngml_av) + 1.96*sd(ngml_av)/sqrt(n()),
            mean_rel=mean(rel), 
            l95_rel=mean(rel) - 1.96*sd(rel)/sqrt(n()),
            u95_rel=mean(rel) + 1.96*sd(rel)/sqrt(n()),
            mean_rna=mean(rna), 
            l95_rna=mean(rna) - 1.96*sd(rna)/sqrt(n()),
            u95_rna=mean(rna) + 1.96*sd(rna)/sqrt(n())) -> modr_4wk_smry

modr_60dpi_meta = suppressMessages(read_tsv('data/samples/PRP210322.tsv'))
modr_60dpi_params = read.table('data/samples/PRP210322_meta.tsv',sep='\t',header=T,quote='',comment.char='')
elisa %>%
  filter(plate==46) %>%
  full_join(modr_60dpi_meta, by=c('sample'='animal')) %>%
  inner_join(modr_60dpi_params, by=c('dose'='dose')) %>%
  select(animal=sample, dose, ngml_av, color, rna) %>%
  arrange(dose, animal) -> modr_60dpi
saline_prot_mean = mean(modr_60dpi$ngml_av[modr_60dpi$dose==0])
modr_60dpi$rel = modr_60dpi$ngml_av / saline_prot_mean


modr_60dpi %>%
  group_by(dose, color) %>%
  summarize(.groups='keep', 
            mean=mean(ngml_av, na.rm=T), 
            l95=mean(ngml_av, na.rm=T) - 1.96*sd(ngml_av, na.rm=T)/sqrt(n()),
            u95=mean(ngml_av, na.rm=T) + 1.96*sd(ngml_av, na.rm=T)/sqrt(n()),
            mean_rel=mean(rel, na.rm=T), 
            l95_rel=mean(rel, na.rm=T) - 1.96*sd(rel, na.rm=T)/sqrt(n()),
            u95_rel=mean(rel, na.rm=T) + 1.96*sd(rel, na.rm=T)/sqrt(n()),
            mean_rna=mean(rna, na.rm=T),
            l95_rna=mean(rna, na.rm=T) - 1.96*sd(rna, na.rm=T)/sqrt(sum(!is.na(rna))),
            u95_rna=mean(rna, na.rm=T) + 1.96*sd(rna, na.rm=T)/sqrt(sum(!is.na(rna)))) -> modr_60dpi_smry

resx=300
png('display_items/figure-4.png',width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(rep(1,4), rep(2,4), rep(3,4),
                         rep(4,6), rep(5,6),
                         rep(6:11, each=2)), nrow=3, byrow=T)
layout(layout_matrix, heights=c(1.25, 1.25, 1))

panel = 1


par(mar=c(4,5,3,1))
panel = 1
xlims = c(0,1.25)
ylims = c(0,1.25)
startplot(xlims, ylims)
axis(side=1, at=0:3/2, labels=percent(0:3/2), cex.axis=1)
mtext(side=1, line=2.5, text=expression('brain'~italic('Prnp')~'mRNA')) # ipsi
axis(side=2, at=0:6/4, labels=percent(0:6/4), cex.axis=1, las=2)
mtext(side=2, line=3, text='brain PrP') # contra
points(modr_2wk$rna, modr_2wk$prp_rel, col=alpha(modr_2wk$color,ci_alpha), pch=19)
segments(x0=modr_2wk_smry$rna_l95, x1=modr_2wk_smry$rna_u95, y0=modr_2wk_smry$prot_mean,  lwd=2, col=modr_2wk_smry$color)
segments(x0=modr_2wk_smry$rna_mean, y0=modr_2wk_smry$prot_l95, y1=modr_2wk_smry$prot_u95, lwd=2, col=modr_2wk_smry$color)
abline(v=1, lwd=.75, lty=3)
abline(h=1, lwd=.75, lty=3)
abline(a=0, b=1, lwd=.25)
modr_2wk$p_kd = 1- modr_2wk$prp_rel
modr_2wk$r_kd = 1- modr_2wk$rna
m = lm(p_kd ~ r_kd + 0, data=modr_2wk)
abline(a=1-m$coefficients['r_kd'], b=m$coefficients['r_kd'],col='blue',lwd=.5)
legend('bottomright',legend=modr_2wk_smry$dose,col=modr_2wk_smry$color,pch=19,text.col=modr_2wk_smry$color,title='dose (µg)',title.col='#000000',cex=0.85)
mtext(side=3, line=0.25, text='naïve mice\n2 weeks')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
write(paste('Figure 4 results: mice 2 weeks post-dose each 1% brain RNA knockdown equals ',paste0(formatC(coefficients(m)['r_kd'],digits=2),'%'),' brain PrP knockdown\n',sep=''),text_stats_path,append=T)

par(mar=c(4,5,3,1))
xlims = c(0,1.25)
ylims = c(0,1.25)
startplot(xlims, ylims)
axis(side=1, at=0:3/2, labels=percent(0:3/2), cex.axis=1)
mtext(side=1, line=2.5, text=expression('brain'~italic('Prnp')~'mRNA')) # ipsi
axis(side=2, at=0:6/4, labels=percent(0:6/4), cex.axis=1, las=2)
mtext(side=2, line=3, text='brain PrP') # contra
points(modr_4wk$rna, modr_4wk$rel, col=alpha(modr_4wk$color,ci_alpha), pch=19)
segments(x0=modr_4wk_smry$l95_rna, x1=modr_4wk_smry$u95_rna, y0=modr_4wk_smry$mean_rel,  lwd=2, col=modr_4wk_smry$color)
segments(x0=modr_4wk_smry$mean_rna, y0=modr_4wk_smry$l95_rel, y1=modr_4wk_smry$u95_rel, lwd=2, col=modr_4wk_smry$color)
abline(v=1, lwd=.75, lty=3)
abline(h=1, lwd=.75, lty=3)
abline(a=0, b=1, lwd=.25)
modr_4wk$p_kd = 1- modr_4wk$rel
modr_4wk$r_kd = 1- modr_4wk$rna
m = lm(p_kd ~ r_kd + 0, data=modr_4wk)
abline(a=1-m$coefficients['r_kd'], b=m$coefficients['r_kd'],col='blue',lwd=.5)
legend('bottomright',legend=modr_4wk_smry$dose,col=modr_4wk_smry$color,pch=19,text.col=modr_4wk_smry$color,title='dose (µg)',title.col='#000000',cex=0.85)
mtext(side=3, line=0.25, text='naïve mice\n4 weeks')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
write(paste('Figure 4 results: mice 4 weeks post-dose each 1% brain RNA knockdown equals ',paste0(formatC(coefficients(m)['r_kd'],digits=2),'%'),' brain PrP knockdown\n',sep=''),text_stats_path,append=T)

xlims = c(0,1.65)
ylims = c(0,1.25)
par(mar=c(4,5,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:3/2, labels=percent(0:3/2), cex.axis=1)
mtext(side=1, line=2.5, text=expression('brain'~italic('Prnp')~'mRNA')) 
axis(side=2, at=0:6/4, labels=percent(0:6/4), cex.axis=1, las=2)
mtext(side=2, line=3, text='brain PrP') 
abline(v=1, lty=3)
abline(h=1, lty=3)
abline(a=0, b=1, lwd=.25)
points(modr_60dpi$rna, modr_60dpi$rel, pch=19, col=alpha(modr_60dpi$color, ci_alpha))
segments(x0=modr_60dpi_smry$l95_rna, x1=modr_60dpi_smry$u95_rna, y0=modr_60dpi_smry$mean_rel, col=modr_60dpi_smry$color, lwd=2)
segments(x0=modr_60dpi_smry$mean_rna, y0=modr_60dpi_smry$l95_rel, y1=modr_60dpi_smry$u95_rel, col=modr_60dpi_smry$color, lwd=2)
legend('bottomright',legend=modr_60dpi_smry$dose,col=modr_60dpi_smry$color,pch=19,text.col=modr_60dpi_smry$color,title='dose (µg)',title.col='#000000',cex=0.85)
mtext(side=3, line=0.25, text='infected mice\n4 weeks')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

cs22$b_rel = cs22$ngml_av / mean(cs22$ngml_av[cs22$dose==0 & !is.na(cs22$dose)])
cs22 = cs22[!is.na(cs22$dose),]
xlims = c(0,1.25)
ylims = c(0,1.25)
par(mar=c(4,5,3,1))
startplot(xlims, ylims)
axis(side=1, at=0:6/4, labels=percent(0:6/4), cex.axis=0.8)
mtext(side=1, line=2.5, text=expression('brain'~italic('Prnp')~'mRNA (% saline)')) # ipsi
axis(side=2, at=0:6/4, labels=percent(0:6/4), cex.axis=0.8, las=2)
mtext(side=2, line=3, text='brain PrP (% saline)') # contra
points(cs22$rna, cs22$b_rel, col=alpha(cs22$color, ci_alpha), pch=19)
cs22 %>% 
  filter(!is.na(dose)) %>%
  group_by(dose, color) %>% 
  summarize( .groups='keep', n=n(), 
            rna_mean=mean(rna), rna_sd=sd(rna), rna_l95=mean(rna)-1.96*sd(rna)/sqrt(n()), rna_u95=mean(rna)+1.96*sd(rna)/sqrt(n()),
            prot_mean=mean(b_rel), prot_sd=sd(b_rel), prot_l95=mean(b_rel)-1.96*sd(b_rel)/sqrt(n()), prot_u95=mean(b_rel)+1.96*sd(b_rel)/sqrt(n())) -> 
  cs22_smry
segments(x0=cs22_smry$rna_l95,  x1=cs22_smry$rna_u95,  y0=cs22_smry$prot_mean,  lwd=2, col=cs22_smry$color)
segments(x0=cs22_smry$rna_mean, y0=cs22_smry$prot_l95, y1=cs22_smry$prot_u95, lwd=2,   col=cs22_smry$color)
abline(v=1, lwd=.75, lty=3)
abline(h=1, lwd=.75, lty=3)
abline(a=0, b=1, lwd=.25)
legend('bottomright',legend=samparams$dose,col=samparams$color,pch=19,text.col=samparams$color,title='dose (µg)',title.col='#000000',cex=0.85)
mtext(side=3, line=0.25, text='naïve rats 4 weeks')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

xlims = c(0,1.25)
ylims = c(0,1.25)
par(mar=c(4,5,3,1))
startplot(xlims, ylims)
axis(side=1, at=0:5/4, labels=percent(0:5/4), cex.axis=0.8)
mtext(side=1, line=2.5, text='brain PrP (% saline)') # contra
axis(side=2, at=0:5/4, labels=percent(0:5/4), cex.axis=0.8, las=2)
mtext(side=2, line=3, text='CSF PrP (% saline)')
points(brain_csf$b_rel, brain_csf$c_rel, pch=19, col=alpha(brain_csf$color, ci_alpha))
brain_csf %>% 
  filter(!is.na(dose)) %>%
  group_by(dose, color) %>% 
  summarize( .groups='keep', n=n(), 
             csf_mean=mean(c_rel), csf_sd=sd(c_rel), csf_l95=mean(c_rel)-1.96*sd(c_rel)/sqrt(n()), csf_u95=mean(c_rel)+1.96*sd(c_rel)/sqrt(n()),
             brain_mean=mean(b_rel), brain_sd=sd(b_rel), brain_l95=mean(b_rel)-1.96*sd(b_rel)/sqrt(n()), brain_u95=mean(b_rel)+1.96*sd(b_rel)/sqrt(n())) -> 
  brain_csf_smry
segments(x0=brain_csf_smry$brain_mean, y0=brain_csf_smry$csf_l95,  y1=brain_csf_smry$csf_u95, lwd=2, col=brain_csf_smry$color)
segments(x0=brain_csf_smry$brain_l95, x1=brain_csf_smry$brain_u95, y0=brain_csf_smry$csf_mean, lwd=2,   col=brain_csf_smry$color)
abline(v=1, lwd=.75, lty=3)
abline(h=1, lwd=.75, lty=3)
abline(a=0, b=1, lwd=.25)
legend('bottomright',legend=samparams$dose,col=samparams$color,pch=19,text.col=samparams$color,title='dose (µg)',title.col='#000000',cex=0.85)
mtext(side=3, line=0.25, text='naïve rats 4 weeks')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

xlims = c(0,1.5)
ylims = c(0,1.5)

par(mar=c(4,1,3,0))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, line=-7.5, at=0:6/4, labels=percent(0:6/4), las=2, lwd=0, lwd.ticks=0, cex.axis=1)
axis(side=2, line=-7, at=0:6/4, labels=NA, las=2, lwd=0, lwd.ticks=1, tck=-0.025, cex.axis=1)
mtext(side=2, line=-5, text='CSF')
legend('bottomleft',legend=samparams$dose,col=samparams$color,pch=19,text.col=samparams$color,title='dose (µg)',title.col='#000000',cex=0.85,bty='n')

par(mar=c(4,0,3,1))
barwidth = 0.33
for (peptide in peptides$peptide[peptides$rat=='yes']) {
  plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
  abline(h=1, lty=3, lwd=0.5)
  abline(v=1, lty=3, lwd=0.5)
  axis(side=1, at=0:3/2, labels=NA)
  axis(side=1, at=c(0,.5,1), labels=c('0%','50%','100%'), cex.axis=0.9)
  mtext(side=1, line=2.5, text='brain')
  axis(side=2, at=c(0,1.5), labels=NA, lwd.ticks=0)
  abline(a=0, b=1, lwd=.25)
  #axis(side=2, at=0:6/4, labels=percent(0:6/4,digits=0), las=2)
  #mtext(side=2, line=3, text='CSF')
  rows = brain_csf_mrm$peptide==peptide
  points(x=brain_csf_mrm$brain_norm_lh[rows], y=brain_csf_mrm$csf_norm_lh[rows], col=alpha(brain_csf_mrm$txcolor[rows],ci_alpha),pch=19)
  smry_subs = brain_csf_mrm_smry[brain_csf_mrm_smry$peptide==peptide,]
  segments(x0=smry_subs$brain_mean, y0=smry_subs$csf_l95,   y1=smry_subs$csf_u95,  lwd=1.25,  col=smry_subs$txcolor)
  segments(x0=smry_subs$brain_l95,  x1=smry_subs$brain_u95, y0=smry_subs$csf_mean, lwd=1.25,  col=smry_subs$txcolor)
  
  mtext(side=3, line=0, text=peptide, cex=0.6)
  mtext(side=3, line=-.75, text=paste0(peptides$codons_ra[peptides$peptide==peptide],'   '), cex=0.6)
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.2, line = 0.75)
  panel = panel + 1
  #  smry_rows = smry27$peptide==peptide
  #  segments(x0=smry27$midx[smry_rows]-barwidth, x1=smry27$midx[smry_rows]+barwidth, y0=smry27$mean[smry_rows], lwd=1, col=smry27$txcolor[smry_rows])
  #  arrows(x0=smry27$midx[smry_rows], y0=smry27$l95[smry_rows],y1=smry27$u95[smry_rows],lwd=1,col=smry27$txcolor[smry_rows], code=3, length=0.05, angle=90)
}


silence_message = dev.off()

brain_csf$b_kd = 1- brain_csf$b_rel
brain_csf$c_kd = 1- brain_csf$c_rel
m = lm(b_kd ~ c_kd + 0, data=brain_csf)

write(paste('Figure 4 results: rats 4 weeks post-dose each 1% CSF PrP knockdown equals ',paste0(formatC(coefficients(m)['c_kd'],digits=3),'%'),' brain PrP knockdown\n',sep=''),text_stats_path,append=T)

#  see https://stats.stackexchange.com/a/51811/153189 for discussion of lm and aov for ancova in R
# lm works but you get individual p values for each peptide, not for the overall fit of peptide:
# ancova_mrm = lm(csf_norm_lh ~ brain_norm_lh + peptide, data=brain_csf_mrm[brain_csf_mrm$peptide %in% peptides$peptide[peptides$rat=='yes'],])
# ancova_pvals = summary(ancova_mrm)$coefficients[c("peptideGENFTETDVK","peptideRPKPGGWNTGGSR","peptideVVEQMCVTQYQK","peptideYPGQGSPGGNR"),'Pr(>|t|)']
ancova_mrm = aov(csf_norm_lh ~ brain_norm_lh + peptide, data=brain_csf_mrm[brain_csf_mrm$peptide %in% peptides$peptide[peptides$rat=='yes'],])
ancova_pval = summary(ancova_mrm)[[1]]['peptide','Pr(>F)']
write(paste('Figure 4 results: ANCOVA p values testing whether MRM peptides behave differently: P',format_p(ancova_pval),'\n',sep=''),text_stats_path,append=T)







cat(file=stderr(), 'Done.\nCreating Figure S4A and C...')

resx=300
pdf('display_items/vector/figure-s4a-c.pdf',width=6.5,height=3.25)

par(mfrow=c(1,2))

path = 'data/recombinant/SHaPrP23-231-elution-figure-data-2020-02-18-MM.csv'
elution_raw = read.table(path,sep=',',header=T,quote='"',comment.char='',fill=T)

elution_uv = elution_raw[4:nrow(elution_raw),1:2]
colnames(elution_uv) = tolower(elution_raw[2,1:2])
elution_uv$min = as.numeric(elution_uv$min)
elution_uv$mau = as.numeric(elution_uv$mau)

elution_concb = elution_raw[4:nrow(elution_raw),3:4]
colnames(elution_concb) = tolower(elution_raw[2,3:4])
elution_concb$min = as.numeric(elution_concb$min)
elution_concb$percentb = as.numeric(elution_concb[,'%'])/100

elution_frac = elution_raw[3:12,5:6]
colnames(elution_frac) = tolower(elution_raw[2,5:6])
elution_frac$min = as.numeric(elution_frac$min)

flowrate = 6 # mL/min flow rate during AKTA elutions

timemax = 60
maumax = 1200
concmax = 1.05
par(mar=c(3,4,3,4))
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,maumax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:timemax, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(0,timemax,by=10), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(0,timemax,by=10), labels=seq(0,timemax,by=10), lwd=0, line=-0.5)
mtext(side=1, line=1.75, text='time (minutes)', font=1)
axis(side=2, at=seq(0,maumax,by=100), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=seq(0,maumax,by=500), lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3.0, text='A280 (mAU)', col='#0001CD', font=1)
points(elution_uv$min, elution_uv$mau, col='#0001CD', type='l', lwd=1)
par(new=T)
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,concmax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=4, at=seq(0,concmax,by=.1), labels=percent(seq(0,concmax,by=.1)), lwd=1, lwd.ticks=1, las=2)
mtext(side=4, line=3.0, text='imidazole gradient (%B)', col='#FF9912', font=1)
# the conc b curve is same for all six, so no need to loop - just plot the final one:
points(elution_concb$min, elution_concb$percentb, col='#FF9912', type='l', lwd=2)
frac_col = '#777777'
segments(x0=elution_frac$min, x1=elution_frac$min, y0=rep(0,nrow(elution_frac)), y1=rep(.075,nrow(elution_frac)), lwd=1.5, col=frac_col)
frac_mids = (elution_frac$min[1:(nrow(elution_frac)-1)] + elution_frac$min[2:(nrow(elution_frac))])/2
text(x=frac_mids, y=rep(0.05, length(frac_mids)), col=frac_col, labels=1:length(frac_mids))
text(x=40, y=.125, labels='fractions', col=frac_col)
panel = 1
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)


sec_raw = read.table('data/recombinant/SHaPrP23-232-excel-SEC-prep71-figure-data.csv',sep=',',header=T,skip=2)
colnames(sec_raw) = gsub('[^a-z0-9_]','_',tolower(colnames(sec_raw)))
sec = sec_raw[,c('min','mau')]

timemax = 111
maumax = 2000
par(mar=c(3,4,3,4))
plot(NA, NA, xlim=c(0,timemax), ylim=c(0,maumax), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:timemax, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(0,timemax,by=10), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(0,timemax,by=30), labels=seq(0,timemax,by=30), lwd=0, line=-0.5)
mtext(side=1, line=1.75, text='time (minutes)', font=1)
axis(side=2, at=seq(0,maumax,by=100), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=seq(0,maumax,by=500), lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3.0, text='A280 (mAU)', col='#0001CD', font=1)
points(sec$min, sec$mau, col='#0001CD', type='l', lwd=1)

panel = 3
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)

silence_message = dev.off()


time_lapsed = Sys.time() - start_time
cat(file=stderr(), paste0('Done.\nAll figures created in ',round(time_lapsed,1),' ',units(time_lapsed),'.\n'))











