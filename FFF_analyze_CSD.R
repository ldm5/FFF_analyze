###
# Compute CSD, taking 2nd derivative over channels. Run this after FFF_analyze.R.
##########

# Constants, functions...
KERNEL <- c(1,-2,1)
fnDiff <- function(x,kernel) convolve(x, kernel, type='filter')

# Convolve our double-derivative kernel over channels at each point in time.
ddtrial_triggered_average_per_channel <- numeric()
for (ix in 1:length(fnDownsampleV(ixshort_time, FAC_DOWNSAMPLE))) ddtrial_triggered_average_per_channel <- cbind(ddtrial_triggered_average_per_channel, fnDiff(trial_triggered_average_per_channel[ix,], KERNEL))

# Plot overhead.
filename_csd_pdf <- sprintf('%s%s.csd_z.pdf',DIR_PLOT,rev(strsplit(DIR_DATA,'/',fixed=TRUE)[[1]])[1])
pdf(file=filename_csd_pdf,width=8.5,height=8.5,paper="special",colormodel="srgb",compress=TRUE)
par(mar=c(3,3,3,3), ps=10)

sc <- 2*max(abs(ddtrial_triggered_average_per_channel))
# Downsample for representation as shaded polygons...
plot(NA, xlim=range(time_trial_s[ixshort_time]) , ylim=c(-length(MAP),1)*sc, main='', axes=F, xlab='', ylab='')
for (iitrace in 1:nrow(ddtrial_triggered_average_per_channel)) {
  yoffset <- -(iitrace-1+floor(length(KERNEL)/2))*sc
  polygon(fnPolybracket(fnDownsampleV(time_trial_s[ixshort_time], FAC_DOWNSAMPLE)), yoffset + c(0,pmax(0,ddtrial_triggered_average_per_channel[iitrace,]),0), col='red', border=NA)
  polygon(fnPolybracket(fnDownsampleV(time_trial_s[ixshort_time], FAC_DOWNSAMPLE)), yoffset + c(0,pmin(0,ddtrial_triggered_average_per_channel[iitrace,]),0), col='blue', border=NA)
  mtext(sprintf('%d',iitrace+1), 2, at=yoffset, cex=0.5, las=1, font=3)
}
lines(rep(TIME_STIM_ON_S, 2), range(c(-length(MAP),1)*sc), lty='dashed', col='black')
axis(1, at=c(TIME_TRACE_START_S,TIME_TRACE_START_S+0.05), labels=FALSE)
mtext('50 ms', 1, at=TIME_TRACE_START_S+0.05/2, line=1, font=3, cex=0.8)
mtext('Electrode number', 2, line=3, font=3, cex=0.8)

# More plot overhead.
dev.off()
print(sprintf('Wrote %s',filename_csd_pdf))
###########

