###
# Test for Granger causation channel pair by channel pair. This runs after FFF_analyze.R.
#####

# Constants, functions...
library('lmtest')
ORDER_MODEL <- 10
fnSquare <- function(ii,jj,sz,brd) polygon(x=ii + c(-sqrt(sz)/2,-sqrt(sz)/2,sqrt(sz)/2,sqrt(sz)/2),
                                           y=jj + c(-sqrt(sz)/2,sqrt(sz)/2,sqrt(sz)/2,-sqrt(sz)/2),
                                           col='gray', border=brd)


# Whiten -- remove the trial-triggered average, and divide by the s.d. across trials.
trace_per_trial_per_channel_zscore <- 0*trace_per_trial_per_channel
for (iichan in 1:dim(trace_per_trial_per_channel)[3]) trace_per_trial_per_channel_zscore[,,iichan] <- (trace_per_trial_per_channel[,,iichan] - apply(trace_per_trial_per_channel[,,iichan], 1, mean)) / apply(trace_per_trial_per_channel[,,iichan], 1, sd)

# The Granger test returns an F statistic which, here, we store in a matrix (fmat) marking pair-wise 'causal' relationships between channels.
fmat <- matrix(0, ncol=length(MAP), nrow=length(MAP))
for (iichan in 1:length(MAP)) {
  for (jjchan in 1:length(MAP)) {
    if (iichan != jjchan) {
message(sprintf('Granger testing channels %d and %d...',iichan,jjchan))
      f_per_trial <- 1:NUM_TRIALS
      for (iitrial in 1:NUM_TRIALS) {
        gt <- grangertest(trace_per_trial_per_channel_zscore[,iitrial,iichan], trace_per_trial_per_channel_zscore[,iitrial,jjchan], order=ORDER_MODEL)
        f_per_trial[iitrial] <- gt$F[2]
      }
      fmat[iichan,jjchan] <- mean(f_per_trial)
    }
  }
}
# Normalize the fmat.
diag(fmat) <- max(fmat, na.rm=TRUE)
frange <- range(fmat)
fmat <- fmat - min(fmat)
fmat <- fmat / max(fmat)
diag(fmat) <- 0

# Plotting overhead.
filename_granger_pdf <- sprintf('%s%s.granger.pdf',DIR_PLOT,rev(strsplit(DIR_DATA,'/',fixed=TRUE)[[1]])[1])
pdf(file=filename_granger_pdf,width=8.5*length(MAP)/32,height=8.5*length(MAP)/32,paper="special",colormodel="srgb",compress=TRUE)
par(mar=c(3,3,3,3), ps=10)

plot(NA, xlim=c(0,length(MAP)+1), ylim=c(0,length(MAP)+1), main='', xlab='', ylab='', axes=FALSE, asp=1)
for (ii in 1:length(MAP)) for (jj in 1:length(MAP)) fnSquare(jj,length(MAP)+1-ii,fmat[ii,jj],gray(0.5))
for (ii in 1:length(MAP)) text(x=ii, y=length(MAP)+0.5, labels=sprintf('%d',ii), font=3, xpd=TRUE)
for (ii in 1:length(MAP)) text(x=-1, y=length(MAP)+1-ii, labels=sprintf('%d',ii), font=3, xpd=TRUE)
mtext('X', side=2, line=1, font=3)
mtext('Y', side=3, line=1, font=3)
mtext(sprintf('Granger test that X causes Y; model order = %d; range of statistic %.2f to %.2f',ORDER_MODEL,frange[1],frange[2]), side=3, line=2, font=3)

# More plotting overhead.
dev.off()
message(sprintf('Wrote %s',filename_granger_pdf))
#####


