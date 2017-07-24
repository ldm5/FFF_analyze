###
# Test for Granger causation between channels. This runs after FFF_analyze.R.
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

# The Granger test returns an F statistic which, here, we store in a matrix marking pair-wise relationships between channels, fmat.
fmat <- matrix(0, ncol=32, nrow=32)
for (ii in 1:32) {
  for (jj in 1:32) {
    if (ii != jj) {
message(sprintf('Granger testing channels %d and %d...',ii,jj))
      p_per_trial <- f_per_trial <- 1:NUM_TRIALS
      for (iitrial in 1:NUM_TRIALS) {
        gt <- grangertest(trace_per_trial_per_channel_zscore[,iitrial,ii], trace_per_trial_per_channel_zscore[,iitrial,jj], order=ORDER_MODEL)
        f_per_trial[iitrial] <- gt$F[2]
        p_per_trial[iitrial] <- gt$Pr[2]
      }
      fmat[ii,jj] <- mean(f_per_trial)
    }
  }
}
# Normalize the fmat; you'll find NAs only on the main diagonal.
diag(fmat) <- max(fmat, na.rm=TRUE)
fmat <- fmat - min(fmat)
fmat <- fmat / max(fmat)
diag(fmat) <- 0

filename_granger_pdf <- sprintf('%s%s.granger.pdf',DIR_PLOT,rev(strsplit(DIR_DATA,'/',fixed=TRUE)[[1]])[1])
pdf(file=filename_granger_pdf,width=8.5,height=8.5,paper="special",colormodel="srgb",compress=TRUE)
par(mar=c(3,3,3,3), ps=10)
plot(NA, xlim=c(0,33), ylim=c(0,33), main='', xlab='', ylab='', axes=FALSE, asp=1)
for (ii in 1:32) for (jj in 1:32) fnSquare(jj,33-ii,fmat[ii,jj],gray(0))
for (ii in 1:32) text(x=ii, y=33.5, labels=sprintf('%d',ii), font=3, xpd=TRUE)
for (ii in 1:32) text(x=-1, y=33-ii, labels=sprintf('%d',ii), font=3, xpd=TRUE)
mtext(sprintf('Granger test that X causes Y, model order = %d',ORDER_MODEL), side=3, line=2)
mtext('X', side=2, line=1)
mtext('Y', side=3, line=1)
mtext(sprintf('Granger test that X causes Y, model order = %d',ORDER_MODEL), side=3, line=2)

dev.off()
message(sprintf('Wrote %s',filename_granger_pdf))
#####

