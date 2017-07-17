###
# This block contains all the constants you might need to configure online.
##########
DIR_DATA <- '/Users/luke/physiolArchive/m664/m664r2#28_2017-03-09_13-36-38_full_field_flash/' # keep trailing slash
DIR_PLOT <- '/tmp/' # keep trailing slash
DIR_R_BIN <- '/Users/luke/physiolArchive/bin/' # keep trailing slash
FLAG_NOTCH_60 <- TRUE
NUM_TRIALS <- 10
##########

###
# Libraries, functions, etc.
##########
library('signal')
source(sprintf('%sRFFFReadFile.R', DIR_R_BIN))
source(sprintf('%sRFFFProcSyncSignal.R', DIR_R_BIN))
fnDownsampleV <- function(x,f) x[seq(from=1, to=length(x), by=f)] # downsample a vector
fnGetMonkeyStringExp <- function(string_in) regmatches(string_in, regexpr("m[0-9]{3}[rl][0-9]{1,3}#[0-9]{1,3}",string_in))
fnDDGetMonkeyString <- function(string_in) regmatches(string_in, regexpr("m[0-9]{3}[rl][0-9]{1,3}",string_in))
##########

###
# Get sync pulses.
##########
FILENAME_SYNC <- sprintf('%s%s',DIR_DATA,'100_ADC1.continuous')
signal_sync <- fnFFFReadFile(FILENAME_SYNC)
ixpulse <- fnFFFProcSyncSignal(signal_sync)
##########

###
# Set-up for all parts of this analysis.
##########
ID_MONK <- substr(fnDDGetMonkeyString(DIR_DATA), start=1, stop=4)
MAP <- c(1,17,16,32,3,19,14,30,9,25,10,20,8,24,2,29,7,26,15,21,11,23,12,28,6,18,13,22,5,27,4,31) # m661, 2, 3, 4
F_SAMP_HZ <- 30000
DURATION_TRIAL_PULSES <- 360
DURATION_TRIAL_S <- 3
TIME_STIM_ON_S <- 0.25
TIME_STIM_OFF_S <- 1.25
TIME_TRACE_START_S <- 0.15
TIME_TRACE_END_S <- 0.5
duration_trial_samples <- DURATION_TRIAL_S*F_SAMP_HZ
time_trial_s <- (0:(duration_trial_samples-1))/F_SAMP_HZ
ixshort_time <- which((time_trial_s >= TIME_TRACE_START_S) & (time_trial_s < TIME_TRACE_END_S))
ixtrial <- ixpulse[seq(from=1, to=length(ixpulse), by=DURATION_TRIAL_PULSES)]
FAC_DOWNSAMPLE <- 100
##########

###
# Collect traces, per trial and electrode. 
##########
trace_per_trial_per_channel <- array(NA, dim=c(length(fnDownsampleV(1:length(ixshort_time), FAC_DOWNSAMPLE)),NUM_TRIALS,length(MAP)))
flag_artifact_per_trial_per_channel <- matrix(TRUE, nrow=length(MAP), ncol=NUM_TRIALS)
trial_triggered_average_per_channel <- matrix(NA, ncol=length(MAP), nrow=length(fnDownsampleV(1:length(ixshort_time), FAC_DOWNSAMPLE)))
for (iich in 1:length(MAP)) {
message(sprintf('Processing channel %d (%d of %d)...',MAP[iich],iich,length(MAP)))

  filename_channel <- sprintf('%s%s',DIR_DATA,sprintf('100_CH%d.continuous',MAP[iich]))
  signal <- fnFFFReadFile(filename_channel)

# Low-pass filter. Compute trial-triggered average.
#signal_test <- 0*signal # debug
#signal_test[ixtrial + TIME_STIM_ON_S*F_SAMP_HZ] <- max(signal) # debug
#signal <- signal_test
  bflow <- butter(3, 100/(F_SAMP_HZ/2), type='low') # 3rd-order Butterworth; roll-off is 18dB per octave
  signal_filtered <- signal::filter(bflow, signal)

  if (FLAG_NOTCH_60) {
    bfnotch <- butter(3, c(50,70)/(F_SAMP_HZ/2), type='stop')
    signal_filtered <- signal::filter(bfnotch, signal_filtered)
  }

# Detect artifacts.
  trace_power_per_trial <- numeric()
  for (iitrial in 1:NUM_TRIALS) {
    this_sig <- signal[ixtrial[iitrial]:(ixtrial[iitrial] + duration_trial_samples-1)]
    this_trace <- this_sig[ixshort_time]
    trace_power_per_trial <- c(trace_power_per_trial,sum(this_trace^2))
  }
  ix_artifact <- which(trace_power_per_trial > quantile(trace_power_per_trial, probs=0.9))

# Collect traces per trial.
  for (iitrial in 1:NUM_TRIALS) {
    this_sig <- signal_filtered[ixtrial[iitrial]:(ixtrial[iitrial] + duration_trial_samples-1)]
    this_trace <- this_sig[ixshort_time]
    trace_per_trial_per_channel[,iitrial,iich] <- fnDownsampleV(this_trace, FAC_DOWNSAMPLE)
  }

  ix_nonartifact <- setdiff(1:NUM_TRIALS,ix_artifact)
  flag_artifact_per_trial_per_channel[iich,ix_nonartifact] <- FALSE
  trial_triggered_average_per_channel[,iich] <- rowMeans(trace_per_trial_per_channel[,ix_nonartifact,iich])
}
##########

##
# 'Per trial' plot of LFPs.
#####
filename_pdf <- sprintf('%s%s.per_trial.pdf',DIR_PLOT,rev(strsplit(DIR_DATA,'/',fixed=TRUE)[[1]])[1])
pdf(file=filename_pdf, width=2*(NUM_TRIALS+1), height=2*length(MAP), paper='special', colormodel='srgb', compress=TRUE)
par(ps=10, mfrow=c(length(MAP),NUM_TRIALS+1), mai=rep(0.1,4), pin=c(2,2))

F_SAMP_HZ_ <- F_SAMP_HZ / FAC_DOWNSAMPLE
ixstim_on <- F_SAMP_HZ_*(TIME_STIM_ON_S - TIME_TRACE_START_S)
sc <- max(abs(trace_per_trial_per_channel))
for (iich in 1:length(MAP)) {
  these_rms60 <- numeric()
  for (iitrial in 1:NUM_TRIALS) {

 # Get this signal, and the amplitude of its line noise.
    this_trace <- trace_per_trial_per_channel[,iitrial,iich]

# If artifact trial, paint it red.
    this_col <- 'black'
    if (as.logical(flag_artifact_per_trial_per_channel[iich,iitrial])) this_col <- 'red'

# Plot it...
    plot(NA, xlim=c(1,length(ixshort_time)/FAC_DOWNSAMPLE), ylim=sc*c(-1,1), xlab='', ylab='', axes=FALSE, main='')
    lines(rep(ixstim_on,2), sc*c(-0.5,0.5), col='gray', lty='dashed')
    lines(c(1,ixstim_on/2), rep(-0.75*sc,2), col='black', lty='solid')
    lines(this_trace, col=this_col)
# Definition: sine wave's r.m.s. = Vpp/(2*sqrt(2))
    these_rms60 <- c(these_rms60,4*Mod(fft(this_trace)[22]/length(this_trace))/(2*sqrt(2)))
    text(x=1, y=sc, labels=sprintf('%d %d',iich,iitrial), font=3, pos=4, offset=0, col='blue')
    text(x=length(ixshort_time)/FAC_DOWNSAMPLE, y=sc, labels=sprintf('%.2f',rev(these_rms60)[1]), font=3, pos=2, offset=0, col='black')
    text(x=ixstim_on/4, y=-0.75*sc, labels=sprintf('%d ms',ceiling(1000*ixstim_on/2/F_SAMP_HZ_)), font=3, pos=1, col='black')
  }

# Plot the mean over trials for this channel...
  this_mean <- trial_triggered_average_per_channel[,iich]
  plot(NA, xlim=c(1,length(ixshort_time)/FAC_DOWNSAMPLE), ylim=sc*c(-1,1), xlab='', ylab='', axes=FALSE, main='')
  lines(rep(ixstim_on,2), sc*c(-0.5,0.5), col='gray', lty='dashed')
  lines(c(1,ixstim_on/2), rep(-0.75*sc,2), col='black', lty='solid')
  lines(this_mean, col='black', lwd=3)
  text(x=1, y=sc, labels=sprintf('%d mean',iich), font=3, pos=4, offset=0, col='blue')
  text(x=length(ixshort_time)/FAC_DOWNSAMPLE, y=sc, labels=sprintf('%.2f, %.2f',sqrt(mean(this_mean^2)),mean(these_rms60)), font=3, pos=2, offset=0, col='black')
  text(x=ixstim_on/4, y=-0.75*sc, labels=sprintf('%d ms',ceiling(1000*ixstim_on/2/F_SAMP_HZ_)), font=3, pos=1, col='black')
}

dev.off()
message(sprintf('Wrote %s',filename_pdf))
#####

###
# Downsample, represent field potentials as shaded polygons.  Note that the
# trial_triggered_average_per_channel is already downsampled, but not ixshort_time.
#####
filename_pdf <- sprintf('%s%s.lfp_z.pdf',DIR_PLOT,rev(strsplit(DIR_DATA,'/',fixed=TRUE)[[1]])[1])
pdf(file=filename_pdf,width=8.5,height=8.5,paper="special",colormodel="srgb",compress=TRUE)
par(mar=c(3,3,3,3), ps=10)

fnPolybracket <- function(x) c(x[1],x,rev(x)[1])
sc <- max(abs(trial_triggered_average_per_channel))
plot(NA, xlim=range(time_trial_s[ixshort_time]) , ylim=c(-length(MAP),1)*sc, main='', axes=F, xlab='', ylab='')
for (iitrace in 1:ncol(trial_triggered_average_per_channel)) {
  yoffset <- -(iitrace-1)*sc
  mtext(sprintf('%d',iitrace), 2, at=yoffset, cex=0.8, las=1, font=3)
  polygon(fnPolybracket(fnDownsampleV(time_trial_s[ixshort_time], FAC_DOWNSAMPLE)), yoffset + c(0,pmax(0,trial_triggered_average_per_channel[,iitrace]),0), col='red', border=NA)
  polygon(fnPolybracket(fnDownsampleV(time_trial_s[ixshort_time], FAC_DOWNSAMPLE)), yoffset + c(0,pmin(0,trial_triggered_average_per_channel[,iitrace]),0), col='blue', border=NA)
}
lines(rep(TIME_STIM_ON_S, 2), range(c(-length(MAP),1)*sc), lty='dashed', col='black')
lines(rep(TIME_STIM_OFF_S, 2), range(c(-length(MAP),1)*sc), lty='dashed', col='black')
axis(1, at=c(TIME_TRACE_START_S,TIME_TRACE_START_S+0.05), labels=FALSE)
mtext('50 ms', 1, at=TIME_TRACE_START_S+0.05/2, line=1, font=3, cex=0.8)
mtext('Electrode number', 2, line=2, font=3, cex=0.8)

dev.off()
message(sprintf('Wrote %s',filename_pdf))
##########
