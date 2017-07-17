# Given a 30kHz signal (probably) of sync pulses, returns the indices of pulses.
fnFFFProcSyncSignal <- function(signal_sync, f_samp_hz=30000, verbose=TRUE) {

  signal_sync_logical <- signal_sync > max(signal_sync)/2
  dsignal_sync_logical <- diff(signal_sync_logical)
  ixpulse <- which(dsignal_sync_logical == 1)
  if (verbose) message(sprintf('Found %d pulses spanning %.2f seconds', length(ixpulse), (rev(ixpulse)[1] - ixpulse[1]) / f_samp_hz))
  if (verbose) message(sprintf('Pulse 1 occurs at %.2f seconds and the last pulse at %.2f seconds', ixpulse[1]/f_samp_hz, rev(ixpulse)[1]/f_samp_hz))

  ixpulse # returned
}

