# Reads an Open Ephys (data format version 0.4) file. Returns a vector of numeric values. 
fnFFFReadFile <- function(filename) {

# Reading Open Ephys header.
LEN_HEADER_BYTES <- 1024
  connectFile <- file(filename, 'rb')
  header_raw <- readBin(connectFile, size=1, n=LEN_HEADER_BYTES, what='raw') # 'raw' is always size 1
  close(connectFile)

# This function was written for Open Ephys data format 0.4. The OE web site
# reads "All headers are 1024 bytes long, and are written as a
# Matlab-compatible string... The current header format uses ASCII encoding,
# and defines Matlab struct with the following fields and values..."
  header_character_string <- readBin(header_raw[1:LEN_HEADER_BYTES], what='character')

# Determine the length, in bytes, of the file. Note that seek() returns the
# byte offset _prior_ to the seek.
  connectFile <- file(filename, 'rb')
  lenFileBytes <- seek(connectFile, where=0, origin='end')
  lenFileBytes <- seek(connectFile, where=0, origin='end')
  close(connectFile)

# Rather than attempting to parse the file header (a Matlab script) it would
# seem sensible instead to eyeball it, and define the requisite constants here
# by hand.
#
# Each file comprises a header, then a number of records. Each record comprises
# header, footer and data.
F_SAMP_HZ <- 30000
LEN_RECORD_HEADER_BYTES <- 12
LEN_RECORD_DATA_BYTES <- 2*1024 # for format version 0.4, 1024 2-byte data samples per record
LEN_RECORD_FOOTER_BYTES <- 10
LEN_RECORD_TOTAL_BYTES <- sum(LEN_RECORD_HEADER_BYTES, LEN_RECORD_DATA_BYTES, LEN_RECORD_FOOTER_BYTES)
  connectFile <- file(filename, 'rb')
  seekWhere <- seek(connectFile, where=LEN_HEADER_BYTES, origin='start') # move to start of record 1
  numRecordsTotal <- (lenFileBytes - LEN_HEADER_BYTES) / LEN_RECORD_TOTAL_BYTES
  signal <- rep(NA, lenFileBytes/2) # pre-allocate a little more than we need; 16-bit samples => can't be more than lenFileBytes/2
  for (iirec in 1:numRecordsTotal) {
    record <- readBin(connectFile, n=LEN_RECORD_TOTAL_BYTES, size=1, what='raw')
    recordTimestamp <- readBin(record[1:8], n=1, size=8, what='integer', signed=TRUE, endian='little')
    recordNumSamples <- readBin(record[9:10], n=1, what='integer', size=2, signed=FALSE, endian='little')
    recordRecordingNumber <- readBin(record[11:12], n=1, what='integer', size=2, signed=FALSE, endian='little') # the purpose of this field isn't immediately clear
# Everything is little-endian, except for the data samples themselves...
    ixdata <- (LEN_RECORD_HEADER_BYTES+1):(LEN_RECORD_HEADER_BYTES+2*recordNumSamples) # an index of bytes
    ixsignal <- ((iirec-1)*recordNumSamples+1):(iirec*recordNumSamples) # an index of integers
    signal[ixsignal] <- readBin(record[ixdata], n=recordNumSamples, what='integer', size=2, signed=TRUE, endian='big')
    #seekWhere <- seek(connectFile, where=0, origin='current') # debugging
  }
  close(connectFile)
  signal <- signal[!is.na(signal)]

  signal # returned
}

