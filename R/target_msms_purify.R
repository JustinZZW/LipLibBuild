#' @title target_msms_purify
#' @description  Purify msms spectra of targets
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param file.name The csv file name of spectra.
#' @param abs.min.count The minimnum value of absolute intensity. Default: 100.
#' @param rela.min.count The minimnum value of relative intensity (1000 as benchmark). Default: 1.


target_msms_purify <- function(file.name=NULL,
                               abs.min.count=100,
                               rela.min.count=1){
  if (is.null(file.name)) {
    stop("Please input csv file name of spectra.\n")
  }
  raw.data <- readLines(file.name)

  # extract target mz
  info <- raw.data[1]
  target.mass <- substr(info,
                        start = regexpr(text = info, pattern = "@")+7,
                        stop = regexpr(text = info, pattern = "\\[")-1)
  target.mass <- as.numeric(target.mass)
  raw.data <- raw.data[-c(1:2)]

  # extract msms information
  final.data <- lapply(raw.data, function(x){
    temp <- unlist(strsplit(x = x, split = ","))
    mz <- round(as.numeric(temp[2]), digits = 4)
    int <- round(as.numeric(temp[3]), digits = 0)
    result <- c(mz, int)
    return(result)
  })

  final.data <- do.call(rbind, final.data)
  final.data <- data.frame(mz=as.numeric(final.data[,1]),
                           int=as.numeric(final.data[,2]),
                           stringsAsFactors = F)

  # purify msms
  # remove mz larger than targeted mz
  temp.idx <- which(final.data$mz <= (target.mass+0.5))
  final.data <- final.data[temp.idx,]

  # remove mz with abs intensity less than cutoff
  temp.idx <- which(final.data$int >= abs.min.count)
  final.data <- final.data[temp.idx,]

  # remove mz with relative intensity less than cutoff
  final.data$int <- final.data$int/max(final.data$int)*1000
  temp.idx <- which(final.data$int >= rela.min.count)
  final.data <- final.data[temp.idx,]

  # remove ring effect
  result <- RemoveRingEffect(spec = final.data,
                             mz.diff.thr = 0.3,
                             int.rel.thr = 0.2)

  file.name <- strsplit(x = file.name, split = ".cs")[[1]][1]
  file.name <- paste(paste(file.name, "Purfi MSMS"),
                     "csv",
                     sep = ".")

  dir.create(file.path(getwd(),"01 Purified MSMS spectra"),
             recursive = TRUE)
  file.name <- file.path(getwd(), "01 Purified MSMS spectra", file.name)


  write.csv(result, file = file.name, row.names = F)
}

RemoveRingEffect <- function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
  nr.ring <- nrow(spec) + 1
  mz <- spec[, 'mz']

  mz.diff <- diff(mz)
  idx.mzdiff <- which(mz.diff <= mz.diff.thr)
  if (length(idx.mzdiff) == 0) {
    return(spec)
  }

  nr.ring.possible <- unique(c(idx.mzdiff, idx.mzdiff + 1))

  # remove ringeffect loop
  while (TRUE) {

    idx.int.max <- which.max(spec[nr.ring.possible, 2])
    nr.int.max <- nr.ring.possible[idx.int.max] # the index of possible Ringeffect ions with maxium intensity
    int.thr <- spec[nr.int.max, 2] * int.rel.thr # the threshold = 0.2*max.int (possible ring)

    mz.diff <- abs(mz[nr.ring.possible[-idx.int.max]] - mz[nr.int.max])
    int <- spec[nr.ring.possible[-idx.int.max], 2]
    nr.ring <- append(nr.ring, nr.ring.possible[-idx.int.max][which(mz.diff <= mz.diff.thr & int <= int.thr)])
    nr.ring.possible <- nr.ring.possible[!nr.ring.possible %in% c(nr.ring, nr.int.max)]
    if (length(nr.ring.possible) == 0) {
      break # break loop untill satisfy the nr.ring.possible==0
    }
  }

  return(spec[-nr.ring, , drop = FALSE])
}
