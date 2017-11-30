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
  output.name <- paste(paste(file.name, "Purfi MSMS"),
                     "csv",
                     sep = ".")
  options(warn = -1)

  dir.create(file.path(getwd(),"01 Purified MSMS spectra"),
             recursive = TRUE)
  output.name <- file.path(getwd(), "01 Purified MSMS spectra", output.name)

  write.csv(result, file = output.name, row.names = F)

  # save Intermediate Data
  dir.create(path = file.path(getwd(), "00 Intermediate Data"))
  result <- list(result=result, precusor.mz=precusor.mz)
  output.name2 <- file.path(getwd(),
                            "00 Intermediate Data",
                            file.name)
  save(result, file = output.name2)
}


# msms_plot ---------------------------------------

#' @title msms_plot
#' @description  plot msms spectra of purified spectra
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spectra.name The csv file name of spectra. Default: NULL
#' @param data.name The R Dataname after purification.
#' @param precusor.mz The precusor mz. Default: NULL

msms_plot <- function(spectra.name=NULL, data.name=NULL, precusor.mz=NULL){

  # load data or read csv file
  if (is.null(spectra.name)) {
    if (!is.null(data.name)) {
      load(data.name)
      raw.data <- result[[result]]
      product.mz <- result[[precusor.mz]]
    } else {
      stop("Please input spectra file name or r data name.\n")
    }
  } else {
    raw.data <- readr::read_csv(spectra.name)
    product.mz <- raw.data$mz
  }


  if (is.null(precusor.mz)) {
    x.range <- c(0, 1.05*max(product.mz))
  } else {
    x.range <- c(0, 1.05*precusor.mz)
  }

  int.plot <- raw.data$int/max(raw.data$int)
  tag <- strsplit(x = spectra.name, split = ".csv")[[1]]
  tag.name <- paste(tag, "pdf", sep = ".")

  pdf(file = tag.name, width = 8, height = 4)
  plot(raw.data$mz, int.plot, type = "h", lwd=2, xlim = x.range, col="dodgerblue", main = tag, xlab="m/z", ylab="Intensity")
  abline(h=0, lwd=1)

  dev.off()
  # text(x = raw.data$mz, y = raw.data$int, labels = raw.data$mz)
}




# RemoveRingEffect -------------------------------------
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
