#' @title Lipid library build
#' @description  Build lipid library based on sdf file and smile file
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param sdf.file The file name of sdf file
#' @param smile.file The file name of smile file
#' @param is.output The output csv file. The default is TRUE.
#'
#' @example LipLibBuild(sdf.file="PC.sdf", smile.file="PC.smile")

LipLibBuild <- function(sdf.file,
                        smile.file,
                        is.output=TRUE){
  raw.data <- readr::read_lines("GPAbbrev.sdf")
  raw.smiles <- readr::read_lines("GPAbbrev.sdf.smiles")

  idx.abbr <- which(raw.data==">  <Abbrev>")+1
  idx.category <- which(raw.data==">  <LM Category>")+1
  idx.main.class <- which(raw.data==">  <LM Main Class>")+1
  idx.sub.class <- which(raw.data==">  <LM Sub Class>")+1
  idx.sn1.chain.length <- which(raw.data==">  <Sn1 Chain Length>")+1
  idx.sn1.double.bonds <- which(raw.data==">  <Sn1 Double Bonds>")+1
  idx.sn2.chain.length <- which(raw.data==">  <Sn2 Chain Length>")+1
  idx.sn2.double.bonds <- which(raw.data==">  <Sn2 Double Bonds>")+1
  idx.sys.name <- which(raw.data==">  <Systematic Name>")+1

  abbr.name <- as.character(raw.data[idx.abbr])
  category <- as.character(raw.data[idx.category])
  main.class <- as.character(raw.data[idx.main.class])
  sub.class <- as.character(raw.data[idx.sub.class])
  sn1.chain.length <- as.numeric(raw.data[idx.sn1.chain.length])
  sn1.double.bonds <- as.numeric(raw.data[idx.sn1.double.bonds])
  sn2.chain.length <- as.numeric(raw.data[idx.sn2.chain.length])
  sn2.double.bonds <- as.numeric(raw.data[idx.sn2.double.bonds])
  sys.name <- as.character(raw.data[idx.sys.name])

  result <- data.frame(abbr.name=abbr.name,
                       smiles=raw.smiles,
                       category=category,
                       main.class=main.class,
                       sub.class=sub.class,
                       sn1.chain.length=sn1.chain.length,
                       sn1.double.bonds=sn1.double.bonds,
                       sn2.chain.length=sn2.chain.length,
                       sn2.double.bonds=sn2.double.bonds,
                       sys.name=sys.name,
                       stringsAsFactors = F)

  if (is.output==TRUE) {
    temp <- paste(getwd(), "library", sep = "/")
    dir.create(path = temp)
    temp <- paste(temp, "lipid_lib.csv", sep = "/")
    write.csv(result, file = temp, row.names = F)
  }
}
