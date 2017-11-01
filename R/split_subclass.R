#' @title split_subclass
#' @description  split a big library into subclasses
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param raw.data The csv file name of raw data.
#' @param raw.MD The csv file name of Molecular Descriptor.
#' @param bulk.std.name The bulk structure name of each classes.
#'
#' @example
#' bulk.std.name <- c("PC(", "PC(O-", "PC(P-", "PC(dO-", "LPC(", "LPC(O-", "LPC(P-")
#' split_subclass(raw.data = "result.csv", raw.MD="Define-descriptors by rcdk.csv", bulk.std.name=bulk.std.name)

split_subclass <- function(raw.data="result.csv",
                           raw.MD="Define-descriptors by rcdk.csv",
                           bulk.std.name){
  raw.data <- readr::read_csv(file = raw.data)
  raw.MD <- readr::read_csv(file = raw.MD)

  load(system.file("list", "lipidmaps_info (171101)", package = "LipLibBuild"))

  tab <- sort(unique(raw.data$sub.class))
  idx.subclass <- lapply(tab, function(x){
    temp <- which(raw.data$sub.class==x)
    return(temp)
  })

  # bulk.std.name <- c("PC(", "PC(O-", "PC(P-", "PC(dO-", "LPC(", "LPC(O-", "LPC(P-")

  lapply(seq(length(idx.subclass)), function(i){
    temp.idx <- idx.subclass[[i]]
    temp.data <- raw.data[temp.idx,]
    temp.MD <- raw.data[temp.idx,]

    temp.num <- formatC(seq(length(temp.idx)), flag = '0', width = 6, mode = "integer")
    temp.num <- paste(temp.data$sub.class, temp.num, sep = "")

    temp.nc <- temp.data$sn1.chain.length+temp.data$sn2.chain.length
    temp.db <- temp.data$sn1.double.bonds+temp.data$sn2.double.bonds
    temp <- paste(temp.nc, temp.db, sep = ":")
    bulk.str <- paste(paste(bulk.std.name[i], temp, sep = ""), ")", sep = "")
    lipid.map.id <- sapply(seq(nrow(temp.data)), function(i){
      temp.id <- which(lipidmaps.info$abbr.name==temp.data$abbr.name[i])
      if (length(temp.id) == 0) {
        temp.id <- ""
      } else {
        temp.id <- lipidmaps.info$lipidmaps.id[temp.id]
      }
      return(temp.id)
    })

    # lipid.map.id <- merge(x = temp.data, y = lipidmaps.info, by = "abbr.name", all.x = T)
    # lipid.map.id <- gsub(pattern = NA, replacement = "", lipid.map.id$lipidmaps.id)

    # lipid.map.id <- do.call(c, lipid.map.id)

    result1 <- data.frame(id=temp.num,
                          temp.data,
                          bulk.structure=bulk.str,
                          lipid.map.id=lipid.map.id,
                          stringsAsFactors = F)

    result2 <- data.frame(id=temp.num,
                          temp.MD,
                          stringsAsFactors = F)

    address <- paste(getwd(), "split_library", sep = "/")
    dir.create(path = address, recursive = T)

    temp <- temp.data$sub.class[1]
    file.name1 <- paste(address, paste(paste(temp, "Info", sep = "_"), "csv", sep = "."), sep = "/")
    file.name2 <- paste(address, paste(paste(temp, "MD", sep = "_"), "csv", sep = "."), sep = "/")

    write.csv(result1, file = file.name1, row.names = F)
    write.csv(result2, file = file.name2, row.names = F)
  })
}
