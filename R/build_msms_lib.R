#' @title build_msms_lib
#' @author Zhiwei Zhou
#' \email{zhozw@@sioc.ac.cn}
#' @param dir.path directory path
#' @param adduct the adduct type. It supports "[M+H]", "[M+Na]", "[M+NH4]", "[M-H]", "[M+HCOO]"

build_msms_lib <- function(dir.path,
                           adduct=c("[M+H]", "[M+Na]", "[M+NH4]",
                                    "[M-H]", "[M+HCOO]")){
  adduct <- match.arg(adduct)

  file.names <- dir(dir.path)
  frag.template <- file.names[grep(x = file.names, pattern = "fragment template")]
  int.template <- file.names[grep(x = file.names, pattern = "intensity template")]

  frag.template <- readxl::read_xls(frag.template)
  col.list <- colnames(frag.template)
  dir.create(file.path(".", "00_result"), recursive = T)

  int.template <- readr::read_csv(int.template)

  if (adduct=="[M+H]") {
    temp.idx <- grep(pattern = "HP", x = col.list)
    if (length(temp.idx)>0){
      temp <- which(frag.template[,1]=="exactmass")
      temp.data <- frag.template[-c(1:temp),]
      names(temp.data) <- frag.template[temp,]
      temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
      temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                              formula=temp.data$Formula)

      cor.idx <- match(colnames(temp.data2), int.template$annotation)
      int.template <- int.template[cor.idx,]

      result <- lapply(seq(ncol(temp.data2)), function(i){
        temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
        int.temp <- int.template$int[i]
        temp <- data.frame(temp, rela.int=int.temp, stringsAsFactors = F)
        return(temp)
      })

      result <- do.call(cbind, result)
      colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

      result <- cbind(temp.data, result)


      output.name <- file.path(".", "00_result", "[M+H].csv")
      write.csv(result, file = output.name, row.names = F)

    } else {
      cat("The adduct type is inconsistent\n\n")
    }
  }


  if (adduct=="[M+Na]") {
    temp.idx <- grep(pattern = "NaP", x = col.list)
    if (length(temp.idx)>0){
      temp <- which(frag.template[,1]=="exactmass")
      temp.data <- frag.template[-c(1:temp),]
      names(temp.data) <- frag.template[temp,]
      temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
      temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                              formula=temp.data$Formula)

      cor.idx <- match(colnames(temp.data2), int.template$annotation)
      int.template <- int.template[cor.idx,]

      result <- lapply(seq(ncol(temp.data2)), function(i){
        temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
        int.temp <- int.template$int[i]
        temp <- data.frame(temp, rela.int=int.temp, stringsAsFactors = F)
        return(temp)
      })

      result <- do.call(cbind, result)
      colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

      result <- cbind(temp.data, result)


      output.name <- file.path(".", "00_result", "[M+Na].csv")
      write.csv(result, file = output.name, row.names = F)

    } else {
      cat("The adduct type is inconsistent\n\n")
    }
  }

  if (adduct=="[M+NH4]") {
    temp.idx <- grep(pattern = "NH4P", x = col.list)
    if (length(temp.idx)>0){
      temp <- which(frag.template[,1]=="exactmass")
      temp.data <- frag.template[-c(1:temp),]
      names(temp.data) <- frag.template[temp,]
      temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
      temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                              formula=temp.data$Formula)

      cor.idx <- match(colnames(temp.data2), int.template$annotation)
      int.template <- int.template[cor.idx,]

      result <- lapply(seq(ncol(temp.data2)), function(i){
        temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
        int.temp <- int.template$int[i]
        temp <- data.frame(temp, rela.int=int.temp, stringsAsFactors = F)
        return(temp)
      })

      result <- do.call(cbind, result)
      colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

      result <- cbind(temp.data, result)


      output.name <- file.path(".", "00_result", "[M+NH4].csv")
      write.csv(result, file = output.name, row.names = F)

    } else {
      cat("The adduct type is inconsistent\n\n")
    }
  }


  if (adduct=="[M-H]") {
    temp.idx <- grep(pattern = "HN", x = col.list)
    if (length(temp.idx)>0){
      temp <- which(frag.template[,1]=="exactmass")
      temp.data <- frag.template[-c(1:temp),]
      names(temp.data) <- frag.template[temp,]
      temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
      temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                              formula=temp.data$Formula)

      cor.idx <- match(colnames(temp.data2), int.template$annotation)
      int.template <- int.template[cor.idx,]

      result <- lapply(seq(ncol(temp.data2)), function(i){
        temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
        int.temp <- int.template$int[i]
        temp <- data.frame(temp, rela.int=int.temp, stringsAsFactors = F)
        return(temp)
      })

      result <- do.call(cbind, result)
      colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

      result <- cbind(temp.data, result)


      output.name <- file.path(".", "00_result", "[M-H].csv")
      write.csv(result, file = output.name, row.names = F)

    } else {
      cat("The adduct type is inconsistent\n\n")
    }
  }


  if (adduct=="[M+HCOO]") {
    temp.idx <- grep(pattern = "HCOON", x = col.list)
    if (length(temp.idx)>0){
      temp <- which(frag.template[,1]=="exactmass")
      temp.data <- frag.template[-c(1:temp),]
      names(temp.data) <- frag.template[temp,]
      temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
      temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                              formula=temp.data$Formula)

      cor.idx <- match(colnames(temp.data2), int.template$annotation)
      int.template <- int.template[cor.idx,]

      result <- lapply(seq(ncol(temp.data2)), function(i){
        temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
        int.temp <- int.template$int[i]
        temp <- data.frame(temp, rela.int=int.temp, stringsAsFactors = F)
        return(temp)
      })

      result <- do.call(cbind, result)
      colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

      result <- cbind(temp.data, result)


      output.name <- file.path(".", "00_result", "[M+HCOO].csv")
      write.csv(result, file = output.name, row.names = F)

    } else {
      cat("The adduct type is inconsistent\n\n")
    }
  }

}
