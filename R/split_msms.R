#' @title split_msms
#' @description split msms fragments into independent csv file according to TJ lipid library. The template should be labeled as "HP"="M+H"; "NaP"="M+Na"; "NH4P"="M+NH4"; "HN"="M-H"; "HCOON"="M+HCOO". The number "1" represents the start column; number "2" represents the end column.
#' @author Zhiwei Zhou
#' {zhouzw@@sioc.ac.cn}
#' @param file.name The csv file name for split.

split_msms <- function(file.name){
  raw.data <- readxl::read_xlsx(file.name)
  col.list <- colnames(raw.data)
  file.address <- strsplit(x = file.name, split = "\\.")[[1]][1]
  file.address <- paste(getwd(), "split_msms_library", file.address, sep = "/")
  dir.create(file.address, recursive = T)

  temp.idx <- grep(pattern = "HP", x = col.list)
  if (length(temp.idx)>0){
    temp <- which(raw.data[,1]=="exactmass")
    temp.data <- raw.data[-c(1:temp),]
    names(temp.data) <- raw.data[temp,]
    temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
    temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                            formula=temp.data$Formula)

    result <- lapply(seq(ncol(temp.data2)), function(i){
      temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
      temp <- data.frame(temp, rela.int="", stringsAsFactors = F)
      return(temp)
    })

    result <- do.call(cbind.data.frame, result)
    colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

    result <- data.frame(temp.data, result, stringsAsFactors = F)

    temp <- paste(file.address, "[M+H].csv", sep = "/")
    write.csv(result, temp, row.names = F)
  }

  #--------------------------------------------------
  temp.idx <- grep(pattern = "NaP", x = col.list)
  if (length(temp.idx)>0){
    temp <- which(raw.data[,1]=="exactmass")
    temp.data <- raw.data[-c(1:temp),]
    names(temp.data) <- raw.data[temp,]
    temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
    temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                            formula=temp.data$Formula)

    result <- lapply(seq(ncol(temp.data2)), function(i){
      temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
      temp <- data.frame(temp, rela.int="", stringsAsFactors = F)
      return(temp)
    })

    result <- do.call(cbind.data.frame, result)
    colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

    result <- data.frame(temp.data, result, stringsAsFactors = F)

    temp <- paste(file.address, "[M+Na].csv", sep = "/")
    write.csv(result, temp, row.names = F)
  }


  #-----------------------------------------------------------
  temp.idx <- grep(pattern = "NH4P", x = col.list)
  if (length(temp.idx)>0){
    temp <- which(raw.data[,1]=="exactmass")
    temp.data <- raw.data[-c(1:temp),]
    names(temp.data) <- raw.data[temp,]
    temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
    temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                            formula=temp.data$Formula)

    result <- lapply(seq(ncol(temp.data2)), function(i){
      temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
      temp <- data.frame(temp, rela.int="", stringsAsFactors = F)
      return(temp)
    })

    result <- do.call(cbind.data.frame, result)
    colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

    result <- data.frame(temp.data, result, stringsAsFactors = F)

    temp <- paste(file.address, "[M+NH4].csv", sep = "/")
    write.csv(result, temp, row.names = F)
  }

  #-------------------------------------------------------------
  temp.idx <- grep(pattern = "HN", x = col.list)
  if (length(temp.idx)>0){
    temp <- which(raw.data[,1]=="exactmass")
    temp.data <- raw.data[-c(1:temp),]
    names(temp.data) <- raw.data[temp,]
    temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
    temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                            formula=temp.data$Formula)

    result <- lapply(seq(ncol(temp.data2)), function(i){
      temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
      temp <- data.frame(temp, rela.int="", stringsAsFactors = F)
      return(temp)
    })

    result <- do.call(cbind.data.frame, result)
    colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

    result <- data.frame(temp.data, result, stringsAsFactors = F)

    temp <- paste(file.address, "[M-H].csv", sep = "/")
    write.csv(result, temp, row.names = F)
  }


  #-----------------------------------------------------------------
  temp.idx <- grep(pattern = "HCOON", x = col.list)
  if (length(temp.idx)>0){
    temp <- which(raw.data[,1]=="exactmass")
    temp.data <- raw.data[-c(1:temp),]
    names(temp.data) <- raw.data[temp,]
    temp.data2 <- temp.data[,temp.idx[1]:temp.idx[2]]
    temp.data <- data.frame(abbr.name=temp.data$Abbrev,
                            formula=temp.data$Formula)

    result <- lapply(seq(ncol(temp.data2)), function(i){
      temp <- round(as.numeric(unlist(temp.data2[,i])), digits = 4)
      temp <- data.frame(temp, rela.int="", stringsAsFactors = F)
      return(temp)
    })

    result <- do.call(cbind.data.frame, result)
    colnames(result)[seq(1, ncol(result), by = 2)] <- colnames(temp.data2)

    result <- data.frame(temp.data, result, stringsAsFactors = F)

    temp <- paste(file.address, "[M+HCOO].csv", sep = "/")
    write.csv(result, temp, row.names = F)
  }
}

