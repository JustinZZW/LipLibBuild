#' @title Info_cal
#' @description  calculate formula, exact mass, and charge
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param raw.data The smiles structure as vector format.
#' @param is.output Whether output the csv file? The default is TRUE.
#' @return a matrix (n*3) including formula, exact mass, and charge
#' @example Info_Cal(raw.data=raw.data)

Info_Cal <- function(raw.data, is.output=TRUE){
  temp <- rcdk::parse.smiles(raw.data)

  cat("Start calculate the Formula and Exact mass\n")
  formula.info <- lapply(seq(length(temp)), function(i){
    if (i %in% seq(length(temp), length.out = 11)) {
      cat(round((i/length(temp))*100, digits = 2)); cat("%"); cat(" ")
    }

    temp.formula <- temp[[i]]
    temp.result <- rcdk::get.mol2formula(temp.formula)

    formula <- temp.result@string
    exact.mass <- round(as.numeric(temp.result@mass), digits = 4)
    charge <- temp.result@charge

    result <- c(formula, exact.mass, charge)
    names(result) <- c("formula", "exact.mass", "charge")

    return(result)
  })

  result <- do.call(rbind, formula.info)
  result <- data.frame(formula=result[,1], exact.mass=result[,2], charge=result[,3])

  if (is.output==TRUE) {
    cat("Start output result\n")
    temp <- paste("info", "csv", sep = ".")
    write.csv(result, temp, row.names = F)
  }

  return(result)

}
