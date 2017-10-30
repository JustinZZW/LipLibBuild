#' @title Info_cal
#' @description  calculate formula, exact mass, and charge
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param raw.data The smiles structure as vector format.
#' @return a matrix (n*3) including formula, exact mass, and charge
#' @example Info_Cal(raw.data=raw.data)

Info_Cal <- function(raw.data){
  temp <- raw.data$smiles[1:5]
  temp <- rcdk::parse.smiles(temp)

  formula.info <- lapply(seq(length(temp)), function(i){
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
  return(result)
}
