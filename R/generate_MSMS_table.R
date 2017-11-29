#' @title generate_MSMS_table
#' @description  generate targeted MS/MS table for each compound
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param com.name Tompound name.
#' @param ext.mass The exact mass of compound.
#' @param adduct The adduct form to dissociation, including "[M+H]", "[M+Na]", "[M+NH4]","[M-H]", "[M+HCOO]"
#' @param ce The collision energy. Default: 20V
#' @param charge The charge form of ion. Default: 1
#' @param rt The rt range of targeted ion. Default: ""
#' @param delta.rt The delta rt range. Default: ""
#' @param iso.width The isolation window, including: "Wide", "Medium", "Narrow"
#'
#' @example


generate_MSMS_table <- function(com.name,
                                ext.mass,
                                adduct=c("[M+H]", "[M+Na]", "[M+NH4]",
                                         "[M-H]", "[M+HCOO]"),
                                ce=20,
                                charge=1,
                                rt="",
                                delta.rt="",
                                iso.width=c("Wide",
                                            "Medium",
                                            "Narrow")
                                ){
  if (length(ext.mass)==0) {
    stop("Please input extact mass.\n")
  } else {
    ext.mass <- as.numeric(ext.mass)
  }

  ext.mass <- sapply(seq(length(adduct)), function(i){
    temp <- switch(adduct[i],
                   "[M+H]"={ext.mass+1.0078},
                   "[M+Na]"={ext.mass+22.9898},
                   "[M+NH4]"={ext.mass+18.0344},
                   "[M-H]"={ext.mass-1.0078},
                   "[M+HCOO]"={ext.mass-44.9977})
  })

  iso.width <- switch (match.arg(iso.width),
                       "Wide"={"Wide (~9 m/z)"},
                       "Medium"={"Medium (~4 m/z)"},
                       "Narrow"={"Narrow (~1.3 m/z)"}
  )


  result <- data.frame(V1=TRUE,
                       V2=ext.mass,
                       V3=charge,
                       V4=rt,
                       V5=delta.rt,
                       V6=iso.width,
                       V7=ce,
                       V8="",
                       stringsAsFactors = F
                       )

  temp <- c('On', 'Prec. m/z', 'Z', 'Ret. Time (min)', 'Time (min) Delta', 'Iso. Width', 'Collision Energy', 'Acquisition Time (ms/spec)')
  result <- rbind(c("TargetedMSMSTable", rep("", 7)), temp, result)

  temp <- paste(com.name, "csv", sep = ".")
  write.csv(result, temp, row.names = F, col.names = F)

}