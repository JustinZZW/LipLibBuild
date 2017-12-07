#' @title average_int
#' @author Zhiwei Zhou
#' @description Averge intensity of same fragment in multiple file.
#' @param dir.path the files direction path
#' @return avg.data list; avg.result: average result for each fragment;
#'     int.result: intensity of all files

average_int <- function(dir.path="."){
  raw.data <- lapply(seq(length(dir(dir.path))), function(i){
    temp <- readr::read_csv(dir()[i])
  })

  all.annotation <- lapply(raw.data, function(x){
    x$annotation
  })

  all.annotation <- do.call(c, all.annotation)
  all.annotation <- unique(all.annotation)

  int.result <- lapply(seq(length(all.annotation)), function(i){
    temp <- all.annotation[i]
    int <- sapply(seq(length(raw.data)), function(i){
      search.data <- raw.data[[i]]$annotation
      idx <- which(search.data==temp)
      if (length(idx) > 0){
        result <- raw.data[[i]]$int[idx]
      } else {
        result <- NA
      }
      return(result)
    })
    return(int)
  })

  names(int.result) <- all.annotation
  avg.result <- lapply(int.result, function(x){
    round(mean(x, na.rm = T), digits = 0)
  })
  avg.result <- do.call(c, avg.result)
  avg.result <- data.frame(annotation=all.annotation,
                           int=avg.result,
                           stringsAsFactors = F)

  write.csv(avg.result, "intensity template.csv", row.names = F)

  avg.data <- list(avg.result=avg.result, int.result=int.result)
  save(avg.data, file = "avg.data")
}
