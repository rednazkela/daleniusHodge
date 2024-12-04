#' 
#' @title Stratification by Dalenius and Hodge 
#'
#' @description This function performs the stratification of a data vector 
#' using the Dalenius and Hodge method. It divides the data into a specified 
#' number of strata with lower and upper limits. 
#' 
#' @param obs A numeric vector containing the data to be stratified. 
#' @param numLevels An integer specifying the number of strata to calculate. 
#'
#' @return A matrix with three columns: 'level', 'lowerLimit', and 'upperLimit'. 
#' Each row of the matrix represents a stratum with its corresponding limits. 
#' 
#' @details 
#' The function sorts the data and divides them into homogeneous strata, 
#' using the cumulative square root of the frequencies to determine 
#' the stratum limits. 
#' 
#' @examples
#' \dontrun{ 
#' # Usage example 
#' data <- rnorm(100) 
#' levels <- 5 
#' strata <- daleniusHodge(data, levels) 
#' print(strata) 
#' } 
#' 
#' @export 
#' 
#' @date 2024-12-04 
#' @author 
#' Alexander Melendez 
#' @maintainer 
#' Alexander Melendez <rednazkela@gmail.com>
daleniusHodge <- function(obs, numLevels) {
  #Ordenar obs
  sortObs <- sort(obs, decreasing = FALSE)
  
  #número de intervalos
  numberIntervals <- floor(4.3*log10(length(sortObs)))
  
  #rangos homogeneos de los intervalos
  intervalRange = (sortObs[length(sortObs)] - sortObs[1]) / numberIntervals

  #Límites inferiores y superiores
  resultDt <- data.table(
    lowerLimit = numeric(), 
    upperLimit = numeric(), 
    indiscFrec = numeric(), 
    discFrec = numeric(), 
    sqrtFrec = numeric(), 
    sumSqrtFrec = numeric()
  )
  
  lowerLimit <- sortObs[1]
  upperLimit <- lowerLimit + intervalRange
  
  for (i in 1:numberIntervals) {
    #Limites y frecuencia indiscriminada
    dtTemp <- data.table(
      lowerLimit = lowerLimit,
      upperLimit = upperLimit,
      indiscFrec = sum(sortObs <= upperLimit, na.rm = TRUE),
      discFrec = 0, 
      sqrtFrec = 0, 
      sumSqrtFrec = 0
    )
    resultDt <- rbind(resultDt, dtTemp)
    
    #Frecuencia discriminada
    if(i == 1) {
      resultDt[.N, discFrec := indiscFrec]
    } else {
      resultDt[.N, discFrec := indiscFrec - resultDt[.N-1, indiscFrec]]
    }
    
    #Raiz cuadrada
    resultDt[.N, sqrtFrec := sqrt(discFrec)]
    
    #Raiz cuadrada acumulada
    if(i == 1) {
      resultDt[.N, sumSqrtFrec := sqrtFrec]
    } else {
      resultDt[.N, sumSqrtFrec := sqrtFrec + resultDt[.N-1, sumSqrtFrec]]
    }
    
    lowerLimit <- upperLimit + 0.00000001
    upperLimit <- lowerLimit + intervalRange
    
    #print(resultDt[.N])
  }
  
  maxSumSqrtFrec = resultDt[.N, sumSqrtFrec]
  frecLevels = maxSumSqrtFrec/numLevels
  #print(maxSumSqrtFrec)
  #print(frecLevels)
  column_names <- paste0("level", 1:numLevels)
  distanceDt <- data.table(matrix(NA, nrow = 0, ncol = numLevels))
  setnames(distanceDt, column_names)
  
  for(i in 1:numberIntervals) {
    dtTemp2 <- data.table()
    for (j in 1:numLevels) {
      temp_col_name <- paste0("level", j)
      temp_col_value <- abs((frecLevels * j) - resultDt[i, sumSqrtFrec])
      dtTemp2 <- cbind(dtTemp2, temp_col_value)
      setnames(dtTemp2, "temp_col_value", temp_col_name)
    }
    #print(dtTemp2)
    distanceDt <- rbind(distanceDt,dtTemp2)
  }
  
  resultLevels = data.table(
    level = character(), 
    lower = numeric(), 
    upper = numeric()
  )
  
  
  #print(distanceDt)
  for (col_name in column_names) {
    row_number <- which.min(distanceDt[[col_name]])
    #print(resultDt)
    resultLevels <- rbind(resultLevels, 
                          data.table(
                            level = col_name,
                            lower = resultDt[1, lowerLimit],
                            upper =resultDt[row_number, upperLimit]
                          ))
    resultDt <- resultDt[row_number + 1:.N]
    distanceDt <- distanceDt[row_number + 1:.N]
  }
  
  return(resultLevels)
}
