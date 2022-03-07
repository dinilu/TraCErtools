#' Load a cropped version of downscaled TraCE21ka .nc files
#'
#' @param file A string with full name (path and file name) of the .nc file to be loaded from the hard disk drive.
#' @param sf A sf object from the sf package to be used as template or mask to crop the downscaled trace .nc file.
#'
#' @return A stars object with the downscaled trace data.
#' @export
#'
#' @examples # TBW
.read_dstrace_subset <- function(file, sf = NULL) {
  data <- stars::read_stars(file)
  if(!is.null(sf)){
    data <- sf::st_crop(data, sf)
  }
  data
}

#' Load downscaled TraCE21ka files by years.
#'
#' @param folder A character string with the path to the folder where the downscaled TraCE21ka files, ordered in folders by variables.
#' @param var A character string with the name of the variables to be loaded.
#' @param y_start A number with the first year to be loaded.
#' @param y_end A number with the last year to be loaded.
#' @param sf A sf object to be used for cropping the object.
#'
#' @return A stars object with the downscaled TraCE21ka data, cropped using the sf object if provided.
#'
#' @details Years should be in the format calibrated Before Present (e.g. 0 for loading year 1950 in the gregorian calendar).
#' @export
#'
#' @examples # TBW
read_dstrace_subset <- function(folder, var, y_start, y_end, sf = NULL){

  if(y_start >= y_end){
    stop("Starting year same or greater than ending year. Provide a y_start value lower than y_end.")
  }

  y_start <- as.numeric(y_start)
  y_end <- as.numeric(y_end)

  if(y_start >= y_end){
    stop("Starting year same or higher than ending year. Provide a y_start value lower than y_end.")
  }

  if(y_start < 0 & y_end > 0){
    y_seq <- c(y_start:-1, 1:y_end)
  } else {
    y_seq <- c(y_start:y_end)
  }

  files <- paste0(folder, "/", var, "/", var, y_seq, ".nc")

  if(length(files) > 1){
    data <- lapply(files, FUN=.read_dstrace_subset, sf)
    data <- do.call(c, data)
  } else {
    data <- .read_dstrace_subset(files, sf)
  }

  names(data) <- var
  
  time_2_calendar_dates(data, y_start, y_end)
}



#' Create a calendar dates by specified time intervals.
#'
#' @param y_start A number with the first year for the required dates.
#' @param y_end A number with the last year for the required dates.
#' @param by A number or string (as in lubridate package) to be used as interval to get the calendar dates. Default to "1 month" which is returning monthly dates (first day of the month).
#'
#' @return A vector with dates from the y_start to the y_end years at the time intervals specified by "by" argument.
#' @export
#'
#' @examples # TBW
calendar_dates <- function(y_start, y_end, by = "1 month"){
  cal_start <- lubridate::ymd("1950-01-01") + lubridate::years(y_start)
  cal_end <- lubridate::ymd("1950-12-01") + lubridate::years(y_end)
  cal_dates <- seq(cal_start, cal_end, by = by)
  cal_dates <- cal_dates[lubridate::year(cal_dates) != 0]
  cal_dates
}



#' Modify time to calendar years
#'
#' Modify a stars object to change the time dimension to calendar years between specified starting and ending years.
#'
#' @param data A stars object with the data to be modified.
#' @param y_start A number with the starting year of the data (in calibrated Before Present format).
#' @param y_end A number with the ending year of the data (in calibrated Before Present format).
#' @param by A number or string specifying the interval of the dates to be used in the new stars object.
#'
#' @return A stars object as in data argument but with changed time dimension.
#' @export
#'
#' @examples # TBW
time_2_calendar_dates <- function(data, y_start, y_end, by = "1 month") {
  cal_dates <- calendar_dates(y_start, y_end, by = by)
  stars::st_dimensions(data)$time$values <- cal_dates
  data
}




#' Read original TraCE21ka datasets 
#'
#' @param folder TBW
#' @param var TBW
#' @param point TBW
#'
#' @return TBW
#' @export
#'
#' @examples # TBW 
read_trace_subset <- function(folder, var, point){
  # folder <- "../Data/TraCE21ka"
  # var <- "TS"
  files <- list.files(paste0(folder, "/", var), full.names = TRUE, pattern=".nc")
  data <- lapply(files, FUN = stars::read_ncdf)
  data <- do.call(c, data)
  
  # gc()
  
  data <- time_2_calendar_dates(data, -22000, 40)
  
  if(var %in% c("TS", "TSMX", "TSMN")){
    data <- kelvin2celsius(data)
  }  
  if(var == "PRECC"){
    data <- flux2mm(data)
  }  
  
  lat <- point[2]
  lats <- stars::st_dimensions(data)$lat$value
  # dif.lats <- lats - lat
  # i <- which.max(dif.lats[dif.lats < 0])
  i <- findInterval(lat, lats)
  
  data <- dplyr::filter(data, lat > floor(lats[i]), lat < ceiling(lats[i+1]))
  # data <- filter(data, lon > -2, lon < 5, lat > 37, lat < 46)
  sf::st_crop(data, point)
} 


#' Transform stars objects from degrees kelvin to degrees celsius
#'
#' @param data A stars object with temperature as degrees kelvin.
#'
#' @return TBW
#' @export
#'
#' @examples # TBW  
kelvin2celsius <- function(data) {
  var <- names(data)
  # data_backup <- data
  data_class <- class(data[[var]])  
  data[[var]]  <- as.numeric(data[[var]]) - 273.15
  # class(data[[var]]) <- class(data_backup[[var]])
  class(data[[var]]) <- data_class
  units <- list(numerator = "\u00B0C", denominator = character(0))
  class(units) <- "symbolic_units"
  attr(data[[var]], "units") <- units 
  data
}


#' Title
#'
#' @param data TBW
#'
#' @return TBW
#' @export
#'
#' @examples # TBW
flux2mm <- function(data){
  var <- names(data)
  data_class <- class(data[[var]]) 
  data[[var]]  <- as.numeric(data[[var]]) * 2592000000
  class(data[[var]]) <-data_class
  units <- list(numerator = "mm", denominator = character(0))
  class(units) <- "symbolic_units"
  attr(data[[var]], "units") <- units 
  data
} 






#' Read downscaled CMIP5 datasets
#'
#' @param folder TBW
#' @param var TBW
#' @param scen TBW
#' @param gcm TBW
#' @param sf TBW
#'
#' @return TBW
#' @export 
#'
#' @examples #TBW  
read_dscmip_subset <- function(folder, var, scen, gcm, sf = NULL){

  files <- paste0(folder, "/", scen, "/", gcm, "/", var, "/", var, "1991-2100.nc")
  
  data <- stars::read_stars(files)
  if(!is.null(sf)){
    data <- sf::st_crop(data, sf)
  } 

  names(data) <- var
  
  data
}


