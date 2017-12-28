#' @title Process and plot output from getWeatherStationData
#' @description
#' \code{procWeather} Take output from getWeatherStationData, get sliding window average
#' temperatures and cumulative precipiation. Use ggplot2 to plot these.
#'
#' @param histWeather Output from getWeatherStationData
#' @param years2plot Numeric, The years that will be plotted as individual lines. The other years in the
#' dataset will be averaged into "historical" data.
#' @param siteLevels Character, The levels to use for the names of the sites. If specified, ensure that
#' all sites are included - these will be the only ones that will be plotted.
#' @param plotit Logical, should a plot be drawn.
#' @param date_max The latest data that must have data from the station.
#' See rnoaa::meteo_nearby_stations.
#' @param limit The maximum number of stations to search. If none of these stations have sufficient
#' data, an error is returned. See rnoaa::meteo_nearby_stations.
#' @param verbose Logical, should status updates be printed?
#' @param ... Not currently in use.
#' @details This function iterates through NOAA weather collection sites, deterimining
#' whether the data collected satisfies the time and datatype constraints specified.
#' @return A data.frame of daily weather data for each location id. If more than one id
#' are included, the data.frames are placed in a named list.
#'
#' @examples
#' \dontrun{
#' ll = data.frame(latitude = c(27.54986,30.38398),
#'   longitude = c(-97.88101,-97.72938),
#'   id = c("KING","PKLE"),
#'   stringsAsFactors = F)
#' data(station_data)
#' histWeather<-getWeatherStationData(lat_lon_df = ll,
#'   station_data = station_data,
#'   date_min = "2000-01-01",
#'   hasCols = c("prcp","tmax","tmin"),
#'   date_max = "2016-12-31")
#' lapply(histWeather, head)
#' out = procWeather(histWeather=histWeather, years2plot = 2012:2016)
#' }

procWeather = function(histWeather, years2plot = c(2016, 2017),
                       siteLevels=NULL, plotit = T){

  histWeather<-do.call(rbind, histWeather)
  histWeather$tmean = with(histWeather, (tmin+tmax)/2)
  histWeather$rain[is.na(histWeather$rain)]<-0

  histAvg = histWeather
  histAvg$year = ifelse(histAvg$year %in% years2plot,
                        as.character(histAvg$year), "historical")
  histAvg<-ddply(histAvg, .(site, year, jd), summarize,
                 rain = mean(rain, na.rm = T),
                 tmean = mean(tmean, na.rm = T))
  tp = ddply(histAvg, .(site, year), mutate,
             cumPrecip = cumsum(rain),
             sw.tmean = rollapply(tmean, width = 20,
                                  partial = F, fill = NA,
                                  FUN = function(x) mean(x, na.rm = T)))

  tpl = reshape2::melt(tp, measure.vars = c("sw.tmean","cumPrecip"))
  tpl$year = factor(tpl$year, levels = c("historical", as.character(years2plot)))
  if(!is.null(siteLevels)){
    tpl$site = factor(tpl$site, levels = siteLevels)
  }
  colfunc <- colorRampPalette(c("darkred","yellow","cyan","darkblue"), space = "Lab")
  cols = c(rgb(0,0,0,.1), colfunc(length(years2plot)))
  lwds = c(3, rep(.5, length(years2plot)))

  tp1 = tpl[tpl$variable == "sw.tmean",]
  tp2 = tpl[tpl$variable == "cumPrecip",]

  p1 = ggplot(tp1, aes(x = jd, y = value, col = year, group = year, size = year))+
    geom_line()+
    scale_color_manual(values = cols, guide = F)+
    scale_size_manual(values = lwds, guide = F)+
    facet_wrap(~site, scale = "free", ncol = 1)+theme_jtl()+
    labs(x = "Julian Days", y = "20-day mean Temp (deg. C)")
  p2 = ggplot(tp2, aes(x = jd, y = value, col = year, group = year, size = year))+
    geom_line()+
    scale_color_manual(values = cols, guide = F)+
    scale_size_manual(values = lwds, guide = F)+
    facet_wrap(~site, scale = "free", ncol = 1)+theme_jtl()+
    labs(x = "Julian Days", y = "Cumulative Precipitation (mm)")

  if(plotit){
    grid.arrange(p1,p2, ncol = 2)
  }
  return(list(histWeather,histAvg,list(p1,p2)))
}



