#' @title Calculate evapo-transpiration stats from monthly weather data
#' @description
#' \code{evapoTranspirationStats} Simple function to clean up ET stats generation
#' from weather data.
#'
#' @param weather A data.frame with the following columns (names must match):
#' "site", "month", "year", "tmin", "tmax", "tmean", "rain".
#' @param siteInfo A dataframe with the column "site" that matches that of weather. This
#' dataset is populated with ET stats
#' @param sites A character vector specifying the names of sites in weathe / siteInfo
#' @param ... Not currently in use.
#' @details See SPEI R package for details. This is just a simple wrapper.
#' @return siteInfo, with three additional columns, spi, spei, pot.et
#'
#' @examples
#' \dontrun{
#' # more here soon.
#' }
#'
#' @import SPEI
#' @export
evapoTranspirationStats<-function(weather, siteInfo, sites){
  do.call(rbind, lapply(sites, function(x){
    tmp.weather<-weather[weather$site == x,]
    tmp.site<-siteInfo[siteInfo$site == x,]
    out<-tmp.weather[with(tmp.weather, order(year, month)),]
    pot.et<-SPEI::thornthwaite(Tave = out$tmean,
                               lat = tmp.site$Latitude, na.rm=T)
    out$pot.et<-as.numeric(pot.et)
    out$spei <- as.numeric(with(out, spei(rain-pot.et,1, na.rm=T))$fitted)
    out$spi <- as.numeric(with(out, spi(rain,1, na.rm=T))$fitted)
    return(out)
  }))
}
