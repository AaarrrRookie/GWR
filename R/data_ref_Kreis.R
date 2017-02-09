#' OSM - Postleitzahlen-Kreis - Referenzdaten
#'
#' Daten aus dem OSM Projekt
#'
#' @docType data
#'
#' @usage data(OSM_Kreis)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references OSM openStreetMap
#' (\href{http://www.geofabrik.de/}{Geofabrik})
#'
#' @source \href{https://www.suche-postleitzahl.org/downloads}{PLZ Archive}
#'
#' @examples
#' data(OSM_Kreis)
#' plz <- attr(OSM, "plz")
#' \donttest{hist(plz)}
"OSM_Kreis"
