#' OSM - Postleitzahlen-Koordination - Referenzdaten
#'
#' Daten aus dem OSM Projekt
#'
#' @docType data
#'
#' @usage data(Hierarchie_Daten)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references Hierarchie_Daten
#'
#' @source \href{https://www.suche-postleitzahl.org/downloads}{PLZ Archive}
#'
#' @examples
#' data(OSM)
#' plz <- attr(OSM, "plz")
#' \donttest{hist(plz)}
"Hierarchie_Daten"
