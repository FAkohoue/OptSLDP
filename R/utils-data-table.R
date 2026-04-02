# utils-data-table.R
# Declare data.table as a proper dependency so data.table 1.14+ cedta() check
# recognises OptSLDP as a data.table-aware package.  This file does nothing at
# runtime; its sole purpose is to generate the importFrom() lines in NAMESPACE
# that data.table requires when ':=' or other special operators are used.

#' @importFrom data.table ':='
#' @importFrom data.table '.N'
#' @importFrom data.table '.SD'
#' @importFrom data.table 'fread'
#' @importFrom data.table 'fwrite'
#' @importFrom data.table 'data.table'
#' @importFrom data.table 'as.data.table'
#' @importFrom data.table 'setDT'
#' @importFrom data.table 'setorder'
#' @importFrom data.table 'rbindlist'
#' @importFrom data.table 'copy'
#' @importFrom tools file_ext
#' @importFrom tools file_path_sans_ext
NULL
