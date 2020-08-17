#' @title Untargeted metabolomic analysis for the clinical screening of inborn errors
#' of metabolism
#'
#' @description Global metabolic profiling obtained by untargeted mass spectrometry-based
#' metabolomic platform for the detection of novel and known inborn errors of
#' metabolism.  This untargeted approach collected z-score values for ~900 unique
#' compounds (including ~500 named human analytes) from human plasma.  Data set
#' contains 190 individual plasma samples (120 confirmed inborn errors of
#' metabolism).  The outcome describes excellent sensitivity and specificity for
#' the detection of a wide rage of metabolic disorders and identified novel
#' biomarkers for some diseases.
#'
#' @name Miller2015
#' @aliases Miller2015
#' @docType data
#' @usage data(Miller2015)
#' @format Miller2015 - The data matrix (metabolite features are rows, patient observations are columns) for 186 untargeted metabolomics 
#'                       patient samples, alongside metabolite annotations.
#' @format diagnoses - A data.frame where all patient IDs are mapped to their given biochemical diagnosis.
#' @keywords datasets
#' @references Miller et al. (2015) J Inherit Metab Dis. 2015; 38: 1029–1039
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4626538/}{PubMed})
#' @source \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4626538/bin/10545_2015_9843_MOESM1_ESM.xls}{Dataset}
#'
#' @examples
#' require(CTD)
#' data(Miller2015)
NULL
