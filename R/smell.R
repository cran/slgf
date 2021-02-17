#' Smell data on olfactory function by age group
#'
#' O'Brien and Heft (1995) studied the University of Pennsylvania Smell
#' Identification Test (UPSIT). 180 subjects of different age groups were
#' asked to describe 40 different odors. Olfactory index was quantified by
#' the Freeman-Tukey modified arcsine transformation on the proportion of
#' correctly identified odors. Subjects were divided into five age groups:
#' group 1 if age 2 or younger; group 2 if between ages 26 and 40; group 3
#' if between ages 41 and 55; group 4 if between ages 56 and 70; and group
#' 5 if older than 75.

#'
#' @docType data
#'
#' @usage data(smell)
#'
#' @format A data frame with 180 rows and 2 variables:
#' \describe{
#'   \item{agecat}{age category, from 1 to 5.}
#'   \item{olf}{olfactory function, measured as the Freeman-Tukey modified arcsine transformation on the proportion of
#' correctly identified odors.}
#' }
#'
#' @references
#' \insertRef{ObrienHeft}{slgf}
#'
#' @source \href{https://documentation.sas.com/?docsetId=statug&docsetTarget=statug_glm_examples10.htm&docsetVersion=15.2&locale=zh-TW}{SAS/STAT 15.2 User's Guide}
#'
#' @keywords datasets
#'
#'
"smell"
