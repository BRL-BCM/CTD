#' Iterated Logarithm (Base 2)
#'
#' This function calculates the number of times the logarithm function (base 2) must be iteratively applied before the result is less than or equal to 1.
#' @param num - An integer.
#' @examples
#' iteratedLogarithm(4)
#' 2
stats.iteratedLog2 = function(num) {
  if (num<=1) {
    return(0)
  } else {
    iterator=1+stats.iteratedLog2(log2(num))
  }
  return (iterator)
}