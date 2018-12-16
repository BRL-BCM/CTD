#' Entropy of a bit-string
#'
#' The entropy of a bitstring (ex: 1010111000) is calculated.
#' @param x - A vector of 0's and 1's.
#' @export stats.entropyFunction
#' @examples
#' stats.entropyFunction(c(1,0,0,0,1,0,0,0,0,0,0,0,0))
#' > 0.6193822
#' stats.entropyFunction(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0))
#' > 1
#' stats.entropyFunction(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
#' > 0
stats.entropyFunction = function(bitString) {
  pT = sum(bitString)/length(bitString)
  pF = 1-pT
  if (pT==1 || pT==0) {
    e = 0
  } else {
    e = -pT*log2(pT)-pF*log2(pF)
  }
  return(e)
}
