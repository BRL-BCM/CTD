#' Entropy of a bit-string
#'
#' The entropy of a bitstring (ex: 1010111000) is calculated.
#' @param bitString - A vector of 0's and 1's.
#' @return e - a floating point percentage, between 0 and 1.
#' @export stat.entropyFunction
#' @examples
#' stat.entropyFunction(c(1,0,0,0,1,0,0,0,0,0,0,0,0))       # Output: 0.6193822
#' stat.entropyFunction(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)) # Output: 1
#' stat.entropyFunction(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)) # Output: 0
stat.entropyFunction = function(bitString) {
    pT = sum(bitString)/length(bitString)
    pF = 1-pT
    if (pT==1 || pT==0) {
        e = 0
    } else {
        e = -pT*log2(pT)-pF*log2(pF)
    }
    return(e)
}
