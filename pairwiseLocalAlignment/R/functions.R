#' Pairwise Local Alignment
#'
#' This function reads 3 numerical parameters which defines the cost of a 
#' match, mismatch or the insertion of a gap between nucleotides 
#' of two sequences aligned. Then it perform the Smith-Waterman algorithm for
#' Pairwise Local Alignment of sequences and return an list collecting the two 
#' aligned sub-sequences.
#'
#' @usage smith_waterman(match_cost, mismatch_cost, gap_cost, seq1, seq2)
#' @param match_cost score of a match between two nucleotides
#' @param mismatch_cost score of a mismatch between two nucleotides
#' @param gap_cost score of the insertion gap between two nucleotides
#' @param seq1 first nucletide sequence
#' @param seq2 second nucletide sequence
#' @importFrom stringr "str_split_fixed"
#' @import testthat
#' @return list with two elements: s1 and s2 which correspond to the two 
#' sub-sequences resulting from the pairwise local alignment
#' @author Alessandro Gennari\cr Politecnico di Milano\cr Maintainer: Alessandro 
#' Gennari\cr E-Mail: <alessandro.gennari@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith–Waterman_algorithm}\cr
#' @examples
#'
#' smith_waterman(10, -5, -6, "GAATC", "CATACG")
#'
#' @export
smith_waterman <- function(match_cost, mismatch_cost, gap_cost, seq1, seq2) {
  if (!is.numeric(match_cost)) {
    stop("match_cost argument is not numeric!")
  }
  if (!is.numeric(mismatch_cost)) {
    stop("mismatch_cost argument is not numeric!")
  }
  if (!is.numeric(gap_cost)) {
    stop("gap_cost argument is not numeric!")
  }
  if (!is.character(seq1)) {
    stop("seq1 argument is not character!")
  }
  if (!is.character(seq2)) {
    stop("seq2 argument is not character!")
  }
  if (match_cost <= 0) {
    warning("match cost is <= that 0!")
  }
  if (match_cost <= mismatch_cost) {
    warning("WARNING: match cost is <= that mismatch cost!")
  }
  rownucleotides <- as.vector(str_split_fixed(seq1, pattern = "", n = nchar(seq1)))
  colnucleotides <- as.vector(str_split_fixed(seq2, pattern = "", n = nchar(seq2)))
  scores <- matrix(data = NA, nrow = nchar(seq1)+1, ncol = nchar(seq2)+1)
  rownames(scores) <- c('', rownucleotides)
  colnames(scores) <- c('', colnucleotides)
  
  directions <- scores    # matrix collecting directions of trace back
  
  # initialization
  scores[1,] <- 0
  scores[,1] <- 0
  traceback_start <- c(0, 0, 0)
  names(traceback_start) <- c('score', 'row', 'column')
  
  # scores computation
  for (i in seq(from = 2, to = nrow(scores))) {
    for (j in seq(from = 2, to = ncol(scores))) {
      diag <- 0
      up <- 0
      left <- 0
      up <- scores[i-1, j] + gap_cost
      left <- scores[i, j-1] + gap_cost
      if (rownames(scores)[i] == colnames(scores)[j]) {
        diag <- scores[i-1, j-1] + match_cost      # score in case of match
      } else {
        diag <- scores[i-1, j-1] + mismatch_cost  # score in case of mismatch
      }
      
      # choose the max score
      scores[i, j] <- max(diag, up, left)
      directions[i, j] <- which(c(diag, up, left) == scores[i, j])[1]
      # in case of equal values, diag is chosen because I put [1]
      # 1 means diag, 2 means up, 3 means left
      
      if (scores[i, j] < 0) {
        scores[i, j] <- 0
      }
      # save the actual max score for traceback
      if (scores[i, j] >= traceback_start[1]) {
        traceback_start[1] <- scores[i, j]
        traceback_start[2] <- i
        traceback_start[3] <- j
      }
    }
  }
  
  # traceback
  i <- traceback_start[[2]]
  j <- traceback_start[[3]]
  seq_ali_1 <- ''
  seq_ali_2 <- ''
  
  while (scores[i, j] != 0) {
    if (directions[i, j] == 1) {  # diagonal
      seq_ali_1 <- paste(rownames(scores)[i], seq_ali_1, sep = '')
      seq_ali_2 <- paste(colnames(scores)[j], seq_ali_2, sep = '')
      j <- j - 1
      i <- i - 1
    } else if (directions[i, j] == 2) {  # up
      seq_ali_1 <- paste(rownames(scores)[i], seq_ali_1, sep = '')
      seq_ali_2 <- paste('-', seq_ali_2, sep = '')
      i <- i - 1
    } else if (directions[i, j] == 3) {  # left
      seq_ali_1 <- paste('-', seq_ali_1, sep = '')
      seq_ali_2 <- paste(colnames(scores)[j], seq_ali_2, sep = '')
      j <- j - 1
    }
  }
  
  return(list(seq_ali_1, seq_ali_2))
}