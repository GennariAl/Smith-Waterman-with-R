## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(stringr)
library(pairwiseLocalAlignment)

## -----------------------------------------------------------------------------
seq1 <- 'CCTATACCTAATTTTCGGCGCATGAGCCGGAATGGTGGGTACCGCTCTAAGCCTCCTCATTC'
seq2 <- 'CATGAGCTGGAATAGTAGGTACCGCCCTAAGCCTCCTAATTCGAGCAGAGCTAGGCCAACCC'

smith_waterman(10, -10, -3, seq1, seq2)

