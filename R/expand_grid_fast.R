
    expand_grid_fast <- function(seq1,seq2) {
      # copied from http://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid
      cbind(Var1 = rep.int(seq1, length(seq2)), 
            Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2)))) 
    }
