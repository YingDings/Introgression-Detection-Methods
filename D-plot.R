Dplot<- function(tree,Doutcome,taxa){
  library(ape)
  tre <- read.tree(text=tree)
  from_one <- to_one <- from_two<- to_two <- vector()
  o <- p <- 1
  outcome_ABBA <- Doutcome
  tip_ABBA <- taxa
  if(length(outcome_ABBA$Z[which(outcome_ABBA$Z > 3)])>0){
    for (j in 1:length(outcome_ABBA$Z[which(outcome_ABBA$Z > 3)])) {
      if(outcome_ABBA$ABBA[which(outcome_ABBA$Z>3)][j]>outcome_ABBA$BABA[which(outcome_ABBA$Z>3)][j]){
        from_one[o] <- as.character(tip_ABBA$p2[j])
        to_one[o] <- as.character(tip_ABBA$p3[j])
        o <- o+1
      }else{
        from_two[p] <- as.character(tip_ABBA$p1[j])
        to_two[p] <- as.character(tip_ABBA$p3[j])
        p <- p+1
      }
    }
    from <- c(from_one,from_two)
    for (k in 1:length(from)) {
      from[k] <- strsplit(as.character(from),"'")[[k]][2]
      from
    }
    to <- c(to_one,to_two)
    for (l in 1:length(to)) {
      to[l] <- strsplit(as.character(to),"'")[[l]][2]
      to
    }
    from_end <- vector()
    for (q in 1:length(from)) {
      from_end[q] <- grep(from[q],tre$tip.label)
    }
    to_end <- vector()
    for (w in 1:length(to)) {
      to_end[w] <- grep(to[w],tre$tip.label)
    }
    result <- evonet(tre,from_end,to_end)
    return(result)
  }
}
