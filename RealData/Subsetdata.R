

proportion.keep <- 0.8

n.ind <- dim(Zuni12.pema.Ch.secondary[[1]])[1]
sample(rownames(Zuni12.pema.Ch.secondary[[1]]))
set.seed(567); keep <- sample(1:n.ind,round(n.ind*proportion.keep))
keep <- sort(keep)



Zuni12.secondary.keep <- list()
Zuni12.secondary.test <- list()
for(i in 1:length(Zuni12.pema.Ch.secondary)){
  Zuni12.secondary.keep[[i]] <- Zuni12.pema.Ch.secondary[[i]][keep,] 
  Zuni12.secondary.keep[[i]] <- Zuni12.pema.Ch.secondary[[i]][-keep,] 
}