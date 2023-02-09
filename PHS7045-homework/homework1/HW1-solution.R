N <- 228
arms <- 4
h0 <- 0.35
h1 <- 0.35
h2 <- 0.35
h3 <- 0.35

# Generate draws from posterior and calculate pmax
postDraws <- function(y,nMcmc,h0,h1,h2,h3,n0,n1,n2,n3){
  
  # Generate draws from posterior
  postDraws0 <- rbeta(n=nMcmc, shape1 = h0 + sum(y[1:n0,"0"]==1), shape2 = (1-h0) + n0 - sum(y[1:n0,"0"]==0))
  postDraws1 <- rbeta(n=nMcmc, shape1 = h1 + sum(y[1:n1,"1"]==1), shape2 = (1-h1) + n1 - sum(y[1:n1,"1"]==0))
  postDraws2 <- rbeta(n=nMcmc, shape1 = h2 + sum(y[1:n2,"2"]==1), shape2 = (1-h2) + n2 - sum(y[1:n2,"2"]==0))
  postDraws3 <- rbeta(n=nMcmc, shape1 = h3 + sum(y[1:n3,"3"]==1), shape2 = (1-h3) + n3 - sum(y[1:n3,"3"]==0))
  
  # Calculate allocation probabilities
  # v1-v3 probability each arm is best
  # v0 see RMatch Viele et al. 2020
  pBest <- pmax(postDraws0,postDraws1,postDraws2,postDraws3)
  
  v1 <- mean(pBest==postDraws1)
  v2 <- mean(pBest==postDraws2)
  v3 <- mean(pBest==postDraws3)
  
  v0 <- min(sum( c(v1,v2,v3) * (c( n1, n2, n3) + 1) / (n0 + 1), max(v1, v2, v3)) )
  
  # Standardize
  V0 <- v0 / (sum(v0,v1,v2,v3))
  V1 <- v1 / (sum(v0,v1,v2,v3))
  V2 <- v2 / (sum(v0,v1,v2,v3))
  V3 <- v3 / (sum(v0,v1,v2,v3))
  
  # Calculate probability each arm is greater than control
  p1 <- mean(postDraws1 > postDraws0)
  p2 <- mean(postDraws2 > postDraws0)
  p3 <- mean(postDraws3 > postDraws0)
  
  pMax <- max(p1,p2,p3)
  out <- c(V0=V0,V1=V1,V2=V2,V3=V3,p1=p1,p2=p2,p3=p3,pMax=pMax)
  return(out)
  
}



design1 <- function(y, nMcmc=10000, h0, h1, h2, h3){
  
  postDraws(y=y,nMcmc=nMcmc,
            n0=nrow(y)/ncol(y),
            n1=nrow(y)/ncol(y),
            n2=nrow(y)/ncol(y),
            n3=nrow(y)/ncol(y),
            h0=h0, h1=h1, h2=h2, h3=h3)
  
}

design2 <- function(N,looks, arms, nInterim,h0,h1,h2,h3){
  
  nt           <- matrix(NA,nrow=looks,ncol=arms)
  colnames(nt) <- 0:(arms-1)
  nt[1,]       <- rep(nInterim / arms,arms)
  size         <- c(rep(nInterim,interimLooks - 1),N - nInterim*interimLooks)
  
  for(i in seq(looks-1)){
    
    alloProbs <- postDraws(y=y, nMcmc = 1000, h0=h0,h1=h1,h2=h2,h3=h3,
                           n0 = nt[i,"0"],
                           n1 = nt[i,"1"],
                           n2 = nt[i,"2"],
                           n3 = nt[i,"3"])
    
    nt[i+1,] <- nt[i,] + c(rmultinom(n = 1, size = size[i], 
                                    prob = alloProbs[c("V0","V1","V2","V3")]))
    
  }

  post <- postDraws(y=y, nMcmc = 1000, h0=h0,h1=h1,h2=h2,h3=h3,
                    n0 = nt[i,"0"],
                    n1 = nt[i,"1"],
                    n2 = nt[i,"2"],
                    n3 = nt[i,"3"])
  
  return(post)

}



# # Generate outcomes
# y0 <- rbinom(N,size=1,prob=h0)
# y1 <- rbinom(N,size=1,prob=h1)
# y2 <- rbinom(N,size=1,prob=h2)
# y3 <- rbinom(N,size=1,prob=h3)
# 
# y <- cbind("0"=y0,"1"=y1,"2"=y2,"3"=y3)
# 
# 
# design1 <- function(y, nMcmc=10000, h0, h1, h2, h3){
#   
#   # Generate draws from posterior
#   postDraws0 <- rbeta(n=nMcmc, shape1 = h0 + sum(y[1:(N/arms),"0"]==1), shape2 = (1-h0) + N/arms - sum(y[1:(N/arms),"0"]==0))
#   postDraws1 <- rbeta(n=nMcmc, shape1 = h1 + sum(y[1:(N/arms),"1"]==1), shape2 = (1-h1) + N/arms - sum(y[1:(N/arms),"1"]==0))
#   postDraws2 <- rbeta(n=nMcmc, shape1 = h2 + sum(y[1:(N/arms),"2"]==1), shape2 = (1-h2) + N/arms - sum(y[1:(N/arms),"2"]==0))
#   postDraws3 <- rbeta(n=nMcmc, shape1 = h3 + sum(y[1:(N/arms),"3"]==1), shape2 = (1-h3) + N/arms - sum(y[1:(N/arms),"3"]==0))
#   
#   # Calculate probability each arm is greater than treatment
#   p1 <- mean(postDraws1 > postDraws0)
#   p2 <- mean(postDraws2 > postDraws0)
#   p3 <- mean(postDraws3 > postDraws0)
#   max(p1,p2,p3)
# }
# 
# design1(y=y,h0=h0,h1=h0,h2=h0,h3=h0)
# 
# 
# 
# design2 <- function(y, nMcmc=1000, h0, h1, h2, h3){
#   
#   interimLooks <- 5
#   nt <- matrix(NA,nrow=interimLooks,ncol=4)
#   colnames(nt) <- 0:3
#   nt[1,] <- rep(10,4)
#   
#   
#   for(i in seq(interimLooks-1)){
#     
#     # Generate draws from posterior
#     postDraws0 <- rbeta(n=nMcmc, shape1 = h0 + sum(y[1:nt[i,"0"],"0"]==1), shape2 = (1-h0) + N/arms - sum(y[1:nt[i,"0"],"0"]==0))
#     postDraws1 <- rbeta(n=nMcmc, shape1 = h1 + sum(y[1:nt[i,"1"],"1"]==1), shape2 = (1-h1) + N/arms - sum(y[1:nt[i,"1"],"1"]==0))
#     postDraws2 <- rbeta(n=nMcmc, shape1 = h2 + sum(y[1:nt[i,"2"],"2"]==1), shape2 = (1-h2) + N/arms - sum(y[1:nt[i,"2"],"2"]==0))
#     postDraws3 <- rbeta(n=nMcmc, shape1 = h3 + sum(y[1:nt[i,"3"],"3"]==1), shape2 = (1-h3) + N/arms - sum(y[1:nt[i,"3"],"3"]==0))
#     
#     # Calculate probability each arm is greater than treatment
#     pMax <- pmax(postDraws0,postDraws1,postDraws2,postDraws3)
#     
#     v1 <- mean(pMax==postDraws1)
#     v2 <- mean(pMax==postDraws2)
#     v3 <- mean(pMax==postDraws3)
#     
#     v0 <- min(sum( c(v1,v2,v3) * 
#                      (nt[i,-which(colnames(nt)=="0")] + 1) / (nt[i,which(colnames(nt)=="0")] + 1)), 
#               max(v1, v2, v3))
#     
#     # Standardize
#     V0 <- v0 / (sum(v0,v1,v2,v3))
#     V1 <- v1 / (sum(v0,v1,v2,v3))
#     V2 <- v2 / (sum(v0,v1,v2,v3))
#     V3 <- v3 / (sum(v0,v1,v2,v3))
#     
#     nt[i+1,] <- nt[i,] + c(rmultinom(n = 1, size = 40, prob = c(V0,V1,V2,V3)))
#     
#   }
#   
#   # Generate draws from posterior
#   postDraws0 <- rbeta(n=nMcmc, shape1 = h0 + sum(y0[1:nt[i,"0"]]==1), shape2 = (1-h0) + N/arms - sum(y0[1:nt[i,"0"]]==0))
#   postDraws1 <- rbeta(n=nMcmc, shape1 = h1 + sum(y1[1:nt[i,"1"]]==1), shape2 = (1-h1) + N/arms - sum(y1[1:nt[i,"1"]]==0))
#   postDraws2 <- rbeta(n=nMcmc, shape1 = h2 + sum(y2[1:nt[i,"2"]]==1), shape2 = (1-h2) + N/arms - sum(y2[1:nt[i,"2"]]==0))
#   postDraws3 <- rbeta(n=nMcmc, shape1 = h3 + sum(y3[1:nt[i,"3"]]==1), shape2 = (1-h3) + N/arms - sum(y3[1:nt[i,"3"]]==0))
#   
#   # Calculate probability each arm is greater than treatment
#   pMax <- pmax(postDraws0,postDraws1,postDraws2,postDraws3)
#   
#   p1 <- mean(pMax==postDraws1)
#   p2 <- mean(pMax==postDraws2)
#   p3 <- mean(pMax==postDraws3)
#   
#   max(p1,p2,p3)
# }
# 
# 
# design2(y=y,h0=h0,h1=h0,h2=h0,h3=h0)
