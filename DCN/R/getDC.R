getDC <-
function(g){                          
  n <- ecount(g)
  if(n > 0){		
    DC <- sum(E(g)$DC) /(sqrt(n))	
  }else{
    DC <- 0
  }	
}
