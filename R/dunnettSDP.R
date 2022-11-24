#' Dunnett multiple comparison test - Step Down Procedure
#'
#' This function allows you to apply Dunnett multiple comparison test - Step Down Procedure.
#' @param d data.table with columns dose and endpoint
#' @keywords dunnett
#' @export
#' @examples
#' dunnettSDP(d)

dunnettSDP <- function(d){
  tmp <- d
  setDT(tmp)
  
  tmp <- tmp[!is.na(get('dose'))] # na.omit
  tmp <- tmp[!is.na(get('endpoint'))] # na.omit
  
  if(nrow(tmp) == 0){
    return(data.table("dose"=numeric(), "q-values"=numeric(), "p.signif"=numeric()))
    # return(tmp)
  }
  
  tmp[, c('N', 'mean'):= list(.N, mean(get('endpoint'))),  keyby = c('dose')]
  sample_statistics <-  tmp[, list('N'=.N, 'mean'=mean(get('endpoint'))),  keyby = c('dose')] # sort by dose
  tmp[, ('(endpoint - mean)^2'):= (get('endpoint') - get('mean'))**2] # auto clustered groups
  
  sssq <- sum(tmp[, get('(endpoint - mean)^2')])
  ns <- sample_statistics[, get('N')]
  n <- sum(ns) - length(ns)
  ssq <- sssq/n # A Multiple Comparison Procedure for Comparing Several Treatments with a Control - Dunnett 1955
  s <- sqrt(ssq) 
  
  x0 <- sample_statistics[1, get('mean')]
  n0 <- sample_statistics[1, get('N')]
  
  sample_statistics[, ('t_i'):=abs(get('mean')- x0)/(s*( sqrt( (1/n0)+(1/get('N')) ) ) ) ]
  N_groups <- nrow(sample_statistics) -1
  n_groups <- sample_statistics[2:nrow(sample_statistics), get('N')]
  
  corr.matrix <- matrix(0, N_groups, N_groups)
  diag(corr.matrix)=rep(1, N_groups)
  for(i in 1:N_groups -1){
    for(j in (i+1):N_groups){
      corr.matrix[i,j]=(1/n0)/(sqrt( (1/n_groups[i]) + (1/n0) )*sqrt( (1/n_groups[j]) + (1/n0) ) )
      corr.matrix[j,i]= corr.matrix[i,j]
    }
  }
  
  res <- pvAdjSDDT(sample_statistics[2:nrow(sample_statistics), get('t_i')], df=n, corr.matrix = corr.matrix)
  
  setDT(res)
  colnames(res)[1] <- 't_i'
  colnames(res)[2] <- 'q-values'
  
  res <- unique(res, by = 't_i')
  
  sample_statistics <- sample_statistics[res, on = .(t_i)]
  
  # mute q-values where N < 3
  sample_statistics[get('N') < 3, ('q-values'):=NA]
  
  sample_statistics <- sample_statistics[, c('dose', 'q-values')]
  sample_statistics[, ('p.signif'):=' ']
  
  # add significance codes
  for(i in 1:nrow(sample_statistics)){
    q_val <- sample_statistics[i, get('q-values')] # significance codes
    if(is.na(q_val)){
      sample_statistics[i, 'p.signif'] <- ''
    }
    else if(q_val <= 0.001){
      sample_statistics[i, 'p.signif'] <- '*'
    }
    else if(q_val <= 0.01){
      sample_statistics[i, 'p.signif'] <- '*'
    }
    else if(q_val <= 0.05){
      sample_statistics[i, 'p.signif'] <- '*'
    }
    else if(q_val <= 0.1){
      sample_statistics[i, 'p.signif'] <- '.'
    }
    else{
      sample_statistics[i, 'p.signif'] <- ''
    }
  }
  
  sample_statistics
}


pvAdjSDDT <- function(teststats, df, corr.matrix){
  k <- length(teststats)
  teststats.ordered <- sort(teststats)
  
  if(length(teststats.ordered) == 0){ # wenn all t_i NaN, i.e. s = 0
    return(data.frame('teststats'= teststats, 'qval' = NA))
  }
  
  qval <- rep(NA, k)
  qval[k] <- tryCatch({
    1 - mvtnorm::pmvt(lower = -abs(rep(teststats.ordered[k], k)), upper = abs(rep(teststats.ordered[k], k)), df = df, corr=corr.matrix)
  }, 
  error = function(cond){
    as.numeric(NA)
  })
  
  for(m in (k-1):1){
    qval[m] <- tryCatch({
      1 - mvtnorm::pmvt(lower = -abs(rep(teststats.ordered[m], m)), upper = abs(rep(teststats.ordered[m], m)), df = df, corr=corr.matrix[1:m, 1:m])
    },
    error = function(cond){
      as.numeric(NA)
    })
    qval[m] <- max(qval[m+1], qval[m], na.rm = TRUE)
  }
  
  res <- data.frame('teststats'= teststats.ordered, 'qval' = qval)
  res
}

