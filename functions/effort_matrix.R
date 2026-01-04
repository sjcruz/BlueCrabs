



effort_matrix = function (E) {
  
  effort = rep(0, ii*months)
  effort[closed_seasons]=NA
  effort[!is.na(effort)] = E #only having the model optimize over months that are open to fishing to decrease optimization time 
  effort = effort/model_ratio #scale effort to model ratio 
  effort=matrix(effort, ncol =ii, byrow = F)
  effort[is.na(effort)] = 0
  
  
  return(effort)
}

################################################################################
# Finding mean, and 95% CI from multiple runs of each scenatio 
################################################################################


read_matrix_from_file <- function(file_list) {
  
  E = as.matrix(scan(file = paste0(results_dir, file_list) , what = numeric(), sep = " ") )
  effort = effort_matrix(E=E)
  
  return(effort)
}

summary_E = function(scenario_name=scenario_name){
  
  effort = data.frame(large_males=NA, females_pot=NA, females_dredge=NA, month=NA)
  sd = data.frame(large_males=NA, females_pot=NA, females_dredge=NA, month=NA)
  
  file_list = list.files(results_dir, pattern =scenario_name )
  matrix_list = lapply(file_list, read_matrix_from_file)
  
  n = length(matrix_list)
  
  if (n==1){
    mean = matrix_list[[1]]
    sd = matrix(0, nrow=nrow(mean), ncol=ncol(mean))
    
  } else {
    result = matrix(0, nrow = nrow(matrix_list[[1]]), ncol = ncol(matrix_list[[1]]))
    
    for (i in 1:n) {
      result = result + matrix_list[[i]]
    }
    
    mean = result/n
    
    sd = sqrt((1/n) * Reduce(`+`, lapply(matrix_list, function(x) (x - mean)^2)))
    
    
  }
  return(list(mean, sd ))
}


random_draws = function(mean_matrix, sd_matrix) {
  result = matrix(0, nrow = nrow(mean_matrix), ncol = ncol(mean_matrix))
  for (i in 1:nrow(mean_matrix)) {
    for (j in 1:ncol(mean_matrix)) {
      random_value = rnorm(1, mean = mean_matrix[i, j], sd = sd_matrix[i, j])
      result[i, j] = pmax(random_value, 0)
    }
  }
  
  
  return(result)
}


