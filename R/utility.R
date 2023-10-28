

measure_eff <- function(all_networks){
	lapply(all_networks, function(y){
		lapply(y, function(x){
			if(length(igraph::V(x)) < 2 | length(igraph::E(x)) < 2){
				0
			}else{
				dis_matrix <- igraph::distances(x)
				nodes_num <- ncol(dis_matrix)
				dis_num <- stats::as.dist(dis_matrix)
				res_single <- sum(1/dis_num)/(nodes_num * (nodes_num - 1))
				res_single
			}
		})
	})
}

# Check whether the input list has the correct format
check_input <- function(input){
	if(!is.list(input)){
		stop("The input data must a list! Please check the data input!")
	}
	for(i in names(input)){
		if(!inherits(input[[i]], "trans_network")){
			stop(i, " in the input list is not a trans_network object! Please check it!")
		}
	}
}

