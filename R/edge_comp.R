#' Generate a microtable object with paired nodes distributions of edges across networks
#'
#' @description
#' Generate a microtable object with paired nodes distributions of edges across networks. Useful for the edge comparisons across different networks.
#' The return otu_table in microtable object has the binary numbers in which 1 represents the presence of the edge in the corresponding network.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of microeco package.
#' @return microtable object
#' @examples
#' data(soil_amp_network)
#' test <- edge_comp(soil_amp_network)
#' # test is a microtable object
#' 
#' @export
edge_comp <- function(network_list){
	check_input(network_list)
	venn_table <- NULL
	for(i in names(network_list)){
		venn_table <- get_edge_pair(network_list[[i]], venn_table, i)
	}
	rownames(venn_table) <- venn_table[, 1]
	venn_table <- venn_table[, -1]
	venn_table[is.na(venn_table)] <- 0
	# create a new microtable object
	tmp <- suppressMessages(microtable$new(otu_table = venn_table))
	tmp
}

# inner function
get_edge_pair <- function(network, raw_table, network_name){
	if(is.null(network$res_edge_table)){
		suppressMessages(network$get_edge_table())
	}
	edge_nodes <- network$res_edge_table[, 1:2] %>% t %>% as.data.frame
	# make the names of paired nodes ordered
	sorted_edge_nodes <- lapply(edge_nodes, function(x){sort(x) %>% paste0(., collapse = " -- ")}) %>% 
		do.call(rbind, .) %>% as.data.frame
	sorted_edge_nodes[, 2] <- 1
	# remove duplicates that may be in directed network
	sorted_edge_nodes %<>% .[!duplicated(.), ]
	colnames(sorted_edge_nodes) <- c("name", network_name)
	if(is.null(raw_table)){
		sorted_edge_nodes
	}else{
		dplyr::full_join(raw_table, sorted_edge_nodes, by = c("name" = "name"))
	}
}



