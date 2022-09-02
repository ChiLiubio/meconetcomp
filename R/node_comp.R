#' Generate a microtable object with node distributions across networks
#'
#' @description
#' Generate a microtable object with node distributions across networks. Useful for the node information comparisons across different networks.
#' The return otu_table in microtable object has the binary number that 1 represent presence of the node in the corresponding network.
#' 
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of microeco package.
#' @param degree default 1; the threshold of degree of nodes remained in the network.
#' @return microtable object
#' @examples
#' \dontrun{
#' data(soil_amp_network)
#' test <- node_comp(soil_amp_network)
#' # test is a microtable object
#' }
#' @export
node_comp <- function(network_list, degree = 1){
	check_input(network_list)
	for(i in names(network_list)){
		if(is.null(network_list[[i]]$res_node_table)){
			network_list[[i]]$get_node_table(node_roles = FALSE)
		}
		tmp <- network_list[[i]]$res_node_table %>% .[.$degree >= degree, "name", drop = FALSE]
		tmp[, i] <- 1
		if(i == names(network_list)[1]){
			venn_table <- tmp
		}else{
			venn_table <- dplyr::left_join(venn_table, tmp, by = c("name" = "name"))
		}
	}
	rownames(venn_table) <- venn_table[, 1]
	venn_table <- venn_table[, -1]
	venn_table[is.na(venn_table)] <- 0
	# create a new microtable object
	tmp <- microtable$new(otu_table = venn_table)
	tmp
}
