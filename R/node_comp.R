#' Generate a microtable object with node distributions across networks
#'
#' @description
#' Generate a microtable object with node distributions across networks. Useful for the node information comparisons across different networks.
#' 
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of \code{microeco} package.
#' @param property default "name"; a colname of \code{res_node_table} in each network; 
#'    the default "name" represents using node presence/absence information in the otu_table of final output, in which
#'    1 represents presence of the node in the corresponding network; 
#'    For other options (such as degree), the results in the output otu_table are the actual values of \code{res_node_table}.
#' @return \code{microtable} object
#' @examples
#' \donttest{
#' data(soil_amp_network)
#' test <- node_comp(soil_amp_network)
#' # test is a microtable object
#' }
#' @export
node_comp <- function(network_list, property = "name"){
	check_input(network_list)
	for(i in names(network_list)){
		if(is.null(network_list[[i]]$res_node_table)){
			suppressMessages(network_list[[i]]$get_node_table(node_roles = FALSE))
		}
		if(! property %in% colnames(network_list[[i]]$res_node_table)){
			stop("Input property is not a colname of res_node_table in network ", i, " !")
		}
		if(property == "name"){
			tmp <- network_list[[i]]$res_node_table %>% .[, property, drop = FALSE]
			tmp[, i] <- 1
		}else{
			tmp <- network_list[[i]]$res_node_table %>% .[, c("name", property), drop = FALSE]
			colnames(tmp)[2] <- i
		}
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
	tmp <- suppressMessages(microtable$new(otu_table = venn_table))
	tmp
}
