#' Calculate properties of sub-networks selected according to features in samples
#'
#' @description
#' Extracting sub-network according to the presence of features in each sample across networks and calculate the sub-network properties.
#' 
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{trans_network} class of \code{microeco} package.
#' @param ... parameters pass to \code{subset_network} function of \code{trans_network} class in microeco package (except the \code{node} parameter).
#' @return \code{data.frame}
#' @examples
#' \donttest{
#' data(soil_amp_network)
#' test <- subnet_property(soil_amp_network)
#' }
#' @export
subnet_property <- function(network_list, ...){
	check_input(network_list)
	res_property <- data.frame()
	for(i in names(network_list)){
		tmp <- data.frame()
		# extract the feature table used for network
		tmp_abund <- network_list[[i]]$data_abund %>% t %>% as.data.frame
		# also check all the nodes
		tmp_i_allnodes <- igraph::V(network_list[[i]]$res_network) %>% names
		for(j in colnames(tmp_abund)){
			tmp1 <- clone(network_list[[i]])
			tmp_nodes <- tmp_abund %>% .[.[, j] != 0, ] %>% rownames %>% .[. %in% tmp_i_allnodes]
			tmp1$res_network <- tmp1$subset_network(node = tmp_nodes, ...)
			suppressMessages(tmp1$cal_network_attr())
			tmp <- rbind(tmp, c(i, j, tmp1$res_network_attr[, 1, drop = TRUE]))
		}
		colnames(tmp) <- c("Network", "Sample", rownames(tmp1$res_network_attr))
		res_property <- rbind(res_property, tmp)
	}
	# rownames(res_property) <- res_property[, 1]
	res_property[, 3:ncol(res_property)] %<>% lapply(., as.numeric)
	res_property
}
