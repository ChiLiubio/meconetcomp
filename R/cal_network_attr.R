#' Calculate network topological property for each network
#'
#' @description
#' Calculate the topological properties of all the networks and merge the results into one table.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{trans_network} class of \code{microeco} package.
#' @return \code{data.frame}
#' @examples
#' data(soil_amp_network)
#' test <- cal_network_attr(soil_amp_network)
#' 
#' @export
cal_network_attr <- function(network_list){
	check_input(network_list)
	# first create a list network_property to place in each result
	network_property <- list()
	for(i in names(network_list)){
		suppressMessages(network_list[[i]]$cal_network_attr())
		network_property[[i]] <- network_list[[i]]$res_network_attr
	}
	# finally merge all the results of the list into a table
	network_property <- do.call(cbind, network_property)
	network_property
}
