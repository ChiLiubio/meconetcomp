#' Get node property tables from networks
#'
#' @description
#' Get node property tables from a list with multiple networks.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of microeco package.
#' @param ... parameter passed to get_node_table function of \code{\link{trans_network}} class.
#' @return list
#' @examples
#' \dontrun{
#' data(soil_amp_network)
#' soil_amp_network %<>% get_node_table(node_roles = FALSE)
#' }
#' @export
get_node_table <- function(network_list, ...){
	check_input(network_list)
	for(i in names(network_list)){
		cat(paste0("Run: ", i, " ...\n"))
		network_list[[i]]$get_node_table(...)
	}
	network_list
}
