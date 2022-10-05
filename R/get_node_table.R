#' Get node property table for each network
#'
#' @description
#' Get node property table for each network in the list with multiple networks.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of \code{microeco} package.
#' @param ... parameter passed to get_node_table function of \code{\link{trans_network}} class.
#' @return \code{list}, with \code{res_node_table} in each network
#' @examples
#' \donttest{
#' data(soil_amp_network)
#' soil_amp_network <- get_node_table(soil_amp_network, node_roles = FALSE)
#' }
#' @export
get_node_table <- function(network_list, ...){
	check_input(network_list)
	for(i in names(network_list)){
		suppressMessages(network_list[[i]]$get_node_table(...))
	}
	network_list
}
