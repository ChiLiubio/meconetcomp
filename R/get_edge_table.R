#' Get edge property table for each network
#'
#' @description
#' Get edge property table for each network in the list with multiple networks.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of microeco package.
#' @return list, with res_edge_table in each network
#' @examples
#' data(soil_amp_network)
#' soil_amp_network <- get_edge_table(soil_amp_network)
#' 
#' @export
get_edge_table <- function(network_list){
	check_input(network_list)
	for(i in names(network_list)){
		suppressMessages(network_list[[i]]$get_edge_table())
	}
	network_list
}
