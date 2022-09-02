#' Get edge property tables from networks
#'
#' @description
#' Get edge property tables from a list with multiple networks.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of microeco package.
#' @return list
#' @examples
#' \dontrun{
#' data(soil_amp_network)
#' soil_amp_network %<>% get_edge_table
#' }
#' @export
get_edge_table <- function(network_list){
	check_input(network_list)
	for(i in names(network_list)){
		cat(paste0("Run: ", i, " ...\n"))
		network_list[[i]]$get_edge_table()
	}
	network_list
}
