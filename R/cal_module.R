#' Assign modules to each network
#'
#' @description
#' Calculating modularity of networks and assign the modules to nodes for each network.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{trans_network} class of \code{microeco} package.
#' @param undirected_method default "cluster_fast_greedy"; the modularity algorithm for undirected network; 
#'       see \code{cal_module} function of \code{trans_network} class for more algorithms.
#' @param directed_method default 'cluster_optimal'; the modularity algorithm for directed network.
#' @param ... other parameters (except for method) passed to \code{cal_module} function of \code{trans_network} class.
#' @return \code{list}, with module attribute in nodes of each network
#' @examples
#' data(soil_amp_network)
#' soil_amp_network <- cal_module(soil_amp_network)
#' 
#' @export
cal_module <- function(network_list, undirected_method = "cluster_fast_greedy", directed_method = 'cluster_optimal', ...){
	check_input(network_list)
	for(i in names(network_list)){
		if(is_directed(network_list[[i]]$res_network)){
			network_list[[i]]$cal_module(method = directed_method, ...)
		}else{
			network_list[[i]]$cal_module(method = undirected_method, ...)
		}
	}
	network_list
}
