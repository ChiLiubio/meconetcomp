#' Save networks as gexf-style files
#'
#' @description
#' Save all networks as gexf-style files, which can be opened by Gephi (\href{https://gephi.org/}{https://gephi.org/}).
#' 
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{\link{trans_network}} class of microeco package.
#' @param dirpath default "."; the directory used to save all the gexf files; the default "." represents the current working directory.
#' @return NULL.
#' @export
save_network <- function(network_list, dirpath = "."){
	check_input(network_list)
	if(!dir.exists(dirpath)){
		dir.create(dirpath)
	}
	for(i in names(network_list)){
		cat(paste0("Run: ", i, " ...\n"))
		network_list[[i]]$save_network(filepath = file.path(dirpath, paste0(i, ".gexf")))
	}
}
