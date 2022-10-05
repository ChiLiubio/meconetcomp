#' Extract subset of network according to the edge intersection of networks
#'
#' @description
#' Extracting a network according to the edge intersection of networks.
#' 
#' @param network_list a list with multiple networks; all the networks should be \code{trans_network} object created from \code{\link{trans_network}} class of microeco package.
#' @param venn default NULL; a \code{microtable} object which must be converted by \code{trans_comm} function of \code{trans_venn} class.
#' @param name default NULL; integer or character; must be a number or one of colnames of the \code{otu_table} in the input \code{venn} parameter.
#' @return \code{trans_network} object, with only the extracted edges in the network
#' @examples
#' \donttest{
#' data(soil_amp_network)
#' # first obtain edge distribution
#' tmp <- edge_comp(soil_amp_network)
#' # obtain edge intersection using trans_venn class
#' tmp1 <- microeco::trans_venn$new(tmp)
#' # convert intersection result to microtable object
#' tmp2 <- tmp1$trans_comm()
#' # extract the intersection of all the three networks ("IW", "TW" and "CW")
#' test <- subset_network(soil_amp_network, venn = tmp2, name = "IW&TW&CW")
#' # test is a trans_network object
#' }
#' @export
subset_network <- function(network_list, venn = NULL, name = NULL){
	# check the input
	check_input(network_list)
	if(is.null(venn)){
		stop("Please provide venn parameter!")
	}else{
		if(!inherits(venn, "microtable")){
			stop("The input venn must a microtable object! Please check the input!")
		}
	}
	if(is.null(name)){
		stop("Please provide name parameter!")
	}else{
		if(length(name) > 1){
			stop("The input name has a length > 1! Please check the input!")
		}
		if(is.numeric(name)){
			if(name > ncol(venn$otu_table)){
				stop("Input name is ", name, " , but venn$otu_table only have ", ncol(venn$otu_table), " columns!")
			}else{
				name <- colnames(venn$otu_table)[name]
			}
		}else{
			if(is.character(name)){
				if(! (name %in% colnames(venn$otu_table))){
					stop("The input name must be one of sample names of input venn!")
				}
			}
		}
	}
	features <- venn$otu_table %>% .[.[, name] > 0, , drop = FALSE] %>% rownames
	# identify the sample name
	if(grepl("&", name)){
		name_split <- strsplit(name, "&") %>% unlist
		name_split_use <- name_split[1]
	}else{
		name_split_use <- name_split <- name
	}
	if(name_split_use %in% names(network_list)){
		network_tmp1 <- network_list[[name_split_use]]
	}else{
		stop("The name ", name_split_use, " is not in names of input network_list!")
	}
	
	tmp2 <- network_tmp1$res_network
	tmp2_edge_table <- igraph::as_data_frame(tmp2, what = "edges")
	# make the names of paired nodes ordered
	edge_nodes <- tmp2_edge_table[, 1:2] %>% t %>% as.data.frame
	sorted_edge_nodes <- lapply(edge_nodes, function(x){sort(x) %>% paste0(., collapse = " -- ")}) %>% unlist
	sub_network <- igraph::delete_edges(tmp2, which(!(sorted_edge_nodes %in% features)))
	igraph::E(sub_network)$weight <- 1
	# delete nodes without edges
	nodes_raw <- igraph::V(sub_network)$name
	edges <- igraph::as_data_frame(sub_network, what = "edges")
	delete_nodes <- nodes_raw %>% .[! . %in% as.character(c(edges[,1], edges[,2]))]
	if(length(delete_nodes) > 0){
		sub_network %<>% igraph::delete_vertices(delete_nodes)
	}
	nodes_use <- igraph::V(sub_network)$name
	
	res <- clone(network_tmp1)
	res$res_network <- sub_network
	
	# abund use mean values
	data_abund_merge <- data.frame()
	for(i in name_split){
		use_network <- network_list[[i]]
		data_abund_merge %<>% rbind(apply(use_network$data_abund[, nodes_use], 2, sum), .)
	}
	colnames(data_abund_merge) <- nodes_use
	rownames(data_abund_merge) <- name_split
	
	res$data_abund <- data_abund_merge
	res$res_cor_p <- NULL
	# return trans_network
	res
}
