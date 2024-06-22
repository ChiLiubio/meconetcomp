#' Extract subset of network according to the edge intersection of networks
#'
#' @description
#' Extracting a network according to the edge intersection of networks.
#' 
#' @param network_list a list with multiple networks; all the networks should be \code{trans_network} object created from \code{trans_network} class of microeco package.
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
		# use the first as the starting network
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
	sorted_edge_nodes <- get_sorted_edge_name(tmp2) %>% rownames
	sub_network <- igraph::delete_edges(tmp2, which(!(sorted_edge_nodes %in% features)))
	# delete nodes without edges in sub_network
	nodes_raw <- igraph::V(sub_network)$name
	edges <- igraph::as_data_frame(sub_network, what = "edges")
	delete_nodes <- nodes_raw %>% .[! . %in% as.character(c(edges[,1], edges[,2]))]
	if(length(delete_nodes) > 0){
		sub_network %<>% igraph::delete_vertices(delete_nodes)
	}
	if(length(name_split) > 1){
		# assign weight with the mean weight across related network
		network_list_related <- list()
		for(i in name_split){
			network_list_related[[i]] <- network_list[[i]]
		}
		if(all(unlist(lapply(network_list_related, function(x){!is.null(igraph::E(x$res_network)$weight)})))){
			# first generate a list to store all the edge tables
			all_edge_tables <- list()
			for(j in names(network_list_related)){
				tmp_table <- get_sorted_edge_name(network_list_related[[j]]$res_network)
				all_edge_tables[[j]] <- tmp_table[features, "weight", drop = FALSE]
				colnames(all_edge_tables[[j]])[1] <- paste0("weight_", j)
			}
			all_edge_tables %<>% do.call(cbind, .)
			# get sub_network edge table
			sub_network_edge_table <- get_sorted_edge_name(sub_network)
			igraph::E(sub_network)$weight <- apply(all_edge_tables, 1, mean)[rownames(sub_network_edge_table)]
			message('The weight of each edge in extracted network is the mean across networks ...')
		}else{
			igraph::E(sub_network)$weight <- 1
			message('Part of networks have no weight attribute of edges. Use 1 as the final weight for all edges in extracted network ...')
		}
	}

	nodes_use <- igraph::V(sub_network)$name
	# select a trans_network object and change files in it
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


# inner function
# make the names of paired nodes in edge ordered
get_sorted_edge_name <- function(network){
	tmp1 <- igraph::as_data_frame(network, what = "edges")
	tmp2 <- tmp1[, 1:2] %>% t %>% as.data.frame
	tmp3 <- lapply(tmp2, function(x){sort(x) %>% paste0(., collapse = " -- ")}) %>% unlist
	rownames(tmp1) <- tmp3
	# output table
	tmp1
}

