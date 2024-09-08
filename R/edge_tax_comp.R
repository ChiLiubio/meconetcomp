#' Taxonomic sum of linked nodes in edges across networks
#'
#' @description
#' Taxonomic sum of linked nodes in edges across networks.
#'
#' @param network_list a list with multiple networks; all the networks should be trans_network object created from \code{trans_network} class of \code{microeco} package.
#' @param taxrank default "Phylum"; Which taxonomic level is used for the sum of nodes in edges.
#' @param label default "+"; "+" or "-" or \code{c("+", "-")}; the edge label used for the selection of edges.
#' @param rel default \code{TRUE}; \code{TRUE} represents using ratio, the denominator is the number of selected edges; 
#'   \code{FALSE} represents the absolute number of the sum of edges.
#' @param sep default " -- "; The separator for two taxonomic names shown in the result.
#' @return \code{data.frame}
#' @examples
#' data(soil_amp_network)
#' test <- edge_tax_comp(soil_amp_network)
#' # test is a microtable object
#' 
#' @export
edge_tax_comp <- function(network_list, taxrank = "Phylum", label = "+", rel = TRUE, sep = " -- "){
	check_input(network_list)
	source_compare <- NULL
	for(i in names(network_list)){
		suppressMessages(network_list[[i]]$cal_sum_links(taxa_level = taxrank))
		if(label == "+"){
			tmp1 <- network_list[[i]]$res_sum_links_pos
		}else{
			if(label == "-"){
				tmp1 <- network_list[[i]]$res_sum_links_neg
			}else{
				if(all(c("+", "-") %in% label)){
					tmp1 <- rbind(network_list[[i]]$res_sum_links_pos, network_list[[i]]$res_sum_links_neg)
				}else{
					stop("Unknown label parameter input!")
				}
			}
		}
		if(is.null(tmp1)){
			next
		}else{
			# convert the result to long format
			tmp2 <- reshape2::melt(tmp1, value.name = i)
			tmp3 <- tmp2[, 1:2] %>% t %>% as.data.frame
			tmp2$name <- lapply(tmp3, function(x){sort(x) %>% paste0(., collapse = sep)}) %>% unlist
			tmp2 <- tmp2[, c("name", i)]
			# remove duplicates
			tmp2 %<>% .[!duplicated(.), ]
			if(is.null(source_compare)){
				source_compare <- tmp2
			}else{
				source_compare <- full_join(source_compare, tmp2, by = c("name" = "name"))
			}
		}
	}
	source_compare[is.na(source_compare)] <- 0
	# calculate the ratio with the edge number as denominator
	if(rel){
		for(i in names(network_list)){
			if(i %in% colnames(source_compare)){
				if(is.null(network_list[[i]]$res_edge_table)){
					suppressMessages(network_list[[i]]$get_edge_table())
				}
				tmp1 <- network_list[[i]]$res_edge_table
				tmp1 %<>% .[.$label %in% label, ]
				source_compare[, i] %<>% {. / nrow(tmp1)}
			}else{
				next
			}
		}
	}
	rownames(source_compare) <- source_compare[, 1]
	source_compare %<>% .[, -1]
	source_compare
}
