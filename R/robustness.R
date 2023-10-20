#' @title Calculate robustness across networks.
#'
#' @description
#' This class is a wrapper for robustness calculation and the visualization.
#'
#' @export
robustness <- R6::R6Class(classname = "robustness",
	public = list(
		#' @param network_list a list with multiple networks; all the networks should be \code{trans_network} object 
		#' 	 created from \code{\link{trans_network}} class of microeco package.
		#' @param link_remove_strategy default "Rand"; 
		#'   \describe{
		#' 	 item{\strong{"Rand"}} {links are randomly removed.}
		#'   }
		#' @param link_remove_perc default seq(0, 1, 0.1).
		#' @param measure default "Eff"; the network functioning measures as the representatives of robustness. 
		#'   \describe{
		#' 	 item{\strong{"Eff"}} {network efficiency.}
		
		#'   }
		#' @param run default 10.
		#' @return \code{res_table}, stored in the object
		initialize = function(network_list, link_remove_strategy = c("Rand")[1], link_remove_perc = seq(0, 1, 0.1), measure = c("Eff")[1], run = 10){
			check_input(network_list)
			res <- list()
			
			for(j in seq_along(network_list)){
				message("Network: ", names(network_list)[j], " ...")
				network <- network_list[[j]]$res_network
				total_edges_number <- length(igraph::E(network))
				for(i in link_remove_perc){
					message("Remove edges percentage: ", i, " ...")
					tmp_res <- c()
					if(i == 1){
						tmp_res_mean <- 0
						tmp_res_sd <- 0
					}else{
						for(k in seq_len(run)){
							if(link_remove_strategy == "Rand"){
								delete_edge_number <- sample(1:total_edges_number, size = round(total_edges_number * i))
							}
							network_del <- igraph::delete_edges(network, delete_edge_number)
							tmp_obj <- clone(network_list[[j]])
							tmp_obj$res_network <- network_del
							network_del <- tmp_obj$subset_network(node = V(network_del)$name, rm_single = TRUE)
							
							if(measure == "Eff"){
								tmp_res %<>% c(., private$measure_eff(network_del))
							}
						}
						tmp_res_mean <- mean(tmp_res)
						tmp_res_sd <- sd(tmp_res)
					}
					tmp_res_perc <- c(names(network_list)[j], i, link_remove_strategy, measure, tmp_res_mean, tmp_res_sd)
					res[[length(res) + 1]] <- tmp_res_perc
				}
			}
			res %<>% do.call(rbind, .) %>% as.data.frame
			colnames(res) <- c("network", "remove_perc", "remove_strategy", "measure", "mean", "sd")
			self$res_table <- res
			message("The result is stored in the object$res_table!")
		},
		#' @description
		#' Plot the distance.
		#'
		#' @param ... parameters pass to \code{plot_alpha} function of \code{trans_alpha} class of \code{microeco} package.
		#' @return \code{ggplot}.
		#' @examples
		#' \donttest{
		#' t1$plot(boxplot_add = "none", add_sig = TRUE)
		#' }
		plot = function(...){
			self$tmp_diff$plot_alpha(measure = "Value", ...)
		}
	),
	private = list(
		measure_eff = function(network){
			dis_matrix <- igraph::distances(network)
			nodes_num <- ncol(dis_matrix)
			dis_num <- as.dist(dis_matrix)
			res_single <- sum(1/dis_num)/(nodes_num * (nodes_num - 1))
			res_single
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
