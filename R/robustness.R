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
		#' @param remove_strategy default "edge_rand"; 
		#'   \describe{
		#' 	 item{\strong{"edge_rand"}} {links are randomly removed.}
		#' 	 item{\strong{"edge_strong"}} {links are removed in decreasing order of weight.}
		#' 	 item{\strong{"edge_weak"}} {links are removed in increasing order of weight.}
		#'   }
		#' @param remove_ratio default seq(0, 1, 0.1).
		#' @param measure default "Eff"; the network functioning measures as the representatives of robustness. 
		#'   \describe{
		#' 	 item{\strong{"Eff"}} {network efficiency.}
		
		#'   }
		#' @param run default 10. Replication number applied for the sampling method.
		#' @return \code{res_table}, stored in the object
		initialize = function(network_list, 
			remove_strategy = c("edge_rand", "edge_strong", "edge_weak", "node_hub")[1], 
			remove_ratio = seq(0, 1, 0.1), 
			measure = c("Eff")[1], 
			run = 10
			){

			require(igraph)
			check_input(network_list)
			if("node_hub" %in% remove_strategy){
				private$check_node_table(network_list)
			}
			res <- list()

			for(n in measure){
				message("measure: ", n, " ...")
				for(m in remove_strategy){
					message("remove_strategy: ", m, " ...")
					for(j in seq_along(network_list)){
						message("Network: ", names(network_list)[j], " ...")
						network <- network_list[[j]]$res_network

						for(i in remove_ratio){
							message("Remove ratio: ", i, " ...")

							if(grepl("edge_", m)){
								total_edges_number <- length(E(network))
								delete_number <- round(total_edges_number * i)
								if(m == "edge_rand"){
									delete_edge_sequence <- replicate(run, sample(1:total_edges_number, size = delete_number), simplify = FALSE)
								}
								if(m == "edge_strong"){
									delete_edge_sequence <- order(E(network)$weight, decreasing = TRUE)[1:delete_number] %>% list
								}
								if(m == "edge_weak"){
									delete_edge_sequence <- order(E(network)$weight, decreasing = FALSE)[1:delete_number] %>% list
								}
								network_del_list <- lapply(delete_edge_sequence, function(x){
									tmp_network <- igraph::delete_edges(network, x)
									tmp_obj <- clone(network_list[[j]])
									tmp_obj$res_network <- tmp_network
									tmp_network <- tmp_obj$subset_network(node = V(tmp_network)$name, rm_single = TRUE)
									tmp_network
								})
								
							}else{
								if(m == "node_hub"){
									tmp_obj <- clone(network_list[[j]])
									total_hub_names <- tmp_obj$res_node_table %>% .[grepl("hub", .$taxa_roles), ] %>% rownames
									if(length(total_hub_names) == 0){
										warning("No network or module hub is found for the network: ", names(network_list)[j], "! The results come from the full network!")
									}
									delete_number <- round(length(total_hub_names) * i)
									delete_node_names <- replicate(run, sample(total_hub_names, size = delete_number), simplify = FALSE)
								}
								network_del_list <- lapply(delete_node_names, function(x){
									igraph::delete_vertices(network, x)
								})
							}
							
							if(n == "Eff"){
								tmp_res <- private$measure_eff(network_del_list)
							}
							
							tmp_res_mean <- mean(tmp_res)
							tmp_res_sd <- sd(tmp_res)
							tmp_res_perc <- c(names(network_list)[j], i, m, n, tmp_res_mean, tmp_res_sd, delete_number)
							res[[length(res) + 1]] <- tmp_res_perc
						}
					}
				}
			}
			res %<>% do.call(rbind, .) %>% as.data.frame %>% microeco::dropallfactors(unfac2num = TRUE)
			colnames(res) <- c("Network", "remove_ratio", "remove_strategy", "measure", "Mean", "SD", "remove_number")
			res[, 1] %<>% factor(., levels = names(network_list))
			self$res_table <- res
			message("The result is stored in the object$res_table!")
		},
		#' @description
		#' Plot the simulation results.
		#'
		#' @param color_values colors used for presentation.
		#' @param show_point default TRUE; whether show the point.
		#' @param point_size default .3; point size value.
		#' @param point_alpha default .6; point alpha value.
		#' @param show_errorbar default TRUE; whether show the errorbar by using the SD result.
		#' @param errorbar_position default position_dodge(0); Position adjustment, either as a string (such as "identity"), 
		#' 	 or the result of a call to a position adjustment function.
		#' @param errorbar_size default 1; errorbar size.
		#' @param errorbar_width default 0.1; errorbar width.
		#' @param add_fitting default FALSE; whether add fitted line.
		#' @param ... parameters pass to ggplot2::geom_line (when add_fitting = FALSE) or ggplot2::geom_smooth (when add_fitting = TRUE).
		#' @return \code{ggplot}.
		#' @examples
		#' \donttest{
		#' t1$plot()
		#' }
		plot = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			show_point = TRUE,
			point_size = 1,
			point_alpha = .6,
			show_errorbar = TRUE,
			errorbar_position = position_dodge(0),
			errorbar_size = 1,
			errorbar_width = 0.1,
			add_fitting = FALSE,
			...
			){
			res_table <- self$res_table
			
			p <- ggplot(res_table, aes(x = remove_ratio, y = Mean, color = Network)) +
				scale_color_manual(values = color_values)
			if(show_point){
				p <- p + geom_point(alpha = point_alpha, size = point_size)
			}
			if(show_errorbar){
				p <- p + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = errorbar_width, position = errorbar_position, size = errorbar_size)
			}
			if(add_fitting == T){
				p <- p + geom_smooth(se = FALSE, ...)
			}else{
				p <- p + geom_line(...)
			}
			p <- p + theme_bw()
			
			p
		}
	),
	private = list(
		measure_eff = function(all_networks){
			lapply(all_networks, function(x){
				if(length(V(x)) < 2 | length(E(x)) < 2){
					0
				}else{
					dis_matrix <- igraph::distances(x)
					nodes_num <- ncol(dis_matrix)
					dis_num <- as.dist(dis_matrix)
					res_single <- sum(1/dis_num)/(nodes_num * (nodes_num - 1))
					res_single
				}
			}) %>% unlist
		},
		check_node_table = function(network_list){
			for(i in names(network_list)){
				if(is.null(network_list[[i]]$res_node_table)){
					network_list[[i]]$get_node_table()
				}
			}
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
