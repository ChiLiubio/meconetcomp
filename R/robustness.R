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
			remove_strategy = c("edge_rand", "edge_strong", "edge_weak", "node_hub", "node_degree")[1], 
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

			for(j in seq_along(network_list)){
				message("Network: ", names(network_list)[j], " ...")
				network <- network_list[[j]]$res_network

				for(i in remove_ratio){
					message("Remove ratio: ", i, " ...")
					network_del <- list()
					
					if(any(grepl("^edge_", remove_strategy))){
						delete_edge_sequence <- list()
						total_edges_number <- length(E(network))
						delete_number <- round(total_edges_number * i)
						if("edge_rand" %in% remove_strategy){
							tmp_delete_edge_sequence <- replicate(run, sample(1:total_edges_number, size = delete_number), simplify = FALSE)
							delete_edge_sequence[[paste0("edge_rand_number_", delete_number)]] <- tmp_delete_edge_sequence
						}
						if("edge_strong" %in% remove_strategy){
							tmp_delete_edge_sequence <- order(E(network)$weight, decreasing = TRUE)[1:delete_number] %>% list
							delete_edge_sequence[[paste0("edge_strong_number_", delete_number)]] <- tmp_delete_edge_sequence
						}
						if("edge_weak" %in% remove_strategy){
							tmp_delete_edge_sequence <- order(E(network)$weight, decreasing = FALSE)[1:delete_number] %>% list
							delete_edge_sequence[[paste0("edge_weak_number_", delete_number)]] <- tmp_delete_edge_sequence
						}
						tmp_network_del_list <- lapply(delete_edge_sequence, function(y){
							lapply(y, function(x){
								tmp_network <- igraph::delete_edges(network, x)
								tmp_obj <- clone(network_list[[j]])
								tmp_obj$res_network <- tmp_network
								tmp_network <- tmp_obj$subset_network(node = V(tmp_network)$name, rm_single = TRUE)
								tmp_network
							})
						})
						network_del %<>% c(., tmp_network_del_list)
					}
					if(any(grepl("^node_", remove_strategy))){
						delete_node_names <- list()
						if("node_hub" %in% remove_strategy){
							tmp_obj <- clone(network_list[[j]])
							total_hub_names <- tmp_obj$res_node_table %>% .[grepl("hub", .$taxa_roles), ] %>% rownames
							if(length(total_hub_names) == 0){
								warning("No network or module hub is found for the network: ", names(network_list)[j], "! The node_hub results come from the full network!")
							}
							delete_number <- round(length(total_hub_names) * i)
							tmp_delete_node_names <- replicate(run, sample(total_hub_names, size = delete_number), simplify = FALSE)
							delete_node_names[[paste0("node_hub_number_", delete_number)]] <- tmp_delete_node_names
						}
						if("node_degree" %in% remove_strategy){
							total_nodes_number <- length(V(network))
							delete_number <- round(total_nodes_number * i)
							node_names_order <- sort(degree(network), decreasing = TRUE) %>% names
							tmp_delete_node_names <- node_names_order[1:delete_number] %>% list
							delete_node_names[[paste0("node_degree_number_", delete_number)]] <- tmp_delete_node_names
						}
						tmp_network_del_list <- lapply(delete_node_names, function(y){
							lapply(y, function(x){
								igraph::delete_vertices(network, x)
							})
						})
						network_del %<>% c(., tmp_network_del_list)
					}
					
					tmp_res <- list()
					if("Eff" %in% measure){
						tmp_res[["Eff"]] <- private$measure_eff(network_del)
					}
					
					tmp_res_table <- lapply(names(tmp_res), function(x){
						tmp_table <- sapply(tmp_res[[x]], function(y){
							y %<>% unlist
							Mean <- mean(y)
							SD <- sd(y)
							c(Mean = Mean, SD = SD)
						}) %>% t %>% as.data.frame
						tmp_table$remove_strategy <- names(tmp_res[[x]])
						tmp_table$measure <- x
						tmp_table
					})
					tmp_res_table %<>% do.call(rbind, .)
					tmp_res_table$remove_ratio <- i
					tmp_res_table$Network <- names(network_list)[j]
					res[[length(res) + 1]] <- tmp_res_table
				}
			}
			res %<>% do.call(rbind, .) %>% as.data.frame %>% microeco::dropallfactors(unfac2num = TRUE)
			res %<>%  tidyr::separate_wider_delim(., cols = "remove_strategy", delim = "_number_", names = c("remove_strategy", "remove_number"))
			res %<>% .[, c("Network", "remove_strategy", "remove_ratio", "remove_number", "measure", "Mean", "SD")] %>% as.data.frame
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
				p <- p + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), linewidth = errorbar_width, position = errorbar_position, size = errorbar_size)
			}
			if(add_fitting == T){
				p <- p + geom_smooth(se = FALSE, ...)
			}else{
				p <- p + geom_line(...)
			}
			p <- p + theme_bw()
			if(length(unique(res_table$remove_strategy)) > 1 | length(unique(res_table$measure)) > 1){
				p <- p + facet_grid(remove_strategy ~ measure, drop = TRUE, scale = "free", space = "fixed")
			}
			
			p
		}
	),
	private = list(
		measure_eff = function(all_networks){
			lapply(all_networks, function(y){
				lapply(y, function(x){
					if(length(V(x)) < 2 | length(E(x)) < 2){
						0
					}else{
						dis_matrix <- igraph::distances(x)
						nodes_num <- ncol(dis_matrix)
						dis_num <- as.dist(dis_matrix)
						res_single <- sum(1/dis_num)/(nodes_num * (nodes_num - 1))
						res_single
					}
				})
			})
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
