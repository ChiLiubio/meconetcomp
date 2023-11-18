#' @title Calculate robustness across networks.
#'
#' @description
#' This class is a wrapper for robustness calculation and visualization.
#'
#' @export
robustness <- R6::R6Class(classname = "robustness",
	public = list(
		#' @param network_list a list with multiple networks; all the networks should be \code{trans_network} object 
		#' 	 created from \code{\link{trans_network}} class of microeco package.
		#' @param remove_strategy default "edge_rand"; 
		#'   \describe{
		#' 	   \item{\strong{"edge_rand"}}{edges are randomly removed.}
		#'     \item{\strong{"edge_strong"}}{edges are removed in decreasing order of weight.}
		#' 	   \item{\strong{"edge_weak"}}{edges are removed in increasing order of weight.}
		#' 	   \item{\strong{"node_rand"}}{nodes are removed randomly.}
		#' 	   \item{\strong{"node_hub"}}{node hubs are removed. The hubs include network hubs and module hubs.}
		#' 	   \item{\strong{"node_degree_high"}}{nodes are removed in decreasing order of degree.}
		#' 	   \item{\strong{"node_degree_low"}}{nodes are removed in increasing order of degree.}
		#'   }
		#' @param remove_ratio default seq(0, 1, 0.1).
		#' @param measure default "Eff"; network robustness measures. 
		#'   \describe{
		#' 	   \item{\strong{"Eff"}}{network efficiency. The average efficiency of the network is defined:
		#' 	         \deqn{Eff = \frac{1}{N(N - 1)} \sum_{i \neq j \in G}\frac{1}{d(i, j)}}
		#' 	   	 where N is the total number of nodes and \emph{d(i,j)} is the shortest path between node i and node j. 
		#' 	   	 When the weight is found in the edge attributes, \eqn{d(i,j)} denotes the weighted shortest path between node i and node j.
		#' 	   	 For more details, please read the references <doi: 10.1007/s11704-016-6108-z> and <doi: 10.1038/s41598-020-60298-7>.
		#' 	   	 }
		#' 	   \item{\strong{"Eigen"}}{natural connectivity <doi: 10.1007/s11704-016-6108-z>.
		#' 	   	 The natural connectivity can be regarded as an average eigenvalue that changes strictly monotonically with the addition or deletion of edges. It is defined:
		#' 	   	     \deqn{\bar{\lambda} = \ln(\frac{1}{N} \sum_{i=1}^{N} e^{\lambda~i~})}
		#' 	   	 where \eqn{\lambda~i~} is the \eqn{i}th eigenvalue of the graph adjacency matrix. The larger the value of \eqn{\bar{\lambda}} is, the more robust the network is.
		#' 	   	 }
		#' 	   \item{\strong{"Pcr"}}{critical removal fraction of vertices (edges) for the disintegration of networks 
		#' 	   	 <doi: 10.1007/s11704-016-6108-z> <doi: 10.1103/PhysRevE.72.056130>.
		#' 	   	 This is a robustness measure based on random graph theory.
		#' 	   	 The critical fraction against random attacks is labeled as \eqn{P_{c}^r}. It is defined:
		#' 	   	     \deqn{P_{c}^r = 1 - \frac{1}{\frac{\langle k^2 \rangle}{\langle k \rangle} - 1}}
		#' 	   	 where \eqn{\langle k \rangle} is the average nodal degree of the original network, and \eqn{\langle k^2 \rangle} is the average of square of nodal degree. 
		#' 	   	 }
		#'   }
		#' @param run default 10. Replication number applied for the sampling method.
		#' @return \code{res_table}, stored in the object.
		#' @examples
		#' tmp <- robustness$new(soil_amp_network, remove_strategy = c("edge_rand"), 
		#'   measure = c("Eff"), run = 3, remove_ratio = c(0.1, 0.5, 0.9))
		#' 
		initialize = function(network_list, 
			remove_strategy = c("edge_rand", "edge_strong", "edge_weak", "node_rand", "node_hub", "node_degree_high", "node_degree_low")[1], 
			remove_ratio = seq(0, 1, 0.1), 
			measure = c("Eff", "Eigen", "Pcr")[1], 
			run = 10
			){
			
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
						total_edges_number <- length(igraph::E(network))
						delete_number <- round(total_edges_number * i)
						if("edge_rand" %in% remove_strategy){
							tmp_delete_edge_sequence <- replicate(run, sample(1:total_edges_number, size = delete_number), simplify = FALSE)
							delete_edge_sequence[[paste0("edge_rand_number_", delete_number)]] <- tmp_delete_edge_sequence
						}
						if("edge_strong" %in% remove_strategy){
							tmp_delete_edge_sequence <- order(igraph::E(network)$weight, decreasing = TRUE)[1:delete_number] %>% list
							delete_edge_sequence[[paste0("edge_strong_number_", delete_number)]] <- tmp_delete_edge_sequence
						}
						if("edge_weak" %in% remove_strategy){
							tmp_delete_edge_sequence <- order(igraph::E(network)$weight, decreasing = FALSE)[1:delete_number] %>% list
							delete_edge_sequence[[paste0("edge_weak_number_", delete_number)]] <- tmp_delete_edge_sequence
						}
						tmp_network_del_list <- lapply(delete_edge_sequence, function(y){
							lapply(y, function(x){
								tmp_network <- igraph::delete_edges(network, x)
								nodes_raw <- igraph::V(tmp_network)$name
								edges <- igraph::as_data_frame(tmp_network, what = "edges")
								delete_nodes <- nodes_raw %>% .[! . %in% as.character(c(edges[,1], edges[,2]))]
								if(length(delete_nodes) > 0){
									tmp_network %<>% igraph::delete_vertices(delete_nodes)
								}
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
						if(any(c("node_rand", "node_degree_high", "node_degree_low") %in% remove_strategy)){
							total_nodes_number <- length(igraph::V(network))
							delete_number <- round(total_nodes_number * i)
						}
						if("node_rand" %in% remove_strategy){
							total_node_names <- igraph::V(network)$name
							tmp_delete_node_names <- replicate(run, sample(total_node_names, size = delete_number), simplify = FALSE)
							delete_node_names[[paste0("node_rand_number_", delete_number)]] <- tmp_delete_node_names
						}
						if("node_degree_high" %in% remove_strategy){
							node_names_order <- sort(igraph::degree(network), decreasing = TRUE) %>% names
							delete_node_names[[paste0("node_degree_high_number_", delete_number)]] <- node_names_order[1:delete_number] %>% list
						}
						if("node_degree_low" %in% remove_strategy){
							node_names_order <- sort(igraph::degree(network), decreasing = FALSE) %>% names
							delete_node_names[[paste0("node_degree_low_number_", delete_number)]] <- node_names_order[1:delete_number] %>% list
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
						tmp_res[["Eff"]] <- measure_eff(network_del)
					}
					if("Eigen" %in% measure){
						tmp_res[["Eigen"]] <- private$measure_eigen(network_del)
					}
					if("Pcr" %in% measure){
						ori_degree <- igraph::degree(network) %>% mean
						tmp_res[["Pcr"]] <- private$measure_pcr(network_del, ori_degree)
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
		#' @param add_fitting default FALSE; whether add fitted smooth line. FALSE denotes add line segment among points.
		#' @param ... parameters pass to ggplot2::geom_line (when add_fitting = FALSE) or ggplot2::geom_smooth (when add_fitting = TRUE).
		#' @return \code{ggplot}.
		#' @examples
		#' \donttest{
		#' tmp$plot(linewidth = 1)
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
				p <- p + facet_grid(measure ~ remove_strategy, drop = TRUE, scale = "free", space = "fixed")
			}
			
			p
		}
	),
	private = list(
		measure_eigen = function(all_networks){
			lapply(all_networks, function(y){
				lapply(y, function(x){
					if(length(igraph::V(x)) < 2 | length(igraph::E(x)) < 2){
						0
					}else{
						am <- igraph::as_adjacency_matrix(x)
						check_res <- tryCatch(pca <- princomp(am), error = function(e) {skip_to_next <- TRUE})
						if(rlang::is_true(check_res)) {
							NA
						}else{
							ev <- (pca$sdev)^2
							log(sum(exp(ev))/length(ev), 10)
						}
					}
				})
			})
		},
		measure_pcr = function(all_networks, od){
			lapply(all_networks, function(y){
				lapply(y, function(x){
					if(length(igraph::V(x)) < 2 | length(igraph::E(x)) < 2){
						NA
					}else{
						all_degree <- igraph::degree(x)
						ksq <- mean(all_degree^2)
						k0 <- ksq/od
						if(k0 == 1){
							NA
						}else{
							1 - (1 / (k0 - 1))
						}
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
