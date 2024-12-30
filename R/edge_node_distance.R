#' @title Perform the distance distribution of paired nodes in edges across networks.
#'
#' @description
#' This class is a wrapper for a series of analysis on the distance values 
#' of paired nodes in edges across networks, including distance matrix conversion, the differential test and the visualization.
#'
#' @export
edge_node_distance <- R6::R6Class(classname = "edge_node_distance",
	public = list(
		#' @param network_list a list with multiple networks; all the networks should be \code{trans_network} object 
		#' 	 created from \code{trans_network} class of \code{microeco} package.
		#' @param dis_matrix default NULL; the distance matrix of nodes, used for the value extraction; 
		#' 	 must be a symmetrical matrix (or data.frame object) with both colnames and rownames (i.e. feature names).
		#' @param label default "+"; "+" or "-" or \code{c("+", "-")}; the edge label used for the selection of edges.
		#' @param with_module default FALSE; whether show the module classification of nodes in the result.
		#' @param module_thres default 2; the threshold of the nodes number of modules remained when \code{with_module = TRUE}.
		#' @return \code{data_table}, stored in the object
		#' @examples
		#' \donttest{
		#' data(soil_amp_network)
		#' data(soil_amp)
		#' # filter useless features to speed up the calculation
		#' node_names <- unique(unlist(lapply(soil_amp_network, function(x){colnames(x$data_abund)})))
		#' filter_soil_amp <- microeco::clone(soil_amp)
		#' filter_soil_amp$otu_table <- filter_soil_amp$otu_table[node_names, ]
		#' filter_soil_amp$tidy_dataset()
		#' # obtain phylogenetic distance matrix
		#' phylogenetic_distance <- as.matrix(cophenetic(filter_soil_amp$phylo_tree))
		#' # choose the positive labels
		#' t1 <- edge_node_distance$new(network_list = soil_amp_network, 
		#' 	 dis_matrix = phylogenetic_distance, label = "+")
		#' }
		initialize = function(network_list, dis_matrix = NULL, label = "+", with_module = FALSE, module_thres = 2){
			check_input(network_list)
			if(is.null(dis_matrix)){
				stop("Please provide dis_matrix parameter!")
			}
			if(!(inherits(dis_matrix, "matrix") | inherits(dis_matrix, "data.frame"))){
				stop("Input dis_matrix must be matrix or data.frame class!")
			}
			if(is.null(colnames(dis_matrix)) | is.null(rownames(dis_matrix))){
				stop("Input dis_matrix must have row names and column names!")
			}
			if(!any(colnames(dis_matrix) %in% rownames(dis_matrix))){
				stop("The colnames of dis_matrix must be same with the rownames! Please check the input dis_matrix!")
			}
			if(!is.logical(with_module)){
				stop("The parameter with_module must be logical!")
			}
			if(with_module){
				if(length(label) > 1){
					stop("The label parameter must be '+' or '-' when with_module = TRUE!")
				}
			}
			res_table <- data.frame()
			for(i in names(network_list)){
				net_obj <- network_list[[i]]
				if(!any(igraph::V(net_obj$res_network)$name %in% colnames(dis_matrix))){
					stop("The node names of network in ", i, " differ from the names in the input dis_matrix! Please check the input data!")
				}
				tmp <- private$get_matrix_value(
					network = net_obj, 
					label = label, 
					dis_matrix = dis_matrix, 
					group_name = i, 
					with_module = with_module, 
					module_thres = module_thres
					)
				res_table %<>% rbind(., tmp)
			}
			res_table$Group %<>% factor(., levels = names(network_list))
			res_table %<>% .[!is.na(.$Value), ]
			self$data_table <- res_table
			self$label <- label
			self$with_module <- with_module
		},
		#' @description
		#' Differential test across networks.
		#'
		#' @param method default "anova"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'anova'}}{Duncan's multiple range test for anova}
		#'     \item{\strong{'KW'}}{KW: Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons, see \code{dunnTest} function in \code{FSA} package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'   }
		#' @param ... parameters passed to \code{cal_diff} function of \code{trans_alpha} class of \code{microeco} package.
		#' @return \code{res_diff} in object. See the Return of \code{cal_diff} function in \code{trans_alpha} class of \code{microeco} package.
		#' @examples
		#' \donttest{
		#' t1$cal_diff(method = "wilcox")
		#' }
		cal_diff = function(method = c("anova", "KW", "KW_dunn", "wilcox", "t.test")[1], ...){
			res <- self$data_table
			res$Measure <- "Value"
			# two cases: only one type of label and two types of labels
			if(length(unique(res$label)) == 1){
				if(!self$with_module){
					suppressMessages(tmp2 <- trans_alpha$new(dataset = NULL))
					tmp2$data_alpha <- res
					tmp2$group <- "Group"
					tmp2$cal_diff(method = method, measure = "Value", ...)
				}else{
					res$raw_Group <- res$Group
					res$Module <- paste0(res$Group, " jointmark ", res$module)
					suppressMessages(tmp2 <- trans_alpha$new(dataset = NULL))
					tmp2$data_alpha <- res
					tmp2$group <- "Module"
					if(method != "anova"){
						message("For multiple labels, only anova can be used!")
					}
					tmp2$cal_diff(method = "anova", measure = "Value", ...)
					split_raw <- strsplit(tmp2$res_diff$Group, split = " jointmark ")
					tmp2$res_diff$by_group <- lapply(split_raw, function(x){x[1]}) %>% unlist
					tmp2$res_diff$Group <- lapply(split_raw, function(x){x[2]}) %>% unlist
					res$by_group <- res$raw_Group
					res$Module <- res$module
					res$Module %<>% factor(., levels = unique(.))
					tmp2$data_alpha <- res
					tmp2$by_group <- "by_group"
				}
			}else{
				res$raw_Group <- res$Group
				# jointmark instead of " - " or "&"
				res$Label <- paste0(res$Group, " jointmark ", res$label)
				suppressMessages(tmp2 <- trans_alpha$new(dataset = NULL))
				tmp2$data_alpha <- res
				tmp2$group <- "Label"
				if(method != "anova"){
					message("For multiple labels, only anova can be used!")
				}
				tmp2$cal_diff(method = "anova", measure = "Value", ...)
				split_raw <- strsplit(tmp2$res_diff$Group, split = " jointmark ")
				tmp2$res_diff$by_group <- lapply(split_raw, function(x){x[1]}) %>% unlist
				tmp2$res_diff$Group <- lapply(split_raw, function(x){x[2]}) %>% unlist
				res$by_group <- res$raw_Group
				res$Label <- res$label
				res$Label %<>% factor(., levels = self$label)
				tmp2$data_alpha <- res
				tmp2$by_group <- "by_group"
			}
			self$res_diff <- tmp2$res_diff
			self$cal_diff_method <- method
			self$tmp_diff <- tmp2
			message('The result is stored in object$res_diff ...')
			invisible(self)
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
		get_matrix_value = function(network, label, dis_matrix, group_name, with_module, module_thres){
			if(!with_module){
				if(is.null(network$res_edge_table)){
					network$get_edge_table()
				}
				tmp <- network$res_edge_table %>% .[.$label %in% label, ]
				if(nrow(tmp) == 0){
					res <- NA
				}else{
					select_value <- private$select_value_matrix(input_table = tmp, input_matrix = dis_matrix)
					res <- data.frame(Group = group_name, Value = select_value, label = tmp$label)
				}
			}else{
				if(!is.numeric(module_thres)){
					stop("The module_thres parameter must be numeric!")
				}
				if(is.null(network$res_node_table)){
					suppressMessages(network$get_node_table(node_roles = FALSE))
				}
				if(! "module" %in% colnames(network$res_node_table)){
					stop("please first use cal_module function to calculate modularity!")
				}
				# check module nodes number
				use_modules <- table(network$res_node_table$module) %>% .[. >= module_thres] %>% names
				res <- NULL
				for(k in use_modules){
					module_nodes <- network$res_node_table %>% .[.$module == k, ] %>% rownames
					t1 <- clone(network)
					t1$res_network <- t1$subset_network(node = module_nodes, rm_single = TRUE)
					suppressMessages(t1$get_edge_table())
					tmp <- t1$res_edge_table %>% .[.$label %in% label, ]
					if(nrow(tmp) == 0){
						next
					}else{
						select_value <- private$select_value_matrix(input_table = tmp, input_matrix = dis_matrix)
						res <- rbind(res, data.frame(Group = group_name, Value = select_value, module = k, label = tmp$label))
					}
				}
			}
			res
		},
		select_value_matrix = function(input_table, input_matrix){
			select_value <- lapply(seq_len(nrow(input_table)), function(x){
				if(all(c(input_table[x, 1], input_table[x, 2]) %in% colnames(input_matrix))){
					input_matrix[input_table[x, 1], input_table[x, 2]]
				}else{
					NA
				}
			}) %>% unlist
			select_value
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
