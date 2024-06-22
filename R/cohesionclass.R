#' Calculate the cohesion of samples for each network
#'
#' @description
#' The cohesion is a method for quantifying the connectivity of microbial communities <doi:10.1038/ismej.2017.91>.
#' It is defined:
#'    \deqn{C_{j}^{pos} = \sum_{i=1}^{n} a_{i} \cdot \bar{r_{i}}_{|r>0}}
#'    \deqn{C_{j}^{neg} = \sum_{i=1}^{n} a_{i} \cdot \bar{r_{i}}_{|r<0}}
#'  where \eqn{C_{j}^{pos}} is the positive cohesion, and \eqn{C_{j}^{neg}} is the negative cohesion.
#'  \eqn{a_{i}} is the relative abundance of species i in sample j.
#'  \eqn{\bar{r_{i}}_{|r>0}} denotes the mean weight (correlation coefficient, interaction strength) of all the edges (related with species i) with positive association.
#'  
#' 
#' @export
cohesionclass <- R6::R6Class(classname = "cohesionclass",
	public = list(
		#' @param network_list a list with multiple networks; all the networks should be \code{trans_network} object 
		#' 	 created from \code{trans_network} class of microeco package.
		#' @return \code{res_list}, stored in the object.
		#'  It includes two tables: res_feature and res_sample. In res_feature, the r_pos and r_neg columns mean the \eqn{\bar{r_{i}}_{|r>0}} and \eqn{\bar{r_{i}}_{|r<0}}.
		#'  In res_sample, the c_pos and c_neg columns denote \eqn{C_{j}^{pos}} and \eqn{C_{j}^{neg}}.
		#' @examples
		#' t1 <- cohesionclass$new(soil_amp_network)
		#' 
		initialize = function(network_list){
			check_input(network_list)
			res_feature <- list()
			res_sample <- list()
			
			for(j in seq_along(network_list)){
				network_name <- names(network_list)[j]
				message("Network: ", network_name, " ...")
				network <- network_list[[j]]
				
				abund_table <- network$data_abund
				rel_abund <- apply(t(abund_table), 2, function(x){x/sum(x)})
				
				if(is.null(network$res_edge_table)){
					suppressMessages(network$get_edge_table())
				}
				edge_nodes <- network$res_edge_table
				if(! "weight" %in% colnames(edge_nodes)){
					stop("The weight column in the res_edge_table is not found!")
				}
				edge_nodes$weight <- ifelse(edge_nodes$label == "+", 1, -1) * edge_nodes$weight
				features <- unique(c(edge_nodes[, 1], edge_nodes[, 2]))
				tmp_res <- lapply(features, function(x){
					tmp_edge <- edge_nodes[edge_nodes[, 1] == x | edge_nodes[, 2] == x, ]
					c(r_pos = mean(tmp_edge[tmp_edge[, "weight"] > 0, "weight"]), r_neg = mean(tmp_edge[tmp_edge[, "weight"] < 0, "weight"]))
				})
				tmp_res %<>% do.call(rbind, .)
				tmp_res_feature <- data.frame(network = network_name, feature = features, tmp_res)
				tmp_res_feature$r_pos[is.nan(tmp_res_feature$r_pos)] <- 0
				tmp_res_feature$r_neg[is.nan(tmp_res_feature$r_neg)] <- 0

				# community
				tmp_res <- lapply(colnames(rel_abund), function(x){
					c(c_pos = sum(rel_abund[tmp_res_feature$feature, x] * tmp_res_feature$r_pos), c_neg = sum(rel_abund[tmp_res_feature$feature, x] * tmp_res_feature$r_neg))
				})
				tmp_res %<>% do.call(rbind, .)
				tmp_res_sample <- data.frame(network = network_name, sample = colnames(rel_abund), tmp_res)

				res_feature[[j]] <- tmp_res_feature
				res_sample[[j]] <- tmp_res_sample
			}
			res_feature %<>% do.call(rbind, .)
			res_sample %<>% do.call(rbind, .)
			res <- list(feature = res_feature, sample = res_sample)
			self$res_list <- res
			message("The result is stored in the object$res_list!")
		},
		#' @description
		#' Differential test.
		#'
		#' @param measure default "c_pos"; "c_pos" or "c_neg" in the \code{res_list$sample}; "r_pos" or "r_neg" in the \code{res_list$feature}.
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
		cal_diff = function(measure = "c_pos", method = c("anova", "KW", "KW_dunn", "wilcox", "t.test")[1], ...){
			measure <- match.arg(measure, c("c_pos", "c_neg", "r_pos", "r_neg"))
			if(measure %in% c("c_pos", "c_neg")){
				data_table <- self$res_list$sample
			}else{
				data_table <- self$res_list$feature
			}
			tmp2 <- private$prepare_data(data_table, measure)
			
			tmp2$cal_diff(method = method, measure = measure, ...)
			
			self$res_diff <- tmp2$res_diff
			self$cal_diff_method <- method
			self$tmp_diff <- tmp2
			message('The result is stored in object$res_diff ...')
		},
		#' @description
		#' Plot the result.
		#'
		#' @param measure default "c_pos"; "c_pos" or "c_neg" in the \code{res_list$sample}; "r_pos" or "r_neg" in the \code{res_list$feature}.
		#' @param ... parameters pass to \code{plot_alpha} function of \code{trans_alpha} class of \code{microeco} package.
		#' @return \code{ggplot}.
		#' @examples
		#' \donttest{
		#' t1$plot(boxplot_add = "none", add_sig = TRUE)
		#' }
		plot = function(measure = "c_pos", ...){
			if(is.null(self$tmp_diff)){
				measure <- match.arg(measure, c("c_pos", "c_neg", "r_pos", "r_neg"))
				if(measure %in% c("c_pos", "c_neg")){
					data_table <- self$res_list$sample
				}else{
					data_table <- self$res_list$feature
				}
				tmp <- private$prepare_data(data_table, measure)
			}else{
				tmp <- self$tmp_diff
				measure <- t1$tmp_diff$measure
			}
			tmp$plot_alpha(measure = measure, ...)
		}
	),
	private = list(
		prepare_data = function(data_table, measure){
			data_table$Measure <- measure
			colnames(data_table)[colnames(data_table) == measure] <- "Value"
			suppressMessages(tmp2 <- trans_alpha$new(dataset = NULL))
			tmp2$data_alpha <- data_table
			tmp2$group <- "network"
			tmp2
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
