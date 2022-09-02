#' Check the input data format
#'
#' @param input default NULL; character or data.frame; matching table used.
check_input <- function(input = NULL){
	if(!is.list(input)){
		stop("The input data must a list! Please check the data input!")
	}
	for(i in names(input)){
		if(!inherits(input[[i]], "trans_network")){
			stop(i, " in the input list is not a trans_network object! Please check it!")
		}
	}
}
