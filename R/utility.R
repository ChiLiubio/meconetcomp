
# Check whether the input list has the correct format
check_input <- function(input){
	if(!is.list(input)){
		stop("The input data must a list! Please check the data input!")
	}
	for(i in names(input)){
		if(!inherits(input[[i]], "trans_network")){
			stop(i, " in the input list is not a trans_network object! Please check it!")
		}
	}
}

