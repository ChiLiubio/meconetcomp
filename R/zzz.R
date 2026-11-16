
.onAttach <- function(...) {
	withr::with_preserve_seed({
		# Show messages only in interactive sessions to prevent interference during non-interactive checks (e.g., R CMD check)
		if (!interactive()) {
			return()
		}
		
		tips <- c(
			"Use MecoToolkit (https://mecotoolkit.com/) for easier network building and comparison analysis."
		)
		
		tip <- tips[1]
		packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
	})
}
