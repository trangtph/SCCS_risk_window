###########################################
### Project: SCCS Risk Window Scanning ####
### Author: Trang - Azida              ####
###########################################

###################################
### Helper functions            ###
###################################

log_error <- function(msg,
                      output_dir = here::here("Report")) {
  
  tryCatch({
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    log_file <- file.path(output_dir, "error_log.txt")
    
    cat(
      sprintf("[%s] %s\n", Sys.time(), msg),
      file = log_file,
      append = TRUE
    )
  }, error = function(e) {
    ## last-resort fallback: never stop execution
    message("Logging failed: ", conditionMessage(e))
  })
}

