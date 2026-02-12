###########################################
### Project: SCCS Risk Window Scanning ####
### Author: Trang - Azida              ####
###########################################

###################################
### Helper functions            ###
###################################


# Log errors ----------------------------------------------------------------

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

# Create the folders to store results ---------------------------------------
create_directory <- function(path) {
  if (!dir.exists(path)) 
    dir.create(path, recursive = TRUE)
}

# Saving results to .csv file -----------------------------------------------
append_to_csv <- function(out_df, file_path) {
  if (!file.exists(file_path)) {
    write.csv(out_df, file_path, row.names = FALSE)
  } else {
    write.table(out_df, file_path, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE)
  }
}




