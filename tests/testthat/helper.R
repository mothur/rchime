# add common utility functions for testing here

# =========================================================================== #
# only used in testing
create_dummy_file <- function(filename, lines = c("test")) {
  remove_file(filename)

  file_conn <- file(filename)
  writeLines(lines, file_conn)
  close(file_conn)
}
# =========================================================================== #
