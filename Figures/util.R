import_sparra_expr <- function(f) {
  sink_marker <- "***************"
  f2 <- readChar(f, file.info(f)$size)
  parse(text = strsplit(gsub("\r", "", f2), sink_marker, fixed = TRUE)[[1]][2])
}
