import_sparra_expr <- function(f) {
  sink_marker <- "***************"
  f2 <- readChar(f, file.info(f)$size)
  parse(text = strsplit(gsub("\r", "", f2), sink_marker, fixed = TRUE)[[1]][2])
}

build_diff <- function(df, xvar, yvar, baseline = "v3") {
  models <- levels(df$Model)

  xrange <- range(df |> pull({{ xvar }}))

  df2 <- tibble("{{ xvar }}" := seq(xrange[1], xrange[2], length = 100)[2:99])

  for(m in models) {
    df2 <- df2 |>
      mutate("{m}" := suppressWarnings(approx(df |> filter(Model == m) |> pull({{xvar}}),
                                              df |> filter(Model == m) |> pull({{yvar}}),
                                              {{xvar}})$y))
  }

  df2
}
