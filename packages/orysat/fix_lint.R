fix_file <- function(path) {
  lines <- readLines(path, warn = FALSE)

  # 1. Remove explicit return(): return(x) -> x
  lines <- gsub("^(\\s*)return\\((.*)\\)\\s*$", "\\1\\2", lines)

  # 2. seq_linter: 1:length(x) -> seq_along(x)
  lines <- gsub("1:length\\(([^)]+)\\)", "seq_along(\\1)", lines)
  # seq_linter: 1:ncol(x) -> seq_len(ncol(x))
  lines <- gsub("1:ncol\\(([^)]+)\\)", "seq_len(ncol(\\1))", lines)
  # seq_linter: 1:nrow(x) -> seq_len(nrow(x))
  lines <- gsub("1:nrow\\(([^)]+)\\)", "seq_len(nrow(\\1))", lines)

  # 3. vector_logic_linter: & -> && and | -> || inside if() conditions
  lines <- gsub("if \\((.+) & (.+)\\)", "if (\\1 && \\2)", lines)
  lines <- gsub("if \\((.+) \\| (.+)\\)", "if (\\1 || \\2)", lines)

  # 4. brace_linter: 'repeat{' -> 'repeat {'
  lines <- gsub("repeat\\{", "repeat {", lines)

  # 5. brace_linter: opening brace on its own line -> merge with previous line
  result <- c()
  for (i in seq_along(lines)) {
    if (length(result) > 0 && grepl("^\\s*\\{\\s*$", lines[i])) {
      result[length(result)] <- paste0(result[length(result)], " {")
    } else {
      result <- c(result, lines[i])
    }
  }
  lines <- result

  writeLines(lines, path)
}

files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in files) {
  cat("Fixing:", basename(f), "\n")
  fix_file(f)
}
cat("Done\n")
