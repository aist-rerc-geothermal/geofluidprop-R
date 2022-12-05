
# printf <- function(...)
#   cat(sprintf(...))

create_initial_Rd_filles_for_swig_functions = function()
{
  swig_generated_R_filename = "R/geofluidprop.R"

  # extract function names
  o = grep("= function", readLines(swig_generated_R_filename), value = TRUE)
  fnames = mapply(function(s) gsub("`(.+)` =.*", "\\1", s), o)

  printf("found %d functions in %s\n", length(fnames), swig_generated_R_filename)

  #print(o)

  for (fname in fnames)
  {
    filename = file.path("man", sprintf("%s.Rd", fname))
    if (file.exists(filename))
      next
    
    lines = c(sprintf("\\name{%s}", fname))
    lines = c(lines, sprintf("\\alias{%s}", fname))
    lines = c(lines, sprintf("\\title{%s}", fname))
    lines = c(lines, sprintf("\\description{%s}", fname))
    fileConn <- file(filename)
    writeLines(lines, fileConn)
    close(fileConn)
  }
}

extend_NAMESPACE = function()
{
  filename = "NAMESPACE"
  str = "exportPattern(\"^[[:alpha:]]+\")\""
  write(str,file=filename,append=TRUE)
  str = "useDynLib(geofluidprop)"
  write(str,file=filename,append=TRUE)
}
