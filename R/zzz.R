##' @importFrom utils packageDescription
.onAttach = function(
    libname,
    pkgname
){
  pkgVersion = packageDescription(pkgname, fields = "Version")
  msg = paste0("Welcome to infiltR package!
=======================================================================
", "You are using ", pkgname, " version ", pkgVersion,
"\n=======================================================================")
  packageStartupMessage(msg)
}
