# Set source root before moving to a separate working-data directory.
source_root <- normalizePath(Sys.getenv("HCC_SOURCE_ROOT", "."), mustWork=TRUE)
workspace <- normalizePath(Sys.getenv("HCC_WORKSPACE", source_root), mustWork=TRUE)
setwd(workspace)
dir.create("outputs", showWarnings=FALSE)
r_library <- Sys.getenv("HCC_R_LIBRARY", "")
if (nzchar(r_library)) .libPaths(c(normalizePath(r_library, mustWork=TRUE), .libPaths()))
jags_modules <- Sys.getenv("HCC_JAGS_MODULE_DIR", "")
if (nzchar(jags_modules)) options(jags.moddir=normalizePath(jags_modules, mustWork=TRUE))
