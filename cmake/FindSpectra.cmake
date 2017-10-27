# - Try to find Spectra lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Spectra)
#
# Once done this will define
#

include(FindPackageHandleStandardArgs)

find_path(SPECTRA_INCLUDE_DIR
    NAMES   SymEigsSolver.h
    PATHS   /Users/npignet/Documents/Outils/spectra/include
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SPECTRA DEFAULT_MSG SPECTRA_INCLUDE_DIR)
