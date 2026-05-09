# Copyright (C) 2026 Frank Koopmans
#
# This file is part of pairscale.
#
# pairscale is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pairscale is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


#' pairscale package declaration
#' @name pairscale
#' @useDynLib pairscale, .registration=TRUE
#' @keywords internal
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

#' cleanup Rcpp code
#'
#' @param libpath library path
#' @noRd
.onUnload = function(libpath = NULL) {
  library.dynam.unload("pairscale", libpath)
}
