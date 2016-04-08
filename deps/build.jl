### build.jl --- Build script for Cuba.jl.

# Copyright (C) 2016  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: numeric integration

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

### Code:

local_dir  = "cuba-shared-object"
local_file = local_dir*".tar.gz"
libcuba = joinpath(dirname(@__FILE__), "libcuba")
@linux? (object="libcuba.so") : (@osx? (object="libcuba.dylib") : (object=""))

# Clean already existing shared object, archive and buil directory in order to
# perform a new clean build.
run(`rm -rf $local_file $local_dir libcuba.so`)

# Download Cuba and build the shared object.
info("Downloading Cuba source...")
download("https://github.com/giordano/cuba/archive/shared-object.tar.gz",
         local_file)
run(`tar xzf $local_file`)
info("Building libcuba...")
cd(local_dir) do
    run(`./configure`)
    run(`make shared`)
    run(`mv $object ..`)
end

# Make sure Julia is able to see the library.
if length(Libdl.find_library([libcuba])) > 0
    info("libcuba successfully installed!")
else
    error("Installation of libcuba failed")
end
