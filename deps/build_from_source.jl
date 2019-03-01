### build.jl --- Build script for Cuba.jl.

# Copyright (C) 2016-2019  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>

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

using Libdl, BinaryProvider

tagfile = "installed_version"

@static if Sys.isunix()
    # SHA hash of the revision to be downloaded from
    # https://github.com/giordano/cuba/tree/julia
    sha="11bfbf509088f168622b8268f49c0a59ee81758b"
    if !isfile(tagfile) || readchomp(tagfile) != sha
        local_dir  = "cuba-$(sha)"
        local_file = "$(local_dir).tar.gz"
        if Sys.isapple()
            object="libcuba.dylib"
        else
            object="libcuba.so"
        end

        # Clean already existing tag file, shared object, and build directory in
        # order to perform a new clean build.  Leave the archive, to avoid
        # downloading it again.
        run(`rm -rf $tagfile $object $local_dir`)

        if !isfile(local_file)
            @info("Downloading Cuba source...")
            download("https://github.com/giordano/cuba/archive/$(sha).tar.gz",
                     local_file)
        end
        run(`tar xzf $local_file`)
        @info("Building libcuba...")
        cd(local_dir) do
            run(`./configure`)
            run(`make shared`)
            run(`mv $object ..`)
        end
        open(tagfile, "w") do file
            println(file, sha)
        end
    end
else
    error("Cannot build CUBA from source for your platform $(triplet(platform_key_abi())).
Try setting ENV[\"JULIA_CUBA_BUILD_SOURCE\"]=\"false\" and building again the package")
end

products = [LibraryProduct(@__DIR__, String["libcuba"], :libcuba)]
write_deps_file(joinpath(@__DIR__, "deps.jl"), products)
