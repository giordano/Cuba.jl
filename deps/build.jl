### build.jl --- Build script for Cuba.jl.

# Copyright (C) 2016  Mosè Giordano

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

using Compat

libcuba = joinpath(dirname(@__FILE__), "libcuba")
tagfile = "installed_version"

if is_windows()
    version="v4.2-20150925" # Version of Cuba library.
    arch=Sys.ARCH
    if !isfile(tagfile) || readdlm(tagfile) != [version, arch]
        info("Downloading Cuba shared library...")
        download("https://bintray.com/giordano/Cuba-Shared-Library/download_file?file_path=libcuba-$(version)-$(Sys.ARCH).dll",
                 "libcuba.dll")
        writedlm(tagfile, [version, arch])
    end
else # GNU/Linux and Mac OS
    # SHA hash of the revision to be downloaded from
    # https://github.com/giordano/cuba/tree/julia
    sha="44da524f3fb9e5a803614c19cbfe7027a372601e"
    if !isfile(tagfile) || readchomp(tagfile) != sha
        local_dir  = "cuba-$(sha)"
        local_file = "$(local_dir).tar.gz"
        if is_linux()
            object="libcuba.so"
        elseif is_apple()
            object="libcuba.dylib"
        else
            object=""
        end

        # Clean already existing tag file, shared object, and build directory in
        # order to perform a new clean build.  Leave the archive, to avoid
        # downloading it again.
        run(`rm -rf $tagfile $object $local_dir`)

        if !isfile(local_file)
            info("Downloading Cuba source...")
            download("https://github.com/giordano/cuba/archive/$(sha).tar.gz",
                     local_file)
        end
        run(`tar xzf $local_file`)
        info("Building libcuba...")
        cd(local_dir) do
            run(`./configure`)
            run(`make shared`)
            run(`mv $object ..`)
        end
        open(tagfile, "w") do file
            println(file, sha)
        end
    end
end

# Make sure Julia is able to see the library.
if isempty(Libdl.find_library([libcuba]))
    if isfile(tagfile)
        # Delete the tagfile in order to force building upon next
        # Pkg.build("Cuba").  Note that Julia 0.4 doesn't have "force" keyword
        # to `rm' function.
        rm(tagfile)
    end
    error("Installation of libcuba failed")
else
    info("libcuba successfully installed!")
end
