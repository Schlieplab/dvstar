# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install dvstar
#
# You can edit this file again by typing:
#
#     spack edit dvstar
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Dvstar(CMakePackage):
    """dvstar is a package for learning variable-length Markov chains (VLMC) from k-mers"""

    # Add a proper url for your package's homepage here.
    homepage = "https://github.com/Schlieplab/dvstar"
    #url = "https://github.com/Schlieplab/dvstar"

    # Add a list of GitHub accounts to
    # notify when the package is updated.
    maintainers("alexander@schlieplab.org")

    # Add the SPDX identifier of the project's license below.
    # See https://spdx.org/licenses/ for a list. Upon manually verifying
    # the license, set checked_by to your Github username.
    license("GPL-3.0-only", checked_by="alexander@schlieplab.org")

    git = "https://github.com/Schlieplab/dvstar.git"
    # Add proper versions and checksums here.
    version("main", branch="main", submodules=True)
    
    # FIXME: Add dependencies if required.
    depends_on("cmake@3.29")
    depends_on("intel-tbb@=2021.12.0")
    
    def install(self, spec, prefix):
        # FIXME: Unknown build system
        print("DVSTAR install", self.prefix, prefix, spec)
        mkdir(prefix.bin)
        src = self.build_directory
        install(join_path(src, 'dvstar'), join_path(prefix.bin, 'dvstar'))
