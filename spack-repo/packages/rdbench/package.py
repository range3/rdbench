# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Rdbench(CMakePackage):
    """2D reaction-diffusion system benchmark using MPI and MPI-IO."""

    homepage = "https://github.com/range3/rdbench"
    url      = "https://github.com/range3/rdbench/archive/v0.1.0-1.tar.gz"
    git      = "https://github.com/range3/rdbench.git"

    maintainers = ['range3']

    version('master', branch='master')
    version('0.1.1', sha256='ea74c1b96b660352814038b779c361792a3ae068461db84a48dc8d11b01edff1', preferred=True)
    version('0.1.0-1', sha256='5239c1702df9dbdca99cef50966e8624d310f7641b8a8405571cadb948e2de15')

    depends_on('mpi')

    def setup_build_environment(self, env):
        env.unset('CPM_SOURCE_CACHE')

    def cmake_args(self):
        args = []
        return args
