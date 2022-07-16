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
    version('0.1.0-1', sha256='5239c1702df9dbdca99cef50966e8624d310f7641b8a8405571cadb948e2de15')

    depends_on('mpi')

    def cmake_args(self):
        args = []
        return args
