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
    version('0.5.1', sha256='43febfa8325143eb068ea25ba9d3f8994e48de4a1b08a8fe45e95c13a5a27a99', preferred=True)
    version('0.4.2', sha256='e832530bb874410c3893f2568ee9f80d8cdc25ac281c31ca9ac42c306cd746a8')
    version('0.4.1', sha256='359a90211b0f07b98e1319366469a44985734f83fc1d39f2df6a0871bcf7183e')
    version('0.4.0', sha256='b9b64a8097b37f7c031ad4fbedbe7e086fd61307d5be616d30fcf8e9802045f5')
    version('0.3.0', sha256='f0e3c86557e44aca44361f061fa88ea4553de2520bf2369c2979e0fc0ae1500e')
    version('0.2.0', sha256='b621c659dff6bb7d6e1d8d04cf7485707b44059746e2fede5534b7b51a784767')
    version('0.1.1', sha256='ea74c1b96b660352814038b779c361792a3ae068461db84a48dc8d11b01edff1')
    version('0.1.0-1', sha256='5239c1702df9dbdca99cef50966e8624d310f7641b8a8405571cadb948e2de15')

    depends_on('mpi')

    def setup_build_environment(self, env):
        env.unset('CPM_SOURCE_CACHE')

    def cmake_args(self):
        args = []
        return args
