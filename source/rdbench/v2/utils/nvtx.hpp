#pragma once

#ifdef RDBENCH_USE_NVTX
#include <nvtx3/nvtx3.hpp>
#define RDBENCH_NVTX_RANGE(name) nvtx3::scoped_range _nvtx_range_(name)
#define RDBENCH_NVTX_MARK(name) nvtx3::mark(name)
#define RDBENCH_NVTX_FUNC_RANGE() NVTX3_FUNC_RANGE()
#else
#define RDBENCH_NVTX_RANGE(name) ((void)0)
#define RDBENCH_NVTX_MARK(name) ((void)0)
#define RDBENCH_NVTX_FUNC_RANGE() ((void)0)
#endif
