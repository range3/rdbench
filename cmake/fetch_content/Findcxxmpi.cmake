include(FetchContent)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.28")
    FetchContent_Declare(
        cxxmpi
        GIT_REPOSITORY https://github.com/range3/cxxmpi.git
        GIT_TAG master
        GIT_REMOTE_UPDATE_STRATEGY REBASE
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(cxxmpi)
else()
    FetchContent_Declare(
        cxxmpi
        GIT_REPOSITORY https://github.com/range3/cxxmpi.git
        GIT_TAG master
        GIT_REMOTE_UPDATE_STRATEGY REBASE
    )
    FetchContent_GetProperties(cxxmpi)
    if(NOT cxxmpi_POPULATED)
        FetchContent_Populate(cxxmpi)
        add_subdirectory("${cxxmpi_SOURCE_DIR}" "${cxxmpi_BINARY_DIR}" EXCLUDE_FROM_ALL)
    endif()
endif()

set(cxxmpi_FOUND 1)
