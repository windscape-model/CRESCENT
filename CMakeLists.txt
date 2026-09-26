cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(CRESCENT CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)

option(CRESCENT_USE_INTERNAL_AMREX "Add AMReX as submodule" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if (CRESCENT_USE_INTERNAL_AMREX)
    set(AMREX_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")

    include(${CMAKE_SOURCE_DIR}/cmake/SetAmrexOptions.cmake)
    list(APPEND CMAKE_MODULE_PATH "${AMREX_LOCATION}/Tools/CMake")

    add_subdirectory(${AMREX_LOCATION})
else()
    find_package(AMReX REQUIRED)
endif()

find_package(MPI REQUIRED)

add_subdirectory(Exec)


set(CRESCENT_PUBLIC_HEADERS
    Source/shared/ProbBase.H
    Source/CRESCENT.H
    Source/external/CRESCENT_ExternalExchange.H
    Source/core/CRESCENT_Spectral.H
    Source/shared/CRESCENT_IO.H
    Source/shared/CRESCENT_Utils.H
    Source/init/CRESCENT_Init.H
    Source/core/CRESCENT_Physics.H
    Source/core/CRESCENT_Dynamics.H
    Source/shared/CRESCENT_Interpolation.H
    Source/nested/CRESCENT_NestedDomains.H
    Source/nested/CRESCENT_NestedMasks.H
    Source/nested/CRESCENT_NestedRK4.H
    Source/nested/CRESCENT_NestedEnergy.H)
set_target_properties(crescent_srclib PROPERTIES PUBLIC_HEADER "${CRESCENT_PUBLIC_HEADERS}")
