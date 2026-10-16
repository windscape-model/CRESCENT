cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(CRESCENT CXX C)

include(CTest)

# GPU (CUDA) enablement. The CUDA language must be enabled at project scope,
# before AMReX is added and the crescent targets are created, so that
# setup_target_for_cuda_compilation() in cmake/BuildCRESCENT.cmake can mark the
# .cpp sources as CUDA. Off by default (CPU build unaffected).
if(CRESCENT_ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_EXTENSIONS OFF)
endif()

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

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()


