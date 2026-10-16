set(crescent_lib_name crescent_srclib)
add_library(${crescent_lib_name} OBJECT)

include(${CMAKE_SOURCE_DIR}/cmake/BuildCRESCENT.cmake)
build_crescent_lib(${crescent_lib_name})

add_executable(crescent "")
build_crescent_exe(crescent ${crescent_lib_name})
