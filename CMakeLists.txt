# set minimum cmake version
cmake_minimum_required(VERSION 3.18)

# project name and language
project(fmm3d LANGUAGES C Fortran)

# verbose makefile
set(CMAKE_VERBOSE_MAKEFILE ON)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Grab Python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
# Grab OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# build static lib libfmm3d.a
# source files for libfmm3d.a
file(GLOB_RECURSE source_list "src/*.f" "src/*.f90")
#message(${source_list})
# remove fast kernels and tree_lr_3d.f
list(FILTER source_list EXCLUDE REGEX ".*_fast\\.f$|.*tree_lr_3d\\.f")
#message(${source_list})
# add fmm3d static lib
add_library(fmm3d STATIC ${source_list})
#compiler options
target_compile_options(fmm3d PRIVATE -fPIC -O3 -march=native -funroll-loops -std=legacy -w)

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

# extensions variables
# hfmm3d_fortran
set(hfmm_module_name "hfmm3d_fortran")
set(hfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap_vec.f
                          ${CMAKE_SOURCE_DIR}/src/Helmholtz/helmkernels.f)
set(f2py_helm_module_c "${hfmm_module_name}module.c")
# lfmm3d_fortran
set(lfmm_module_name "lfmm3d_fortran")
set(lfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap_vec.f
                          ${CMAKE_SOURCE_DIR}/src/Laplace/lapkernels.f)
set(f2py_lap_module_c "${lfmm_module_name}module.c")
# emfmm3d_fortran
set(emfmm_module_name "emfmm3d_fortran")
set(emfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap.f
                           ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap_vec.f
                           ${CMAKE_SOURCE_DIR}/src/Helmholtz/helmkernels.f
                           ${CMAKE_SOURCE_DIR}/src/Maxwell/emfmm3d.f90)
set(f2py_em_module_c "${emfmm_module_name}module.c")
# stfmm3d_fortran
set(stfmm_module_name "stfmm3d_fortran")
set(stfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap.f
                           ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap_vec.f
                           ${CMAKE_SOURCE_DIR}/src/Laplace/lapkernels.f
                           ${CMAKE_SOURCE_DIR}/src/Stokes/stfmm3d.f
                           ${CMAKE_SOURCE_DIR}/src/Stokes/stokkernels.f)
set(f2py_st_module_c "${stfmm_module_name}module.c")

# Generate extensions' sources
# hfmm3d_fortran
add_custom_target(
  hfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_helm_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_helm_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${hfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "hfmm3d_fortran"
                   ${hfmm_fortran_src_file}
                   --lower
  DEPENDS ${hfmm_fortran_src_file} # Fortran source
)
# lfmm3d_fortran
add_custom_target(
  lfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_lap_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_lap_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${lfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "lfmm3d_fortran"
                   ${lfmm_fortran_src_file}
                   --lower
  DEPENDS ${lfmm_fortran_src_file} # Fortran source
)
# emfmm3d_fortran
add_custom_target(
  emfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_em_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_em_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${emfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "emfmm3d_fortran"
                   ${emfmm_fortran_src_file}
                   --lower
  DEPENDS ${emfmm_fortran_src_file} # Fortran source
)
# stfmm3d_fortran
add_custom_target(
  stfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_st_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_st_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${stfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "stfmm3d_fortran"
                   ${stfmm_fortran_src_file}
                   --lower
  DEPENDS ${stfmm_fortran_src_file} # Fortran source
)

# Set up extensions targets
# hfmm3d_fortran
Python_add_library(hfmm3d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_helm_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${hfmm_fortran_src_file}" # Fortran source(s)
)
# lfmm3d_fortran
Python_add_library(lfmm3d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_lap_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${lfmm_fortran_src_file}" # Fortran source(s)
)
# emfmm3d_fortran
Python_add_library(emfmm3d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_em_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${emfmm_fortran_src_file}" # Fortran source(s)
)
# stfmm3d_fortran
Python_add_library(stfmm3d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_st_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${stfmm_fortran_src_file}" # Fortran source(s)
)

# Dependencies for extensions
# hfmm3d_fortran
target_link_libraries(hfmm3d_fortran PRIVATE Python::NumPy fmm3d)
add_dependencies(hfmm3d_fortran hfmm_genpyf)
target_include_directories(hfmm3d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")
# lfmm3d_fortran
target_link_libraries(lfmm3d_fortran PRIVATE Python::NumPy fmm3d)
add_dependencies(lfmm3d_fortran lfmm_genpyf)
target_include_directories(lfmm3d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")
# emfmm3d_fortran
target_link_libraries(emfmm3d_fortran PRIVATE Python::NumPy fmm3d)
add_dependencies(emfmm3d_fortran emfmm_genpyf)
target_include_directories(emfmm3d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")
# stfmm3d_fortran
target_link_libraries(stfmm3d_fortran PRIVATE Python::NumPy fmm3d)
add_dependencies(stfmm3d_fortran stfmm_genpyf)
target_include_directories(stfmm3d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")

add_subdirectory(python)
