# set minimum cmake version
cmake_minimum_required(VERSION 3.18)

# project name and language
project(fmm3d LANGUAGES C Fortran)

# Allow wheel builds to override CPU-specific tuning for portable binaries.
set(FMM3D_ARCH_FLAGS "-march=native" CACHE STRING "Architecture-specific compiler flags for FMM3D")
separate_arguments(FMM3D_ARCH_FLAGS_LIST NATIVE_COMMAND "${FMM3D_ARCH_FLAGS}")
option(FMM3D_FAST_KER "Build with optimized direct kernels and vectorized C++ kernels" OFF)

if(FMM3D_FAST_KER)
  enable_language(CXX)
endif()

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
set(fmm3d_openmp_components Fortran)
if(FMM3D_FAST_KER)
  list(APPEND fmm3d_openmp_components CXX)
endif()
find_package(OpenMP REQUIRED COMPONENTS ${fmm3d_openmp_components})

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
file(GLOB_RECURSE source_list CONFIGURE_DEPENDS "src/*.f" "src/*.f90")
# remove tree_lr_3d.f and select either the reference or fast direct kernels
list(FILTER source_list EXCLUDE REGEX ".*tree_lr_3d\\.f$")
if(FMM3D_FAST_KER)
  list(
    FILTER source_list EXCLUDE REGEX
      ".*/Helmholtz/helmkernels\\.f$|.*/Helmholtz/hndiv\\.f$|.*/Laplace/lapkernels\\.f$|.*/Laplace/lndiv\\.f$"
  )
  list(APPEND source_list "${CMAKE_SOURCE_DIR}/vec-kernels/src/libkernels.cpp")
else()
  list(FILTER source_list EXCLUDE REGEX ".*_fast\\.f$")
endif()

set(FMM3D_HELM_KERNEL_SOURCE "${CMAKE_SOURCE_DIR}/src/Helmholtz/helmkernels.f")
set(FMM3D_LAP_KERNEL_SOURCE "${CMAKE_SOURCE_DIR}/src/Laplace/lapkernels.f")
if(FMM3D_FAST_KER)
  set(FMM3D_HELM_KERNEL_SOURCE "${CMAKE_SOURCE_DIR}/src/Helmholtz/helmkernels_fast.f")
  set(FMM3D_LAP_KERNEL_SOURCE "${CMAKE_SOURCE_DIR}/src/Laplace/lapkernels_fast.f")
endif()

# add fmm3d static lib
add_library(fmm3d STATIC ${source_list})
#compiler options
target_compile_options(
  fmm3d
  PRIVATE
    -fPIC
    -O3
    -funroll-loops
    ${FMM3D_ARCH_FLAGS_LIST}
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy;-w>
    $<$<COMPILE_LANGUAGE:CXX>:-w>
)
target_link_libraries(fmm3d PUBLIC OpenMP::OpenMP_Fortran)
if(FMM3D_FAST_KER)
  target_compile_features(fmm3d PRIVATE cxx_std_11)
  target_compile_definitions(fmm3d PRIVATE SCTL_PROFILE=-1)
  target_include_directories(fmm3d PRIVATE "${CMAKE_SOURCE_DIR}/vec-kernels/include")
  target_link_libraries(fmm3d PUBLIC OpenMP::OpenMP_CXX)
endif()

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)
cmake_print_variables(FMM3D_FAST_KER)

# extensions variables
# hfmm3d_fortran
set(hfmm_module_name "hfmm3d_fortran")
set(hfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap_vec.f
                          ${FMM3D_HELM_KERNEL_SOURCE})
set(f2py_helm_module_c "${hfmm_module_name}module.c")
# lfmm3d_fortran
set(lfmm_module_name "lfmm3d_fortran")
set(lfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap_vec.f
                          ${FMM3D_LAP_KERNEL_SOURCE})
set(f2py_lap_module_c "${lfmm_module_name}module.c")
# emfmm3d_fortran
set(emfmm_module_name "emfmm3d_fortran")
set(emfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap.f
                           ${CMAKE_SOURCE_DIR}/src/Helmholtz/hfmm3dwrap_vec.f
                           ${FMM3D_HELM_KERNEL_SOURCE}
                           ${CMAKE_SOURCE_DIR}/src/Maxwell/emfmm3d.f90)
set(f2py_em_module_c "${emfmm_module_name}module.c")
# stfmm3d_fortran
set(stfmm_module_name "stfmm3d_fortran")
set(stfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap.f
                           ${CMAKE_SOURCE_DIR}/src/Laplace/lfmm3dwrap_vec.f
                           ${FMM3D_LAP_KERNEL_SOURCE}
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

set(FMM3D_PYTHON_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/fmm3dpy")
file(MAKE_DIRECTORY "${FMM3D_PYTHON_PACKAGE_DIR}")
configure_file("${CMAKE_SOURCE_DIR}/python/fmm3dpy/__init__.py" "${FMM3D_PYTHON_PACKAGE_DIR}/__init__.py" COPYONLY)
configure_file("${CMAKE_SOURCE_DIR}/python/fmm3dpy/fmm3d.py" "${FMM3D_PYTHON_PACKAGE_DIR}/fmm3d.py" COPYONLY)
foreach(fmm3d_python_target hfmm3d_fortran lfmm3d_fortran emfmm3d_fortran stfmm3d_fortran)
  set_target_properties(
    ${fmm3d_python_target}
    PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${FMM3D_PYTHON_PACKAGE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${FMM3D_PYTHON_PACKAGE_DIR}"
      ARCHIVE_OUTPUT_DIRECTORY "${FMM3D_PYTHON_PACKAGE_DIR}"
  )
endforeach()

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

include(CTest)
if(BUILD_TESTING)
  set(FMM3D_FORTRAN_TEST_COMPILE_OPTIONS -O3 -std=legacy -w)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    list(APPEND FMM3D_FORTRAN_TEST_COMPILE_OPTIONS -fallow-argument-mismatch)
  endif()

  function(fmm3d_add_fortran_driver_test target_name output_name source_file suite_label)
    set(output_dir "${CMAKE_CURRENT_BINARY_DIR}/test/${suite_label}")
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/ctest_runs/${suite_label}/${target_name}")
    file(MAKE_DIRECTORY "${output_dir}")
    file(MAKE_DIRECTORY "${run_dir}")

    add_executable(${target_name} "${source_file}")
    set_target_properties(
      ${target_name}
      PROPERTIES
        OUTPUT_NAME "${output_name}"
        RUNTIME_OUTPUT_DIRECTORY "${output_dir}"
    )
    target_compile_options(${target_name} PRIVATE ${FMM3D_FORTRAN_TEST_COMPILE_OPTIONS})
    target_link_libraries(${target_name} PRIVATE fmm3d)

    add_test(NAME "fortran_${suite_label}_${target_name}" COMMAND $<TARGET_FILE:${target_name}>)
    set_tests_properties(
      "fortran_${suite_label}_${target_name}"
      PROPERTIES
        LABELS "fortran;${suite_label}"
        WORKING_DIRECTORY "${run_dir}"
    )
  endfunction()

  fmm3d_add_fortran_driver_test(helm_helmrouts3d int2-test-helmrouts3d "${CMAKE_SOURCE_DIR}/test/Helmholtz/test_helmrouts3d.f" Helmholtz)
  fmm3d_add_fortran_driver_test(helm_hfmm3d int2-test-hfmm3d "${CMAKE_SOURCE_DIR}/test/Helmholtz/test_hfmm3d.f" Helmholtz)
  fmm3d_add_fortran_driver_test(helm_hfmm3d_scale int2-test-hfmm3d-scale "${CMAKE_SOURCE_DIR}/test/Helmholtz/test_hfmm3d_scale.f" Helmholtz)
  fmm3d_add_fortran_driver_test(helm_hfmm3d_vec int2-test-hfmm3d-vec "${CMAKE_SOURCE_DIR}/test/Helmholtz/test_hfmm3d_vec.f" Helmholtz)
  fmm3d_add_fortran_driver_test(helm_hfmm3d_mps int2-test-hfmm3d-mps "${CMAKE_SOURCE_DIR}/test/Helmholtz/test_hfmm3d_mps.f90" Helmholtz)

  fmm3d_add_fortran_driver_test(lap_laprouts3d int2-test-laprouts3d "${CMAKE_SOURCE_DIR}/test/Laplace/test_laprouts3d.f" Laplace)
  fmm3d_add_fortran_driver_test(lap_lfmm3d int2-test-lfmm3d "${CMAKE_SOURCE_DIR}/test/Laplace/test_lfmm3d.f" Laplace)
  fmm3d_add_fortran_driver_test(lap_lfmm3d_scale int2-test-lfmm3d-scale "${CMAKE_SOURCE_DIR}/test/Laplace/test_lfmm3d_scale.f" Laplace)
  fmm3d_add_fortran_driver_test(lap_lfmm3d_vec int2-test-lfmm3d-vec "${CMAKE_SOURCE_DIR}/test/Laplace/test_lfmm3d_vec.f" Laplace)

  fmm3d_add_fortran_driver_test(stokes_stokkernels int2-test-stokkernels "${CMAKE_SOURCE_DIR}/test/Stokes/test_stokkernels.f" Stokes)
  fmm3d_add_fortran_driver_test(stokes_stokkernels_rotlet_doublet int2-test-stokkernels-rotlet-doublet "${CMAKE_SOURCE_DIR}/test/Stokes/test_stokkernels_rotlet_doublet.f" Stokes)
  fmm3d_add_fortran_driver_test(stokes_stfmm3d int2-test-stfmm3d "${CMAKE_SOURCE_DIR}/test/Stokes/test_stfmm3d.f" Stokes)

  fmm3d_add_fortran_driver_test(maxwell_emfmm3d int2-test-emfmm3d "${CMAKE_SOURCE_DIR}/test/Maxwell/test_emfmm3d.f" Maxwell)

  if(WIN32)
    set(FMM3D_PYTHONPATH_SEP ";")
  else()
    set(FMM3D_PYTHONPATH_SEP ":")
  endif()
  set(FMM3D_TEST_PYTHONPATH "${CMAKE_CURRENT_BINARY_DIR}")
  if(DEFINED ENV{PYTHONPATH} AND NOT "$ENV{PYTHONPATH}" STREQUAL "")
    set(FMM3D_TEST_PYTHONPATH "${FMM3D_TEST_PYTHONPATH}${FMM3D_PYTHONPATH_SEP}$ENV{PYTHONPATH}")
  endif()

  add_test(
    NAME python_hfmm
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${FMM3D_TEST_PYTHONPATH}"
            ${Python_EXECUTABLE} -m pytest -q "${CMAKE_SOURCE_DIR}/python/test/test_hfmm.py"
  )
  add_test(
    NAME python_lfmm
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${FMM3D_TEST_PYTHONPATH}"
            ${Python_EXECUTABLE} -m pytest -q "${CMAKE_SOURCE_DIR}/python/test/test_lfmm.py"
  )
endif()
