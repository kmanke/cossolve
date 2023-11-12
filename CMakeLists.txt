cmake_minimum_required(VERSION 3.25)

# Generate compile_commands
set(CMAKE_EXPORT_COMPILE_COMMANDS true)

# Do not allow in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

project(cossolve VERSION 0.0.0)
	
# CXX standard information
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# CXX compile options
add_compile_options("-I/usr/include/eigen3")

# Build options
option(USE_MKL "Use the intel MKL library for linear algebra. If not, umfpack will be used.")

# Copy the compile_commands.json file over to the source tree so it can be used by lsp
execute_process(COMMAND cp "${CMAKE_BINARY_DIR}/compile_commands.json" "${CMAKE_SOURCE_DIR}/")
add_library(cossolve SHARED)

# If mkl is specified, use that. Otherwise, use UMFPACK.
if (USE_MKL)
   if (NOT DEFINED MKL_ROOT)
      message(STATUS "No MKL root specified. Assuming default.")
      set(MKL_ROOT "/opt/intel/oneapi/mkl/latest")
   endif()
   
   find_library(MKL_INTEL_LP64 "mkl_intel_lp64" "${MKL_ROOT}/lib/intel64")
   if (NOT MKL_INTEL_LP64)
      message(FATAL_ERROR "mkl_intel_ilp64 not found.")
   endif()

   find_library(MKL_GNU_THREAD "mkl_gnu_thread" "${MKL_ROOT}/lib/intel64")
   if (NOT MKL_GNU_THREAD)
      message(FATAL_ERROR "mkl_gnu_thread not found.")
   endif()

   find_library(MKL_CORE "mkl_core" "${MKL_ROOT}/lib/intel64")
   if (NOT MKL_CORE)
      message(FATAL_ERROR "mkl_core not found.")
   endif()

   find_library(MKL_BLACS "mkl_blacs_intelmpi_lp64" "${MKL_ROOT}/lib/intel64")
   if (NOT MKL_BLACS)
      message(FATAL_ERROR "mkl_blacs_intelmpi_lp64 not found.")
   endif()

   find_library(LIBGOMP "gomp" "/usr/lib/gcc/x86_64-pc-linux-gnu/*/")
   if (NOT LIBGOMP)
      message(FATAL_ERROR "libgomp not found.")
   endif()

   target_link_libraries(cossolve PUBLIC "${MKL_INTEL_LP64}" "${MKL_GNU_THREAD}" "${MKL_CORE}" "${LIBGOMP}")
   target_include_directories(cossolve PUBLIC "${MKL_ROOT}/include")
   add_definitions(-DUSE_MKL -DEIGEN_USE_MKL_ALL)
   add_compile_options("-fopenmp" "-lpthread" "-Wl,--no-as-needed")
else()
   find_library(UMFPACK "umfpack")
   if (NOT UMFPACK)
      message(FATAL_ERROR "umfpack not found.")
   endif()
   
   target_link_libraries(cossolve PUBLIC "${UMFPACK}")
endif()

target_include_directories(cossolve PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(cossolve PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_sources(cossolve PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src/LieGroupOps.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MatrixOps.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/capi.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/StaticSolver.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/SystemMatrix.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MatrixConstructor.cpp"
)
