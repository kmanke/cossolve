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

# Here we define all of the modules that should be built and linked
#if(true)
#  add_subdirectory(modules/messaging)
#endif()

# Copy the compile_commands.json file over to the source tree so it can be used by lsp
execute_process(COMMAND cp "${CMAKE_BINARY_DIR}/compile_commands.json" "${CMAKE_SOURCE_DIR}/")
add_library(cossolve SHARED)

# Find cholmod
find_library(MKL "mkl")
find_library(UMFPACK "umfpack")
target_link_libraries(cossolve PUBLIC "${UMFPACK}")

target_include_directories(cossolve PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(cossolve PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_sources(cossolve PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src/se3.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/capi.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Solver.cpp"
)
