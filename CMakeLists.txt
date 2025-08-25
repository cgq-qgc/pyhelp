cmake_minimum_required(VERSION 3.18)
project(pyhelp LANGUAGES Fortran C)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Build the extension as pyhelp.HELP3O using standard CMake
set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.FOR)

add_library(HELP3O MODULE ${SOURCES})

# Get numpy include directory using Python
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Numpy include dirs
include_directories(${NumPy_INCLUDE_DIR})

# Linker flags and Python extension properties
set_target_properties(HELP3O PROPERTIES
    PREFIX ""
    SUFFIX ".so"
    OUTPUT_NAME "HELP3O"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/pyhelp"
    LINK_FLAGS "-static -static-libgfortran -static-libgcc"
)

target_link_libraries(HELP3O PRIVATE Python::Module)

# Install to pyhelp package directory
install(TARGETS HELP3O LIBRARY DESTINATION pyhelp)
