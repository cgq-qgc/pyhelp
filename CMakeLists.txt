
cmake_minimum_required(VERSION 3.18)
project(pyhelp LANGUAGES Fortran)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Use f2py to build the HELP3O Python extension from Fortran
execute_process(
    COMMAND ${Python_EXECUTABLE} -m numpy.f2py -c -m HELP3O ${CMAKE_SOURCE_DIR}/pyhelp/HELP3O.FOR
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Install the built module to the pyhelp package directory
install(FILES ${CMAKE_BINARY_DIR}/HELP3O*.so DESTINATION pyhelp)
