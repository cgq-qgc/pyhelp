
cmake_minimum_required(VERSION 3.18)
project(pyhelp LANGUAGES Fortran)

find_package(Python COMPONENTS Interpreter REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Custom command to build HELP3O.so using f2py
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.so
    COMMAND bash -c "rm -f ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.so"
    COMMAND ${Python3_EXECUTABLE} -m numpy.f2py -c ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.FOR -m HELP3O --fcompiler=gnu95 --f90flags="-fallow-argument-mismatch -fno-range-check"
    COMMAND bash -c "mv ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O*.so ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.so"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.FOR
    COMMENT "Building HELP3O.so with f2py, cleaning old file, and renaming output"
)

add_custom_target(build_help3o ALL
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.so
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/pyhelp/HELP3O.so DESTINATION pyhelp)
