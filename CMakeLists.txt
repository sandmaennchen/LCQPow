##
##	This file is part of LCQPow.
##
##	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
##	Copyright (C) 2020 - 2021 by Jonas Hall et al.
##
##	LCQPow is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	LCQPow is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with LCQPow; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

# cmake minimum required version
cmake_minimum_required(VERSION 3.20.3)

# Project name
project(lcqpow)

## Options ------------------------------------------------------------------------------
option(
    BUILD_EXAMPLES
    "Option to build examples"
    ON
)

option(
    BUILD_DOCUMENTATION
    "Option to build documentation"
    ON
)

option(
    UNIT_TESTS
    "Option to build and run unit tests"
    ON
)

option(
    PROFILING
    "Option to enable profiling (gprof)."
    OFF
)

option(
    QPOASES_SCHUR
    "Use the Schur Complement Method with MA57 Solver"
    OFF
)

## Paths & Packages ---------------------------------------------------------------------

# 1) Install Prefix
set(
    CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}"
)

# 2) Locate Matlab and modify default values of interface & Schur
find_package(Matlab)

if (${Matlab_FOUND})
    set(
        BUILD_MATLAB_INTERFACE ON
        CACHE BOOL
        "Build matlab interface"
    )

    get_filename_component(
        MATLAB_LIBDIR
        ${Matlab_MEX_LIBRARY}
        DIRECTORY
    )

else()
    set(
        BUILD_MATLAB_INTERFACE OFF
        CACHE BOOL
        "Build matlab interface"
    )

    if (${QPOASES_SCHUR})
        message(
            FATAL_ERROR
            "qpOASES with Schur Complement method is currently only supported in combination with matlab, which can not be located."
        )
    endif()
endif()

## Print Options ------------------------------------------------------------------------

message(
    "EXAMPLES           ${BUILD_EXAMPLES}\n"
    "MATLAB INTERFACE   ${BUILD_MATLAB_INTERFACE}\n"
    "DOCUMENTATION      ${BUILD_DOCUMENTATION}\n"
    "UNIT_TESTS         ${UNIT_TESTS}\n"
    "PROFILING          ${PROFILING}\n"
    "QPOASES_SCHUR      ${QPOASES_SCHUR}\n"
)

## ADD ALL EXTERNAL PROJECTS ------------------------------------------------------------
# external modules
include(ExternalProject)

# 1) qpOASES
set(
    DEF_SOLVER
    "SOLVER_NONE"
)

if (${QPOASES_SCHUR})

    set(
        LIB_SOLVER
        "${MATLAB_LIBDIR}/libmwma57.so"
    )

    set(
        DEF_SOLVER
        "SOLVER_MA57"
    )

    set(
        USE_SOLVER
        "MA57"
    )

    set(
        LINKHSL
        "-Wl,-rpath=${MATLAB_LIBDIR}"
    )

    set(
        LAPACK_LIB
        "${MATLAB_LIBDIR}/libmwlapack.so"
    )

    set(
        BLAS_LIB
        "${MATLAB_LIBDIR}/libmwblas.so"
    )

    set(
        qpOASES_MAKE_ARGS
        ${qpOASES_MAKE_ARGS} LIB_SOLVER=${LIB_SOLVER} LIB_BLAS=${BLAS_LIB} LIB_LAPACK=${LAPACK_LIB} LINKHSL=${LINKHSL}
    )

    message(
        "Located the following dependencies for Schur Comeplement Method:\n"
        "MA57   ${LIB_SOLVER}\n"
        "BLAS   ${BLAS_LIB}\n"
        "LAPACK ${LAPACK_LIB}"
    )
endif()

set(
    qpOASES_CPP_FLAGS
    "-Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D${DEF_SOLVER} -D__USE_LONG_FINTS__ -D__NO_COPYRIGHT__"
)

set(
    qpOASES_MAKE_ARGS
    ${qpOASES_MAKE_ARGS} CPPFLAGS=${qpOASES_CPP_FLAGS}
)

ExternalProject_Add(
    qpoases
    DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/qpOASES/. ${CMAKE_BINARY_DIR}/external/src/qpoases
    PREFIX external
    CONFIGURE_COMMAND mkdir -p bin;
    BUILD_IN_SOURCE ON
    BUILD_COMMAND make ${qpOASES_MAKE_ARGS} src;
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/qpoases/bin/libqpOASES.so ${CMAKE_INSTALL_PREFIX}/lib/libqpOASES.so
)

set(qpoases_lib "-L${CMAKE_BINARY_DIR}/lib -lqpOASES")
set(qpoases_include "${CMAKE_BINARY_DIR}/external/src/qpoases/include")

# 2) OSQP
ExternalProject_Add(
    osqp
    DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/osqp/. ${CMAKE_BINARY_DIR}/external/src/osqp
    PREFIX external
    CMAKE_ARGS
        -DDLONG=OFF
        -DBUILD_SHARED_LIBS=ON
    BUILD_COMMAND cmake --build .
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/osqp-build/out/libosqp.so ${CMAKE_INSTALL_PREFIX}/lib/libosqp.so
)

set(osqp_lib "-L${CMAKE_BINARY_DIR}/lib -losqp")
set(osqp_include "${CMAKE_BINARY_DIR}/external/src/osqp/include")

# 3) googletest
ExternalProject_Add(
    gtest
    DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/googletest/. ${CMAKE_BINARY_DIR}/external/src/gtest
    PREFIX external
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/gtest-build/lib/libgtest.a ${CMAKE_INSTALL_PREFIX}/lib/
)

set(gtest_lib "-L${CMAKE_BINARY_DIR}/external/src/gtest-build/out -lgtest")
set(gtest_include "${CMAKE_BINARY_DIR}/external/src/gtest/googletest/include")

## CREATE LCQPow libraries ----------------------------------------------------------
# Include project headers
include_directories(${PROJECT_NAME} PUBLIC include)

# Add compile options
add_compile_options(
    -Wall
    -pedantic
    -Wshadow
    -Wfloat-equal
    -O3
    -Wconversion
    -Wsign-conversion
    -fPIC
    -DLINUX
    -D__USE_LONG_FINTS__
    -D${DEF_SOLVER}
)

if (${PROFILING})
    add_compile_options(-pg)
endif()

# Save auxiliar source files to variable
aux_source_directory(src SRC_FILES)

# create static lib
add_library(${PROJECT_NAME}-static STATIC ${SRC_FILES})
set_target_properties(
    ${PROJECT_NAME}-static
    PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Add the external projects as depenedncies for our project
add_dependencies(
    ${PROJECT_NAME}-static
    qpoases
    osqp
)

if (${QPOASES_SCHUR})
    target_link_libraries(
        ${PROJECT_NAME}-static
        PRIVATE ${Matlab_LIBRARIES}
    )
endif()


# Add include directories of dependencies: qpOASES, OSQP
include_directories(${PROJECT_NAME}
    PRIVATE ${qpoases_include}
    PRIVATE ${osqp_include}
)

# create shared lib
add_library(${PROJECT_NAME}-shared SHARED ${SRC_FILES})
set_target_properties(
    ${PROJECT_NAME}-shared
    PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Add the external projects as depenedncies for our project
add_dependencies(
    ${PROJECT_NAME}-shared
    qpoases
    osqp
)

if (${QPOASES_SCHUR})
    target_link_libraries(
        ${PROJECT_NAME}-shared
        PRIVATE ${Matlab_LIBRARIES}
    )
endif()

## Build examples -----------------------------------------------------------------------
if (${BUILD_EXAMPLES})
    aux_source_directory(examples EXAMPLE_FILES)

    FOREACH(ELEMENT ${EXAMPLE_FILES})
        # get filename w/o dir and extension
        get_filename_component(EXAMPLE_NAME ${ELEMENT} NAME_WE)

        # generate executable target
        add_executable(${EXAMPLE_NAME} ${ELEMENT})

        # link libraries
        target_link_libraries(
            ${EXAMPLE_NAME}
            PUBLIC ${PROJECT_NAME}-shared
            PRIVATE ${qpoases_lib} ${osqp_lib}
        )

        if (${QPOASES_SCHUR})
            target_link_libraries(
                ${EXAMPLE_NAME}
                PRIVATE ${Matlab_LIBRARIES}
            )
        endif()

        # specify output directory
        set_target_properties(
            ${EXAMPLE_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples"
        )
    endforeach()
endif()

## Build Matlab interface ---------------------------------------------------------------

if (${BUILD_MATLAB_INTERFACE})

    if (${QPOASES_SCHUR})
        set(
            MA57_lib
            -lmwma57
        )
    endif()

    set(
        MATLAB_INTERFACE_DIR
        ${CMAKE_SOURCE_DIR}/interfaces/matlab
    )

    set(
        MATLAB_INTERFACE_NAME
        LCQPow
    )

    set(
        MATLAB_INTERFACE_DESTINATION
        ${CMAKE_BINARY_DIR}/lib
    )


    matlab_add_mex(
        NAME matlab_interface
        SRC ${MATLAB_INTERFACE_DIR}/${MATLAB_INTERFACE_NAME}.cpp
        OUTPUT_NAME ${MATLAB_INTERFACE_DESTINATION}/${MATLAB_INTERFACE_NAME}
        LINK_TO "-llcqpow -lqpOASES -losqp -lmwblas -lmwlapack -lmwma57"
    )

    target_link_directories(
        matlab_interface
        PUBLIC ${CMAKE_BINARY_DIR}/lib
    )

    target_compile_options(
        matlab_interface
        PUBLIC -D${DEF_SOLVER}
        PUBLIC -D__USE_LONG_FINTS__
        PUBLIC -D__MATLAB__
        PUBLIC -D__NO_COPYRIGHT__
        PUBLIC -D__cpluplus
        PUBLIC -O
        PUBLIC -largeArrayDims
        PUBLIC -lmwblas
        PUBLIC -lmwlapack
        PUBLIC -lmwma57
        PUBLIC -llcqpow
        PUBLIC -lqpOASES
        PUBLIC -losqp
    )

    target_include_directories(
        matlab_interface
        PUBLIC include
        PUBLIC ${osqp_include}
        PUBLIC ${qpoases_include}
    )

    file(
        COPY ${MATLAB_INTERFACE_DIR}/${MATLAB_INTERFACE_NAME}.m
        DESTINATION ${MATLAB_INTERFACE_DESTINATION}
    )

endif()


## Build unit tests ---------------------------------------------------------------------
if (${UNIT_TESTS})
    # Add include directories of dependencies: qpOASES, OSQP
    include_directories(${PROJECT_NAME}
        PRIVATE ${gtest_include}
        PRIVATE ${qpoases_include}
        PRIVATE ${osqp_include}
    )

    # Add Unit testing source file
    add_executable(RunUnitTests test/RunUnitTests.cpp)

    add_dependencies(
        RunUnitTests
        qpoases
        osqp
        gtest
    )

    set_target_properties(
        RunUnitTests
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
    )

    message(
        INFO
        ${qpoases_lib} ${osqp_lib} ${gtest_lib}
    )

    # Add libraries: gtest, qpOASES, OSQP
    target_link_libraries(
        RunUnitTests
        PUBLIC ${PROJECT_NAME}-shared
        PRIVATE ${qpoases_lib} ${osqp_lib} ${gtest_lib}
        PRIVATE -lpthread
    )
endif()

# ## Build documentation
# if (${BUILD_DOCUMENTATION})
#     find_package(
#         Doxygen
#         REQUIRED dot
#     )

#     if (DOXYGEN_FOUND)
#         message(
#             STATUS
#             "Found doxygen: Building documentation."
#         )

#         doxygen_add_docs(
#             doxy
#             "include"
#         )

#         set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../docs")

#         set(DOXYGEN_GENERATE_HTML YES)
#         set(DOXYGEN_GENERATE_MAN YES)

#         set(DOXYGEN_PROJECT_BRIEF "LCQPow -- A Linear Complementarity Quadratic Programming Solver")
#     else()
#         message(
#             WARNING
#             "Could not find doxygen: Not building documentation."
#         )
#     endif()
# endif()

## End of cmake file --------------------------------------------------------------------
