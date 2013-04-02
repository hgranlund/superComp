PROJECT(poisson)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# comment if you only going to write fortran code
ENABLE_LANGUAGE(C)

INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-std=c99" HAVE_C99)

IF (HAVE_C99)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
ELSE(HAVE_C99)
  MESSAGE(WARNING "C99 support not available, program may fail to compile!")
ENDIF(HAVE_C99)

# needed no matter since the fst is written in fortran
ENABLE_LANGUAGE(Fortran)

# We require MPI
FIND_PACKAGE(MPI REQUIRED)

# Add the MPI include path so that our program will find mpi.h
SET(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})

SET(DEPLIBS ${DEPLIBS} ${MPI_C_LIBRARIES})

# Send the include directories to the compiler
INCLUDE_DIRECTORIES(${INCLUDES})

ADD_LIBRARY(common common.c)
ADD_EXECUTABLE(poisson oving6-poison-mpi.c fst.f)
SET_TARGET_PROPERTIES(poisson PROPERTIES LINKER_LANGUAGE C)

TARGET_LINK_LIBRARIES(poisson common ${DEPLIBS})
