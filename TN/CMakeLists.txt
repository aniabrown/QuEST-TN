# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----- TN LIBRARY BUILD SYSTEM --------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Builds TN as a shared library libTN.so


# -----------------------------------------------------------------------------
# ----- LIBRARY USER SETTINGS -------------------------------------------------
# -----------------------------------------------------------------------------


# *****************************************************************************
# ***** NO CHANGES SHOULD BE REQUIRED FROM THE USER BEYOND THIS POINT *********
# *****************************************************************************

# -----------------------------------------------------------------------------
# ----- FIND PACKAGES ---------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# ----- SET COMPILER FLAGS ----------------------------------------------------
# -----------------------------------------------------------------------------

# ----- C COMPILER FLAGS --------------------------------------------------

# set C flags that are common between compilers and build types
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -std=c99"
)

# Use -O2 for all but debug mode by default 
if (NOT("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
        -O2"
    )
endif()

# Set c flags for release
set(CMAKE_C_FLAGS_RELEASE 
    "-O2"
)

# Set c flags to use in debug mode
set(CMAKE_C_FLAGS_DEBUG 
    "-g"
)

# TODO standardize
# set C compiler flags based on compiler type
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -Wall"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -Wall"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable cpu-dispatch"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
endif()

# ----- C++ COMPILER FLAGS --------------------------------------------------

# set C++ flags that are common between compilers and build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -std=c++98 -Wall"
)

# Use -O2 for all but debug mode by default 
if (NOT("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} \
        -O2"
    )
endif()

# Set c++ flags for release
set(CMAKE_CXX_FLAGS_RELEASE
    "-O2"
)

# Set c++ flags to use in debug mode
set(CMAKE_CXX_FLAGS_DEBUG 
    "-g"
)

# set C++ compiler flags based on compiler type
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
endif()

# -----------------------------------------------------------------------------
# ----- BUILD LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------


add_library(TN SHARED
    QuEST_tn.c
    )

# ----- Location of header files ----------------------------------------------

target_include_directories(TN
    PUBLIC .
)


# ----- Definitions -----------------------------------------------------------

# -----------------------------------------------------------------------------
# ----- LINK LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

# ----- QuEST -----------------------------------------------------------------

add_subdirectory(QuEST)

target_link_libraries(TN QuEST)



