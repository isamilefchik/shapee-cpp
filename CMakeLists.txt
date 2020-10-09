######################################################################
# @author      : Isa Milefchik (isavmilefchik@gmail.com)
# @file        : CMakeLists
# @created     : Saturday Jun 06, 2020 23:15:51 PDT
######################################################################

# set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required(VERSION 2.8)
project(shapee++)

enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)
include_directories(/opt/intel/mkl/include)
# include_directories(extern)

# find_package( Boost 1.72 COMPONENTS program_options REQUIRED )
# include_directories( ${Boost_INCLUDE_DIR} )
 
file(GLOB SOURCES "src/*.cpp")

# set(ENV{MKLROOT}, "/opt/intel/")
# set(ENV{BLAS_LIB_DIR}, "/opt/intel/mkl/lib")
# set(ENV{LAPACK_LIB_DIR}, "/usr/lib")
set(BLA_VENDOR Intel10_64_dyn)
find_package(LAPACK REQUIRED)
 
add_executable(shapee++ ${SOURCES})

# find_package(OpenMP REQUIRED)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
# endif()

if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

# set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/10.0.0_3/lib")
# set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/10.0.0_3/include")
set(OPENMP_INCLUDES "/usr/local/opt/libomp/include/")
find_package(OpenMP REQUIRED)

# if(OpenMP_CXX_FOUND)
#     target_link_libraries(shapee PUBLIC OpenMP::OpenMP_CXX)
# endif()

# if (OPENMP_FOUND)
#     message("FOUND OPENMP")
#     message("${OpenMP_CXX_LIB_NAMES}")
#     include_directories("${OPENMP_INCLUDES}")
#     link_directories(${OPENMP_LIBRARIES})
#     # target_link_libraries(shapee "${OpenMP_CXX_LIB_NAMES}")
#     target_link_libraries(shapee libomp libiomp5)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif(OPENMP_FOUND)

message("${OPENMP_INCLUDES}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(OpenMP_Targets "-L/usr/local/Cellar/llvm/10.0.0_3/lib/" "omp")

include_directories("${OPENMP_INCLUDES}")
target_link_libraries(shapee++ ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${OpenMP_Targets})
# target_link_libraries(shapee LINK_PUBLIC ${Boost_LIBRARIES})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
