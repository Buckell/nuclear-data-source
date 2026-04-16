include(FetchContent)

cmake_minimum_required(VERSION 3.28)
project(nuclear_data_source)

set(CMAKE_CXX_STANDARD 20)

add_executable(nuclear_data_source main.cpp)

set(Python_USE_STATIC_LIBS ON)
find_package(Python3 COMPONENTS Development REQUIRED)

target_link_libraries(nuclear_data_source PRIVATE Python3::Python)

FetchContent_Declare(
    ZipLib
    GIT_REPOSITORY https://github.com/DreamyCecil/ZipLib.git
)
FetchContent_MakeAvailable(ZipLib)
include_directories(${ZipLib_SOURCE_DIR}/Source/ZipLib)

add_library(ZipLib::ZipLib UNKNOWN IMPORTED)

set_target_properties(ZipLib::ZipLib PROPERTIES
    IMPORTED_LOCATION "${ZipLib_BINARY_DIR}/Source/ZipLib/ZipLib.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${ZipLib_SOURCE_DIR}/Source/ZipLib"
)

target_link_libraries(nuclear_data_source PRIVATE ZipLib::ZipLib)

add_library(ZipLib::bzip2 UNKNOWN IMPORTED)

set_target_properties(ZipLib::bzip2 PROPERTIES
    IMPORTED_LOCATION "${ZipLib_BINARY_DIR}/Source/ZipLib/extlibs/bzip2/bzip2.lib"
)

target_link_libraries(nuclear_data_source PRIVATE ZipLib::bzip2)

add_library(ZipLib::lzma UNKNOWN IMPORTED)

set_target_properties(ZipLib::lzma PROPERTIES
    IMPORTED_LOCATION "${ZipLib_BINARY_DIR}/Source/ZipLib/extlibs/lzma/lzma.lib"
)

target_link_libraries(nuclear_data_source PRIVATE ZipLib::lzma)

add_library(ZipLib::zlib UNKNOWN IMPORTED)

set_target_properties(ZipLib::zlib PROPERTIES
        IMPORTED_LOCATION "${ZipLib_BINARY_DIR}/Source/ZipLib/extlibs/zlib/zlib.lib"
)

target_link_libraries(nuclear_data_source PRIVATE ZipLib::zlib)
