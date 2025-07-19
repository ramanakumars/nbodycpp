cmake_minimum_required(VERSION 3.16)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

project(quadtree VERSION 0.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(quadtree ${SRC_FILES})

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)

include(FetchContent)
if(WIN32)
    FetchContent_Declare(SFML
        GIT_REPOSITORY https://github.com/SFML/SFML.git
        GIT_TAG 2.6.x)
    FetchContent_MakeAvailable(SFML)
    add_custom_command(
        TARGET quadtree
        COMMENT "Copy OpenAL DLL"
        PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${SFML_SOURCE_DIR}/extlibs/bin/$<IF:$<EQUAL:${CMAKE_SIZEOF_VOID_P},8>,x64,x86>/openal32.dll $<TARGET_FILE_DIR:quadtree>        
        VERBATIM)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    target_link_libraries(quadtree PUBLIC OpenMP::OpenMP_CXX)
endif()

TARGET_LINK_LIBRARIES(quadtree PRIVATE sfml-graphics)
TARGET_LINK_LIBRARIES(quadtree PRIVATE sfml-window)
TARGET_LINK_LIBRARIES(quadtree PRIVATE sfml-system)

target_include_directories(quadtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

install(TARGETS quadtree DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
