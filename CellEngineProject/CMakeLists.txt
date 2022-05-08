cmake_minimum_required(VERSION 3.0.0)

project(CellEngineProject)

#set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Ofast")

LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

find_package(OpenGL)

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#set (CMAKE_DEBUG_POSTFIX "_d")
#
#if(WIN32)
#	set(COMMON_LIBS sb7 optimized glfw3 debug glfw3_d ${GLFW_LIBRARIES} ${OPENGL_LIBRARIES})
#elseif (UNIX)
#	find_package(PkgConfig REQUIRED)
#	pkg_check_modules(GLFW REQUIRED glfw3)
#	set(COMMON_LIBS sb7 glfw3 X11 Xrandr Xinerama Xi Xxf86vm Xcursor GL rt dl)
#else()
#	set(COMMON_LIBS sb7)
#endif()
#
#set(COMMON_LIBS ${COMMON_LIBS} ${EXTRA_LIBS})

add_library(sb7
		srcsb7/sb7.cpp
		srcsb7/sb7color.cpp
		srcsb7/sb7ktx.cpp
		srcsb7/sb7object.cpp
		srcsb7/sb7shader.cpp
		srcsb7/sb7textoverlay.cpp
		srcsb7/gl3w.c)

include_directories(include)
include_directories(extern/glfw-3.0.4/include)
include_directories(
	"../../Boost/*"
	"../Common/Platform"
	"../Common/MathUtils"
	"../Common/Include"
	"../Common/Logger"
	"../Common/FileUtils"
	"../Common/StringUtils"
	"../Common/DateTimeUtils"
	"../Common/TerminalColorsUtils")

add_executable(
	CellEngineProject
	../Common/Logger/Logger.cpp
	../Common/MathUtils/ArcBall.cpp
	../Common/FileUtils/FileUtils.cpp
	../Common/StringUtils/StringUtils.cpp
	../Common/DateTimeUtils/DateTimeUtils.cpp
	../Common/MathUtils/AdditionalFunctions.cpp
	CellEnginePDBDataFile.cpp
	CellEngineCIFDataFile.cpp
	CellEngineOpenGLVisualizer.cpp)

target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread -lglew32 -lsb7 -lglfw3 -lopengl32 -lglu32)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)