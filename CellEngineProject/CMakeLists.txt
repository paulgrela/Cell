cmake_minimum_required(VERSION 3.0.0)

project(CellEngineProject)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Ofast")
set(CMAKE_CXX_FLAGS "-fmodules-ts")

LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

find_package(OpenGL)

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(COMMON_LIBS ${COMMON_LIBS} ${EXTRA_LIBS})

add_library(sb7
		srcsb7/sb7.cpp
		srcsb7/sb7color.cpp
		srcsb7/sb7ktx.cpp
		srcsb7/sb7object.cpp
		srcsb7/sb7shader.cpp
		srcsb7/sb7textoverlay.cpp
		srcsb7/gl3w.c)

add_library(imgui
		extern/imgui/imgui.cpp
		extern/imgui/imgui_demo.cpp
		extern/imgui/imgui_draw.cpp
		extern/imgui/imgui_tables.cpp
		extern/imgui/imgui_widgets.cpp
		extern/imgui/backends/imgui_impl_glfw.cpp
		extern/imgui/backends/imgui_impl_opengl3.cpp)

install(TARGETS sb7 DESTINATION lib/sb7)

include_directories(glm)
include_directories(include)

include_directories(extern/imgui)
include_directories(extern/imgui/backends)

include_directories(extern/glfw-3.3.7.bin.WIN64/include)
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
	../Common/Logger/Logger.h
	../Common/Logger/Logger.cpp
	../Common/MathUtils/ArcBall.h
	../Common/MathUtils/ArcBall.cpp
	../Common/FileUtils/FileUtils.h
	../Common/FileUtils/FileUtils.cpp
	../Common/StringUtils/StringUtils.h
	../Common/StringUtils/StringUtils.cpp
	../Common/DateTimeUtils/DateTimeUtils.h
	../Common/DateTimeUtils/DateTimeUtils.cpp
	CellEngineTypes.h
	CellEngineColors.h
	CellEngineColors.cpp
	CellEngineAtom.h
	CellEngineAtom.cpp
	CellEngineParticle.h
	CellEngineReaction.h
	CellEngineImGuiMenu.cpp
	CellEnginePDBDataFile.h
	CellEnginePDBDataFile.cpp
	CellEngineCIFDataFile.h
	CellEngineCIFDataFile.cpp
	CellEngineCIFDataFileReaderOfVoxelSimulationSpace.h
	CellEngineCIFDataFileReaderOfFullAtomSimulationSpace.h
	CellEngineVoxelSimulationSpace.h
	CellEngineVoxelSimulationSpace.cpp
	CellEngineOpenGLVisualiser.h
	CellEngineOpenGLVisualiser.cpp
	CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h
	CellEngineOpenGLVisualiserOfVoxelSimulationSpace.cpp
	CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h
	CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.cpp
	CellEngineConfigurationFileReaderWriter.h
	CellEngineConfigurationFileReaderWriter.cpp
	CellEngineWellStirredChemicalReactionsSimulation.h
	CellEngineWellStirredChemicalReactionsSimulation.cpp)

if(WIN32)
target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread -lglew32 -lsb7 -lglfw3 -lopengl32 -lglu32 -limgui)
else()
target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread sb7 ${OPENGL_LIBRARIES} glfw3 X11 Xrandr Xinerama Xi Xxf86vm Xcursor GL rt dl imgui)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)