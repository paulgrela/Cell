cmake_minimum_required(VERSION 3.0.0)

project(CellEngineProject)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 23)
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

set(COMMON_LIBS ${COMMON_LIBS} ${EXTRA_LIBS})

add_library(
	sb7
	srcsb7/sb7.cpp
	srcsb7/sb7color.cpp
	srcsb7/sb7ktx.cpp
	srcsb7/sb7object.cpp
	srcsb7/sb7shader.cpp
	srcsb7/sb7textoverlay.cpp
	srcsb7/gl3w.c)

add_library(
	imgui
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
	"../Common/DoublyLinkedList"
	"../Common/TerminalColorsUtils")

add_executable(
	CellEngineProject
	../Common/Logger/Logger.h
	../Common/Logger/Logger.cpp
	../Common/MathUtils/ArcBall.h
	../Common/MathUtils/ArcBall.cpp
	../Common/MathUtils/Combinatorics.cpp
	../Common/MathUtils/Combinatorics.h
	../Common/FileUtils/FileUtils.h
	../Common/FileUtils/FileUtils.cpp
	../Common/StringUtils/StringUtils.h
	../Common/StringUtils/StringUtils.cpp
	../Common/DateTimeUtils/DateTimeUtils.h
	../Common/DateTimeUtils/DateTimeUtils.cpp
	../Common/DoublyLinkedList/DoublyLinkedList.h
	../Common/DoublyLinkedList/DoublyLinkedList.cpp
	CellEngineTypes.h
	CellEngineUseful.h
	CellEngineConstants.h
	CellEngineColors.h
	CellEngineColors.cpp
	CellEngineAtom.h
	CellEngineAtom.cpp
	CellEngineParticle.h
	CellEngineParticleKind.h
	CellEngineParticlesKindsManager.h
	CellEngineChemicalReaction.h
	CellEngineChemicalReactionsManager.h
	CellEngineChemicalReactionsEngine.h
	CellEngineChemicalReactionsEngine.cpp
	CellEngineRandomDeviceEngine.cpp
	CellEngineRandomDeviceEngine.h
	CellEngineDataFile.h
	CellEngineFilmOfStructures.h
	CellEngineImGuiMenu.cpp
	CellEnginePDBDataFileReader.h
	CellEnginePDBDataFileReader.cpp
	CellEngineCIFDataFileReader.h
	CellEngineCIFDataFileReader.cpp
	CellEngineParticlesDataFile.cpp
	CellEngineParticlesDataFile.h
	CellEngineBuildParticlesDataOperations.h
	CellEngineIllinoisDataCreator.cpp
	CellEngineIllinoisDataCreator.h
	CellEngineCompiledDataCreator.cpp
	CellEngineCompiledDataCreator.h
	CellEngineBasicVoxelsOperations.h
	CellEngineNucleicAcidsBasicOperations.h
	CellEngineParticlesVoxelsOperations.h
	CellEngineBasicParticlesOperations.cpp
	CellEngineBasicParticlesOperations.h
	CellEngineParticlesVoxelsShapesGenerator.h
	CellEngineParticlesVoxelsShapesGenerator.cpp
	CellEngineRealRandomParticlesInVoxelSpaceGenerator.cpp
	CellEngineRealRandomParticlesInVoxelSpaceGenerator.h
	CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.cpp
	CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h
	CellEngineTestParticlesInVoxelSpaceGenerator.cpp
	CellEngineTestParticlesInVoxelSpaceGenerator.h
	CellEngineChemicalReactionsInBasicSimulationSpace.cpp
	CellEngineChemicalReactionsInBasicSimulationSpace.h
	CellEngineChemicalReactionsInSimulationSpace.cpp
	CellEngineChemicalReactionsInSimulationSpace.h
	CellEngineVoxelSimulationSpaceStatistics.cpp
	CellEngineVoxelSimulationSpaceStatistics.h
	CellEngineSimulationSpaceStatistics.cpp
	CellEngineSimulationSpaceStatistics.h
	CellEngineParticlesBinaryDataFileReaderWriter.cpp
	CellEngineParticlesBinaryDataFileReaderWriter.h
	CellEngineDataBuilderForVoxelSimulationSpace.h
	CellEngineDataBuilderForFullAtomSimulationSpace.h
	CellEngineSimulationSpace.h
	CellEngineSimulationSpace.cpp
	CellEngineVoxelSimulationSpace.cpp
	CellEngineVoxelSimulationSpace.h
	CellEngineNucleicAcidsComplexOperations.cpp
	CellEngineNucleicAcidsComplexOperations.h
	CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.cpp
	CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.h
	CellEngineChemicalReactionsInVoxelSimulationSpace.cpp
	CellEngineChemicalReactionsInVoxelSimulationSpace.h
	CellEngineOpenGLVisualiser.h
	CellEngineOpenGLVisualiser.cpp
	CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h
	CellEngineOpenGLVisualiserOfVoxelSimulationSpace.cpp
	CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h
	CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.cpp
	CellEngineConfigurationFileReaderWriter.h
	CellEngineConfigurationFileReaderWriter.cpp
	CellEngineWellStirredChemicalReactionsSimulation.h
	CellEngineWellStirredChemicalReactionsSimulation.h
	CellEngineWellStirredChemicalReactionsSimulation.cpp)

if(WIN32)
target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread -lglew32 -lsb7 -lglfw3 -lopengl32 -lglu32 -limgui)
else()
target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread sb7 ${OPENGL_LIBRARIES} glfw3 X11 Xrandr Xinerama Xi Xxf86vm Xcursor GL rt dl imgui)
endif()

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include(CTest)
enable_testing()

set (GMOCK_DIR "/usr/src/googletest/googlemock")
add_subdirectory(${GMOCK_DIR} ${CMAKE_BINARY_DIR}/gmock)
include_directories("/usr/src/googletest/googlemock/include")

add_executable(
	CellEngineRunTests
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
	../Common/DoublyLinkedList/DoublyLinkedList.h
	../Common/DoublyLinkedList/DoublyLinkedList.cpp
	CellEngineTypes.h
	CellEngineUseful.h
	CellEngineColors.h
	CellEngineColors.cpp
	CellEngineAtom.h
	CellEngineAtom.cpp
	CellEngineParticle.h
	CellEngineChemicalReaction.h
	CellEngineChemicalReactionsEngine.h
	CellEngineChemicalReactionsEngine.cpp
	CellEngineWellStirredChemicalReactionsSimulation.h
	CellEngineWellStirredChemicalReactionsSimulation.cpp
	tests/CellEngineRunTests.cpp)

if(WIN32)
target_link_libraries(CellEngineRunTests -static-libgcc -static-libstdc++ -lpthread -lglew32 -lsb7 -lglfw3 -lopengl32 -lglu32 -limgui)
else()
target_link_libraries(CellEngineRunTests -static-libgcc -static-libstdc++ -lpthread sb7 ${OPENGL_LIBRARIES} glfw3 X11 Xrandr Xinerama Xi Xxf86vm Xcursor GL rt dl imgui gtest gmock)
endif()

add_test(CellEngineRunTests CellEngineRunTests)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)