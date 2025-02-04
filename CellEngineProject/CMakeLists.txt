cmake_minimum_required(VERSION 3.0.0)

set(USING_MODULES "FALSE")

project(CellEngineProject)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Ofast")
set(CMAKE_CXX_FLAGS "-fmodules-ts")

LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL)

link_directories("${OPENMP_LIBRARIES}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

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
	CellEngineMacros.h
	CellEngineUseful.h
	CellEngineConstants.h
	CellEngineColors.h
	CellEngineColors.cpp
	CellEngineAminoAcids.h
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
	CellEngineMain.cpp
	CellEngineImGuiMenu.h
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
	CellEngineExecutionTimeStatistics.h
	CellEngineBasicVoxelsOperations.h
	CellEngineNucleicAcidsBasicOperations.h
	CellEngineParticlesFullAtomOperations.h
	CellEngineParticlesVoxelsOperations.h
	CellEngineBasicParallelExecutionData.h
	CellEngineBasicParticlesOperations.cpp
	CellEngineBasicParticlesOperations.h
	CellEngineParticlesVoxelsShapesGenerator.h
	CellEngineParticlesVoxelsShapesGenerator.cpp
	CellEngineParticlesFullAtomShapesGenerator.cpp
	CellEngineParticlesFullAtomShapesGenerator.h
	CellEngineRealRandomParticlesInVoxelSpaceGenerator.cpp
	CellEngineRealRandomParticlesInVoxelSpaceGenerator.h
	CellEngineRealRandomParticlesInFullAtomSpaceGenerator.cpp
	CellEngineRealRandomParticlesInFullAtomSpaceGenerator.h
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
	CellEngineFullAtomSimulationSpace.cpp
	CellEngineFullAtomSimulationSpace.h
	CellEngineNucleicAcidsComplexOperations.cpp
	CellEngineNucleicAcidsComplexOperations.h
	CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.cpp
	CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.h
	CellEngineChemicalReactionsInVoxelSimulationSpace.cpp
	CellEngineChemicalReactionsInVoxelSimulationSpace.h
	CellEngineChemicalReactionsInFullAtomSimulationSpace.cpp
	CellEngineChemicalReactionsInFullAtomSimulationSpace.h
	CellEngineOpenGLVisualiser.h
	CellEngineOpenGLVisualiser.cpp
	CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h
	CellEngineOpenGLVisualiserOfVoxelSimulationSpace.cpp
	CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h
	CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.cpp
	CellEngineConfigurationFileReaderWriter.h
	CellEngineConfigurationFileReaderWriter.cpp
	CellEngineGenesPromotersAndStartCodonsFinder.cpp
	CellEngineGenesPromotersAndGenesStartCodonsFinder.h
	CellEngineSimulationParallelExecutionManager.cpp
	CellEngineSimulationParallelExecutionManager.h
	CellEngineWellStirredChemicalReactionsSimulation.h
	CellEngineWellStirredChemicalReactionsSimulation.h
	CellEngineWellStirredChemicalReactionsSimulation.cpp)

if(USING_MODULES STREQUAL "TRUE")
target_sources(CellEngineProject
		PUBLIC
		FILE_SET all_my_modules TYPE CXX_MODULES FILES
		CellEngineColors.cpp)
endif()

if(WIN32)
target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread -lglew32 -lsb7 -lglfw3 -lopengl32 -lglu32 -limgui)
else()
target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread sb7 -fopenmp ${OPENGL_LIBRARIES} glfw3 X11 Xrandr Xinerama Xi Xxf86vm Xcursor GL rt dl imgui)
endif()

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)
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