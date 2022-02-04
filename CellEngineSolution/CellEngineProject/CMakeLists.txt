cmake_minimum_required(VERSION 3.0.0)
project(CellEngineProject)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Ofast")

include_directories(
        "../../boost/*"
        "../../Common/Platform"
	"../../Common/MathUtils"
        "../../Common/Include"
        "../../Common/Logger"
        "../../Common/FileUtils"
        "../../Common/StringUtils"
        "../../Common/DateTimeUtils"
        "../../Common/TerminalColorsUtils"
)

add_executable(
        CellEngineProject
        ../../Common/Logger/Logger.cpp
	../../Common/MathUtils/ArcBall.cpp
        ../../Common/FileUtils/FileUtils.cpp
        ../../Common/StringUtils/StringUtils.cpp
        ../../Common/DateTimeUtils/DateTimeUtils.cpp
	../../Common/MathUtils/AdditionalFunctions.cpp
        CellEnginePDBDataFile.cpp
	CellEnginePDBWindowGL.cpp
	CellEngineWindowGL.cpp
)

target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++ -lpthread -lopengl32 -lglu32)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)