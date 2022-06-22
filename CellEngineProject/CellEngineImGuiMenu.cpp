
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <sb7.h>
#include <sb7color.h>
#include <vmath.h>
#include <object.h>
#include <shader.h>
#include <sb7textoverlay.h>

#include <string>
#include <memory>

#include "Logger.h"
#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineDataFile.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineConfigurationFileReaderWriter.h"

using namespace std;

struct APPINFO
{
    char Title[128];
    int WindowWidth;
    int WindowHeight;
    int MajorVersion;
    int MinorVersion;
    int Samples;
    union
    {
        struct
        {
            unsigned int FullScreen : 1;
            unsigned int VSync : 1;
            unsigned int Cursor : 1;
            unsigned int Stereo : 1;
            unsigned int Debug : 1;
            unsigned int Robust : 1;
        };
        unsigned int All;
    }
    Flags;
};

APPINFO Info;

void InitializeLoggerManagerParameters()
{
    try
    {
        using namespace string_utils;

        LoggersManagerObject.InitializeFilesNames({ "AllMessages" });
        LoggersManagerObject.InitializeSelectiveWordsFunctions({ [](const string& s) { return true; } });
        LoggersManagerObject.InitializeLoggerManagerDataForTask("CELL_RESULTS", ".\\", string("Logs." + GetActualDateTimeStandardCPP(".", ".", ".", ".", ".")), true, 0, function<void(const uint64_t& CurrentThreadId, const uint64_t FileNumber, const string& MessageStr)>());
        LoggersManagerObject.InitializePrintingParameters(CellEngineConfigurationFileReaderWriterObject.PrintLogToConsole,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogToFiles,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogLineNumberToConsole,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogDateTimeToConsole,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogProcessIdToConsole,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogProcessPriorityLevelToConsole,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogThreadIdToConsole,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogLineNumberToFile,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogDateTimeToFile,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogProcessIdToFile,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogProcessPriorityLevelToFile,
                                                          CellEngineConfigurationFileReaderWriterObject.PrintLogThreadIdToFile,
                                                          CellEngineConfigurationFileReaderWriterObject.MaximalNumberOfLinesInOneFile);
    }
    CATCH("initializing logger manager parameters")
}

void ReadInitConfiguration()
{
    try
    {
        InitializeLoggerManagerParameters();
        LoggersManagerObject.Log(STREAM("START CELL"));

        uint64_t ExecuteCellStateId;
        if (__argc > 1)
            ExecuteCellStateId = stoi(__argv[1]);
        else
            LoggersManagerObject.Log(STREAM("Lack of cell id to execute in program parameters"));

        CellEngineConfigurationFileReaderWriterObject.ReadChessConfigurationFile("CellEngineProjectConfig.xml", ExecuteCellStateId);
    }
    CATCH("reading of data file")
}

unique_ptr<CellEngineOpenGLVisualiser> CellEngineOpenGLVisualiserPointer;

void CellEngineOpenGLVisualiserThreadFunction(int XPosWindow, int YPosWindow, int WidthWindow, int HeightWindow)
{
    CellEngineOpenGLVisualiserPointer = make_unique<CellEngineOpenGLVisualiser>();
    CellEngineOpenGLVisualiserPointer->Run(XPosWindow, YPosWindow, WidthWindow, HeightWindow);
}

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

int main(int argc, const char ** argv)
{
    ReadInitConfiguration();

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    Info.WindowWidth = 1600;
    Info.WindowHeight = 1200;
    Info.MajorVersion = 4;
    Info.MinorVersion = 3;
    Info.Samples = 0;
    Info.Flags.All = 0;
    Info.Flags.Cursor = 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, Info.MajorVersion);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, Info.MinorVersion);

    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

    if (Info.Flags.Robust)
        glfwWindowHint(GLFW_CONTEXT_ROBUSTNESS, GLFW_LOSE_CONTEXT_ON_RESET);

    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, Info.Samples);
    glfwWindowHint(GLFW_STEREO, Info.Flags.Stereo ? GL_TRUE : GL_FALSE);

    GLFWwindow* window = glfwCreateWindow(CellEngineConfigurationFileReaderWriterObject.WidthMenuWindow, CellEngineConfigurationFileReaderWriterObject.HeightMenuWindow, "Dear ImGui GLFW+OpenGL3 example", NULL, NULL);
    if (window == NULL)
        return 1;

    glfwSetWindowPos(window, CellEngineConfigurationFileReaderWriterObject.XTopMenuWindow, CellEngineConfigurationFileReaderWriterObject.YTopMenuWindow);

    if (!Info.Flags.Cursor)
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    gl3wInit();

    glfwMakeContextCurrent(window);

    glfwSwapInterval(1); // Enable vsync

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    const char* glsl_version = "#version 130";
    ImGui_ImplOpenGL3_Init(glsl_version);

    bool show_demo_window = true;
    bool show_another_window = true;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    thread CellEngineOpenGLVisualiserThreadObject(CellEngineOpenGLVisualiserThreadFunction, CellEngineConfigurationFileReaderWriterObject.XTopMainWindow, CellEngineConfigurationFileReaderWriterObject.YTopMainWindow, CellEngineConfigurationFileReaderWriterObject.WidthMainWindow, CellEngineConfigurationFileReaderWriterObject.HeightMainWindow);

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        static float f = 0.0f;
        static int counter = 0;

        ImGui::Begin("Hello, world!");

        ImGui::Text("This is some useful text.");
        ImGui::Checkbox("Demo Window", &show_demo_window);
        ImGui::Checkbox("Another Window", &show_another_window);

        ImGui::SliderFloat("float", &f, 0.0f, 1.0f);
        ImGui::ColorEdit3("clear color", (float*)&clear_color);

        if (ImGui::Button("Button"))
        {
            CellEngineOpenGLVisualiserPointer->OnMouseWheel(1);
            counter++;
        }
        ImGui::SameLine();
        ImGui::Text("counter = %d", counter);

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        ImGui::End();

        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }

        ImGui::Render();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    CellEngineOpenGLVisualiserThreadObject.detach();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
