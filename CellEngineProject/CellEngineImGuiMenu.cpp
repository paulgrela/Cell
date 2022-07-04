
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <string>
#include <memory>

#include "Logger.h"
#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
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
        LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);
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

void DrawPlusMinusScalarButton(float& VariableToChange, const int64_t Step, const int64_t MinValue, const int64_t MaxValue, const string& Description, int& IDButton)
{
    ImGui::PushID(IDButton);
    ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0 / 7.0f, 0.6f, 0.6f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0 / 7.0f, 0.7f, 0.7f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0 / 7.0f, 0.8f, 0.8f));
    if (ImGui::Button(" - "))
        if (VariableToChange - Step >= MinValue)
            VariableToChange -= Step;
    ImGui::PopStyleColor(3);
    ImGui::PopID();
    IDButton++;

    ImGui::SameLine();
    ImGui::PushID(IDButton);
    ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(3 / 7.0f, 0.6f, 0.6f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(3 / 7.0f, 0.7f, 0.7f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(3 / 7.0f, 0.8f, 0.8f));
    if (ImGui::Button(" + "))
        if (VariableToChange + Step <= MaxValue)
            VariableToChange += Step;
    ImGui::PopStyleColor(3);
    ImGui::PopID();
    IDButton++;

    ImGui::SameLine();
    ImGui::Text("%s", string(to_string(VariableToChange) + " [" + Description + "]").c_str());
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

    GLFWwindow* window = glfwCreateWindow(CellEngineConfigDataObject.WidthMenuWindow, CellEngineConfigDataObject.HeightMenuWindow, "Dear ImGui GLFW+OpenGL3 example", NULL, NULL);
    if (window == NULL)
        return 1;

    glfwSetWindowPos(window, CellEngineConfigDataObject.XTopMenuWindow, CellEngineConfigDataObject.YTopMenuWindow);

    if (!Info.Flags.Cursor)
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    gl3wInit();

    glfwMakeContextCurrent(window);

    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    //ImGui::StyleColorsLight();
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    const char* glsl_version = "#version 130";
    ImGui_ImplOpenGL3_Init(glsl_version);

    bool show_demo_window = true;
    bool show_another_window = true;
    ImVec4 BackgroundColor = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    thread CellEngineOpenGLVisualiserThreadObject(CellEngineOpenGLVisualiserThreadFunction, CellEngineConfigDataObject.XTopMainWindow, CellEngineConfigDataObject.YTopMainWindow, CellEngineConfigDataObject.WidthMainWindow, CellEngineConfigDataObject.HeightMainWindow);

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        //static bool closef = false;
        //bool* p_open = &closef;

//        ImGuiWindowFlags window_flags = 0;
//        window_flags |= ImGuiWindowFlags_NoMove;
//        window_flags |= ImGuiWindowFlags_NoResize;
//        ImGui::Begin("Cell Engine Visualiser", NULL, window_flags);
        ImGui::Begin("Cell Engine Visualiser");


//        static int e = 0;
//        ImGui::RadioButton("radio a", &e, 0); ImGui::SameLine();
//        ImGui::RadioButton("radio b", &e, 1); ImGui::SameLine();
//        ImGui::RadioButton("radio c", &e, 2);



        ImGui::Text("STATUS");

        int IDButton = 1;

        if (ImGui::CollapsingHeader("View Move"))
        {
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewXMoveShortStep, 1, 1, 10, "View X Move Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewYMoveShortStep, 1, 1, 10, "View Y Move Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewZMoveShortStep, 1, 1, 10, "View Z Move Short Step", IDButton);

            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewX, CellEngineConfigDataObject.ViewXMoveShortStep, -3000, 3000, "View X Change Using Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewY, CellEngineConfigDataObject.ViewYMoveShortStep, -3000, 3000, "View Y Change Using Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewZ, CellEngineConfigDataObject.ViewZMoveShortStep, -3000, 3000, "View Z Change Using Short Step", IDButton);

            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewXMoveLongStep, 10, 0, 100, "View X Move Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewYMoveLongStep, 10, 0, 100, "View Y Move Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewZMoveLongStep, 10, 0, 100, "View Z Move Long Step", IDButton);

            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewX, CellEngineConfigDataObject.ViewXMoveLongStep, -3000, 3000, "View X Change Using Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewY, CellEngineConfigDataObject.ViewYMoveLongStep, -3000, 3000, "View Y Change Using Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewZ, CellEngineConfigDataObject.ViewZMoveLongStep, -3000, 3000, "View Z Change Using Long Step", IDButton);
        }
        if (ImGui::CollapsingHeader("Camera Move"))
        {
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraXMoveShortStep, 1, 1, 10, "Camera X Move Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraYMoveShortStep, 1, 1, 10, "Camera Y Move Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraZMoveShortStep, 1, 1, 10, "Camera Z Move Short Step", IDButton);

            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraXPosition, CellEngineConfigDataObject.CameraXMoveShortStep, -3000, 3000, "Camera X Position Change Using Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraYPosition, CellEngineConfigDataObject.CameraYMoveShortStep, -3000, 3000, "Camera Y Position Change Using Short Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraZPosition, CellEngineConfigDataObject.CameraZMoveShortStep, -3000, 3000, "Camera Z Position Change Using Short Step", IDButton);

            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraXMoveLongStep, 10, 0, 100, "Camera X Move Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraYMoveLongStep, 10, 0, 100, "Camera Y Move Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraZMoveLongStep, 10, 0, 100, "Camera Z Move Long Step", IDButton);

            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraXPosition, CellEngineConfigDataObject.CameraXMoveLongStep, -3000, 3000, "Camera X Position Change Using Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraYPosition, CellEngineConfigDataObject.CameraYMoveLongStep, -3000, 3000, "Camera Y Position Change Using Long Step", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.CameraZPosition, CellEngineConfigDataObject.CameraZMoveLongStep, -3000, 3000, "Camera Z Position Change Using Long Step", IDButton);
        }
        if (ImGui::CollapsingHeader("Rotation Angle"))
        {
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.RotationAngle1, 1, -360, 360, "Rotation Angle 1", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.RotationAngle2, 1, -360, 360, "Rotation Angle 2", IDButton);
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.RotationAngle3, 1, -360, 360, "Rotation Angle 3", IDButton);
        }


        ImGui::Checkbox("Log parameters of rendering to file", &CellEngineConfigDataObject.LogParametersOfRenderingToFile);

        bool UseStencilBuffer = true;
        CellEngineConfigDataObject.NumberOfStencilBufferLoops == 1 ? UseStencilBuffer = false : UseStencilBuffer = true;
        ImGui::Checkbox("Show details of picked atom", &UseStencilBuffer);
        UseStencilBuffer == true ? CellEngineConfigDataObject.NumberOfStencilBufferLoops = 3 : CellEngineConfigDataObject.NumberOfStencilBufferLoops = 1;

        ImGui::Checkbox("Print atom description on screen", &CellEngineConfigDataObject.PrintAtomDescriptionOnScreen);

        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops == 3)
        {
            ImGui::Text("ATOM DATA:");
            ImGui::Text("%s", string(CellEngineConfigDataObject.AtomDescriptionStr1 + " " + CellEngineConfigDataObject.AtomDescriptionStr2).c_str());
            ImGui::Text("%s", string(CellEngineConfigDataObject.AtomDescriptionStr3 + " " + CellEngineConfigDataObject.AtomDescriptionStr4).c_str());
        }

        ImGui::Text("%s", CellEngineConfigDataObject.TimeParametersOfRenderingStr.c_str());
        ImGui::Text("%s", CellEngineConfigDataObject.NumberOfRenderedAtomsParametersOfRenderingStr.c_str());

        ImGui::ColorEdit3("Background Color", (float*)&BackgroundColor);

        ImGui::End();

        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;

            if (ImGui::CollapsingHeader("Configuration"))
            {
                if (ImGui::TreeNode("Particles Kinds"))
                {
                    for (int ParticleKindIndex = 0; ParticleKindIndex < CellEngineConfigDataObject.ParticlesKinds.size(); ParticleKindIndex++)
                        ImGui::Checkbox(string(to_string(CellEngineConfigDataObject.ParticlesKinds[ParticleKindIndex].Identifier) + " " + CellEngineConfigDataObject.ParticlesKinds[ParticleKindIndex].NameFromDataFile).c_str(), &CellEngineConfigDataObject.ParticlesKinds[ParticleKindIndex].Visible);
                    ImGui::TreePop();
                }

                if (ImGui::TreeNode("Atoms"))
                {
                    for (int AtomKindIndex = 0; AtomKindIndex < CellEngineConfigDataObject.AtomsKinds.size(); AtomKindIndex++)
                        ImGui::ColorEdit3(string(CellEngineConfigDataObject.AtomsKinds[AtomKindIndex].Name + " Atom Color").c_str(), (float*)&CellEngineConfigDataObject.AtomsKinds[AtomKindIndex].Color);
                    ImGui::TreePop();
                }
            }


            static char   s8_v  = 127;
            char s8 = 1;

            ImGui::InputScalar("input long value with text s8", ImGuiDataType_S8, &s8_v, &s8, NULL, "%d");

            ImGui::End();
        }







        ImGui::Render();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(BackgroundColor.x * BackgroundColor.w, BackgroundColor.y * BackgroundColor.w, BackgroundColor.z * BackgroundColor.w, BackgroundColor.w);
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
