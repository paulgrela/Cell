
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

void DrawPlusMinusScalarButton(float& VariableToChange, const float Step, const float MinValue, const float MaxValue, const string& Description, int& IDButton)
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

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    const char* glsl_version = "#version 130";
    ImGui_ImplOpenGL3_Init(glsl_version);

    thread CellEngineOpenGLVisualiserThreadObject(CellEngineOpenGLVisualiserThreadFunction, CellEngineConfigDataObject.XTopMainWindow, CellEngineConfigDataObject.YTopMainWindow, CellEngineConfigDataObject.WidthMainWindow, CellEngineConfigDataObject.HeightMainWindow);

    static bool ModifiableWindow = false;

    while (!glfwWindowShouldClose(window))
    {
        if (CellEngineConfigDataObject.ImGuiLightVersion == true)
            ImGui::StyleColorsLight();
        else
            ImGui::StyleColorsDark();

        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (CellEngineConfigDataObject.ImGuiDemoWindowMenu == true)
            ImGui::ShowDemoWindow();

        ImGuiWindowFlags window_flags = 0;
        window_flags |= ImGuiWindowFlags_NoMove;
        window_flags |= ImGuiWindowFlags_NoResize;

        if (ModifiableWindow == false)
            ImGui::Begin("Cell Engine Visualiser", nullptr, window_flags);
        else
            ImGui::Begin("Cell Engine Visualiser");

        ImGui::Text("STATUS");

        if (ImGui::CollapsingHeader("Menu"))
        {
            ImGui::Checkbox("Menu Style Of Colors", &CellEngineConfigDataObject.ImGuiLightVersion);
            ImGui::Checkbox("Modifiable Window", &ModifiableWindow);
        }

        int IDButton = 1;

        if (ImGui::CollapsingHeader("View Move", ImGuiTreeNodeFlags_DefaultOpen))
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

        if (ImGui::CollapsingHeader("Size Of Atoms"))
        {
            DrawPlusMinusScalarButton(CellEngineConfigDataObject.SizeStep, 0.01, 0, 10, "Size Of Atoms Change Step", IDButton);

            static float SizeOfAtom3Axis = CellEngineConfigDataObject.SizeOfAtomX;
            DrawPlusMinusScalarButton(SizeOfAtom3Axis, CellEngineConfigDataObject.SizeStep, 0, 10, "Size Of Atoms 3 Axis", IDButton);
            CellEngineConfigDataObject.SizeOfAtomX = SizeOfAtom3Axis;
            CellEngineConfigDataObject.SizeOfAtomY = SizeOfAtom3Axis;
            CellEngineConfigDataObject.SizeOfAtomZ = SizeOfAtom3Axis;
        }

        if (ImGui::CollapsingHeader("Film"))
        {
            ImGui::PushID(IDButton);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(3 / 7.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(3 / 7.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(3 / 7.0f, 0.8f, 0.8f));
            if (ImGui::Button(" START "))
                CellEngineDataFileObjectPointer->StartFilmOfStructures();
            ImGui::PopStyleColor(3);
            ImGui::PopID();
            IDButton++;

            ImGui::SameLine();
            if (ImGui::Button(" NEXT "))
                CellEngineDataFileObjectPointer->ShowNextStructure();

            ImGui::SameLine();
            if (ImGui::Button(" PREV "))
                CellEngineDataFileObjectPointer->ShowPrevStructure();

            ImGui::SameLine();
            ImGui::PushID(IDButton);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0 / 7.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0 / 7.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0 / 7.0f, 0.8f, 0.8f));
            if (ImGui::Button(" STOP "))
                CellEngineDataFileObjectPointer->StopFilmOfStructures();
            ImGui::PopStyleColor(3);
            ImGui::PopID();
            IDButton++;
        }

        if (ImGui::CollapsingHeader("Shape of particles - atoms"))
        {
            ImGui::RadioButton("Sphere", &CellEngineConfigDataObject.ChosenShapeOfAtoms, 1);
            ImGui::SameLine();
            ImGui::RadioButton("Cube", &CellEngineConfigDataObject.ChosenShapeOfAtoms, 2);
            ImGui::SameLine();
            ImGui::RadioButton("Torus", &CellEngineConfigDataObject.ChosenShapeOfAtoms, 3);
        }


        if (ImGui::CollapsingHeader("Picked Atom Detailed Information", ImGuiTreeNodeFlags_DefaultOpen))
        {
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
        }

        if (ImGui::CollapsingHeader("Rendering Information", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Text("%s", CellEngineConfigDataObject.TimeParametersOfRenderingStr.c_str());
            ImGui::Text("%s", CellEngineConfigDataObject.NumberOfRenderedAtomsParametersOfRenderingStr.c_str());
            ImGui::Checkbox("Log parameters of rendering to file", &CellEngineConfigDataObject.LogParametersOfRenderingToFile);
        }

        if (ImGui::CollapsingHeader("Background"))
        {
            for (uint64_t BackgroundColorIndex = 1; BackgroundColorIndex <= 3; BackgroundColorIndex++)
            {
                ImVec4 BackgroundColor = ImVec4(CellEngineConfigDataObject.BackgroundColors[BackgroundColorIndex].X(), CellEngineConfigDataObject.BackgroundColors[BackgroundColorIndex].Y(), CellEngineConfigDataObject.BackgroundColors[BackgroundColorIndex].Z(), 1.00f);
                ImGui::ColorEdit3(string("Background Color " + to_string(BackgroundColorIndex)).c_str(), (float*)&BackgroundColor);
                CellEngineConfigDataObject.BackgroundColors[BackgroundColorIndex] = vmath::vec3(BackgroundColor.x, BackgroundColor.y, BackgroundColor.z);
            }
            const char* BackgroundColorComboBoxItems[] = { "Background Color 1", "Background Color 2", "Background Color 3" };
            static int BackgroundColorComboBoxItemsIndex = static_cast<int>(CellEngineConfigDataObject.ChosenBackgroundColor - 1);
            ImGui::Combo( " Chosen Background Color", &BackgroundColorComboBoxItemsIndex, BackgroundColorComboBoxItems, IM_ARRAYSIZE(BackgroundColorComboBoxItems));
            CellEngineConfigDataObject.ChosenBackgroundColor = BackgroundColorComboBoxItemsIndex + 1;
        }

        ImGui::End();

        if (ModifiableWindow == false)
            ImGui::Begin("Details Of Cell Drawing", nullptr, window_flags);
        else
            ImGui::Begin("Details Of Cell Drawing");

        if (ImGui::CollapsingHeader("Visibility parameters of particles", ImGuiTreeNodeFlags_DefaultOpen))
        {
            const char* DensityOfDrawedAtomsComboBoxItems[] = { "1", "10", "100", "AUTOMATIC" };
            int DensityOfDrawedAtomsItemIndex;
            switch (CellEngineConfigDataObject.LoadOfAtomsStep)
            {
                case 1 : DensityOfDrawedAtomsItemIndex = 0; break;
                case 10 : DensityOfDrawedAtomsItemIndex = 1; break;
                case 100 : DensityOfDrawedAtomsItemIndex = 2; break;
                default : break;
            }
            if (CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep == true)
                DensityOfDrawedAtomsItemIndex = 3;
            ImGui::Combo( " Density Of Drawed Atoms", &DensityOfDrawedAtomsItemIndex, DensityOfDrawedAtomsComboBoxItems, IM_ARRAYSIZE(DensityOfDrawedAtomsComboBoxItems));
            switch (DensityOfDrawedAtomsItemIndex)
            {
                case 0 : CellEngineConfigDataObject.LoadOfAtomsStep = 1; CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = false; break;
                case 1 : CellEngineConfigDataObject.LoadOfAtomsStep = 10; CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = false; break;
                case 2 : CellEngineConfigDataObject.LoadOfAtomsStep = 100; CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = false; break;
                case 4 : CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = true; break;
                default : break;
            }

            ImGui::Checkbox("Automatic Change Of Size Of Atom", &CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom);
            ImGui::Checkbox("Show Details In Atom Scale", &CellEngineConfigDataObject.ShowDetailsInAtomScale);
            ImGui::Checkbox("Draw Bonds Between Atoms", &CellEngineConfigDataObject.DrawBondsBetweenAtoms);
            ImGui::Checkbox("Draw Bonds Between Particles Centers", &CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters);
            ImGui::Checkbox("Render Objects", &CellEngineOpenGLVisualiserPointer->RenderObjectsBool);

            const char* TypesOfVisibilityComboBoxItems[] = { "ALL", "ONLY DNA", "SELECTED" };
            static int TypesOfVisibilityComboBoxCurrentItemIndex = 0;
            ImGui::Combo( " Types of Visibility", &TypesOfVisibilityComboBoxCurrentItemIndex, TypesOfVisibilityComboBoxItems, IM_ARRAYSIZE(TypesOfVisibilityComboBoxItems));
            switch (TypesOfVisibilityComboBoxCurrentItemIndex)
            {
                case 0 : CellEngineOpenGLVisualiserPointer->SetVisibilityOfAllParticles(true); break;
                case 1 : CellEngineOpenGLVisualiserPointer->SetVisibilityOfParticlesExcept(694, false); break;
                default : break;
            }

            if (TypesOfVisibilityComboBoxCurrentItemIndex == 2)
                if (ImGui::TreeNode("Particles Kinds"))
                {
                    for (auto& ParticlesKind : CellEngineConfigDataObject.ParticlesKinds)
                        ImGui::Checkbox(string(to_string(ParticlesKind.Identifier) + " " + ParticlesKind.NameFromDataFile).c_str(), &ParticlesKind.Visible);
                    ImGui::TreePop();
                }

            if (ImGui::TreeNode("Atoms"))
            {
                for (auto& AtomsKind : CellEngineConfigDataObject.AtomsKinds)
                    ImGui::ColorEdit3(string(AtomsKind.Name + " Atom Color").c_str(), (float*)&AtomsKind.Color);

                for (auto& ParticleCenter : CellEngineDataFileObjectPointer->GetParticlesCenters())
                {
                    ParticleCenter.AtomColor = CellEngineConfigDataObject.GetAtomKindDataForAtom(ParticleCenter.Name[0])->Color;
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        //DODATKOWY WARUNEK
                        for (auto& AtomObject : CellEngineDataFileObjectPointer->GetAllAtoms()[ParticleCenter.AtomIndex])
                            AtomObject.AtomColor = CellEngineConfigDataObject.GetAtomKindDataForAtom(AtomObject.Name[0])->Color;
                }

                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Type of generated colors"))
            {
                static int MakeColorsTypeData = static_cast<int>(CellEngineConfigDataObject.MakeColorsTypeObject);
                ImGui::RadioButton("Draw Color For Every Atom", &MakeColorsTypeData, 1);
                ImGui::RadioButton("Draw Color For Every Particle", &MakeColorsTypeData, 2);
                ImGui::RadioButton("Draw Random ColorFor Every Particle", &MakeColorsTypeData, 3);
                CellEngineConfigDataObject.MakeColorsTypeObject = static_cast<CellEngineConfigData::MakeColorsType>(MakeColorsTypeData);
                ImGui::TreePop();
            }
        }

        ImGui::End();


        ImGui::Render();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        ImVec4 BackgroundColor = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
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
