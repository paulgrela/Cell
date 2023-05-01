
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "DestinationPlatform.h"

#include <string>
#include <memory>

#include "Logger.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineUseful.h"
#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"
#include "CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h"
#include "CellEngineConfigurationFileReaderWriter.h"

#include "CellEngineWellStirredChemicalReactionsSimulation.h"

using namespace std;

class CellEngineImGuiMenu
{
public:
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
    APPINFO Info{};

public:
    static void InitializeLoggerManagerParameters()
    {
        try
        {
            const char* LogDirectory = "." OS_DIR_SEP;
            LoggersManagerObject.InitializeFilesNames({ "AllMessages" });
            LoggersManagerObject.InitializeSelectiveWordsFunctions({ [](const string& s) { return true; } });
            LoggersManagerObject.InitializeLoggerManagerDataForTask("CELL_RESULTS", LogDirectory, string("Logs." + GetActualDateTimeStandardCPP(".", ".", ".", ".", ".")), true, 0, function<void(const UnsignedInt& CurrentThreadId, const UnsignedInt FileNumber, const string& MessageStr)>());
            LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);
        }
        CATCH("initializing logger manager parameters")
    }

    static void ReadInitConfiguration(int argc, const char** argv)
    {
        try
        {
            InitializeLoggerManagerParameters();
            LoggersManagerObject.Log(STREAM("START CELL"));

            UnsignedInt ExecuteCellStateId = 1;
            if (argc > 1)
                ExecuteCellStateId = stoi(argv[1]);
            else
                LoggersManagerObject.Log(STREAM("Lack of cell id to execute in program parameters"));

            CellEngineConfigurationFileReaderWriterObject.ReadCellConfigurationFile("CellEngineProjectConfig.xml", ExecuteCellStateId);
        }
        CATCH("reading of data file")
    }

    static void ColorButton(const char* Text, float& VariableToChange, const float Step, const float MinValue, const float MaxValue, const float ColorParam, int& IDButton, void (*FunctionToExecute)(float& VariableToChange, const float Step, const float MinValue, const float MaxValue))
    {
        try
        {
            ImGui::PushID(IDButton);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(ColorParam / 7.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(ColorParam / 7.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(ColorParam / 7.0f, 0.8f, 0.8f));
            if (ImGui::Button(Text))
                FunctionToExecute(VariableToChange, Step, MinValue, MaxValue);
            ImGui::PopStyleColor(3);
            ImGui::PopID();
            IDButton++;
        }
        CATCH("drawing color button");
    }

    static void DrawPlusMinusScalarButton(float& VariableToChange, const float Step, const float MinValue, const float MaxValue, const string& Description, int& IDButton)
    {
        try
        {
            ColorButton(" - ", VariableToChange, Step, MinValue, MaxValue, 0, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue){ if (VariableToChange - Step >= MinValue) VariableToChange -= Step; });
            ImGui::SameLine();
            ColorButton(" + ", VariableToChange, Step, MinValue, MaxValue, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue){ if (VariableToChange + Step <= MaxValue) VariableToChange += Step; });
            ImGui::SameLine();
            ImGui::Text("%s", string(to_string(VariableToChange) + " [" + Description + "]").c_str());
        }
        CATCH("drawing plus minus scalar button");
    }

    static void glfw_error_callback(int Error, const char* Description)
    {
        LoggersManagerObject.Log(STREAM("Glfw Error nr " << Error << " : " << Description << endl));
    }

    GLFWwindow* PrepareImGuiMenuGLFWData()
    {
        GLFWwindow* ImGuiMenuWindow;

        try
        {
            glfwSetErrorCallback(glfw_error_callback);
            if (!glfwInit())
                #ifdef WINDOWS_PLATFORM
                ExitProcess(1);
                #endif
                #ifdef UNIX_PLATFORM
                exit(1);
                #endif

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

            ImGuiMenuWindow = glfwCreateWindow(CellEngineConfigDataObject.WidthMenuWindow, CellEngineConfigDataObject.HeightMenuWindow, "Dear ImGui GLFW+OpenGL3 example", NULL, NULL);
            if (ImGuiMenuWindow == nullptr)
                #ifdef WINDOWS_PLATFORM
                ExitProcess(1);
                #endif
                #ifdef UNIX_PLATFORM
                exit(1);
                #endif

            glfwSetWindowPos(ImGuiMenuWindow, CellEngineConfigDataObject.XTopMenuWindow, CellEngineConfigDataObject.YTopMenuWindow);

            if (!Info.Flags.Cursor)
                glfwSetInputMode(ImGuiMenuWindow, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

            gl3wInit();

            glfwMakeContextCurrent(ImGuiMenuWindow);

            glfwSwapInterval(1);

            IMGUI_CHECKVERSION();
            ImGui::CreateContext();
            ImGuiIO& io = ImGui::GetIO();

            ImGui_ImplGlfw_InitForOpenGL(ImGuiMenuWindow, true);
            const char* glsl_version = "#version 130";
            ImGui_ImplOpenGL3_Init(glsl_version);
        }
        CATCH("preparing imgui menu glfw data");

        return ImGuiMenuWindow;
    };

public:
    static void MenuParametersMenu(bool& ModifiableWindow)
    {
        try
        {
            if (ImGui::CollapsingHeader("Menu Parameters"))
            {
                ImGui::Checkbox("Menu Type Of Colors - Light Style", &CellEngineConfigDataObject.ImGuiLightVersion);
                ImGui::Checkbox("Modifiable Window", &ModifiableWindow);
            }
        }
        CATCH("executing menu parameters menu");
    }

    static void ViewMoveMenu(int& IDButton)
    {
        try
        {
            if (ImGui::CollapsingHeader("View Move", ImGuiTreeNodeFlags_DefaultOpen))
            {
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewXMoveShortStep, 1, 1, 10, "View X Move Short Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewYMoveShortStep, 1, 1, 10, "View Y Move Short Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewZMoveShortStep, 1, 1, 10, "View Z Move Short Step", IDButton);

                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewPositionX, CellEngineConfigDataObject.ViewXMoveShortStep, -3000, 3000, "View X Change Using Short Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewPositionY, CellEngineConfigDataObject.ViewYMoveShortStep, -3000, 3000, "View Y Change Using Short Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewPositionZ, CellEngineConfigDataObject.ViewZMoveShortStep, -3000, 3000, "View Z Change Using Short Step", IDButton);

                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewXMoveLongStep, 10, 0, 100, "View X Move Long Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewYMoveLongStep, 10, 0, 100, "View Y Move Long Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewZMoveLongStep, 10, 0, 100, "View Z Move Long Step", IDButton);

                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewPositionX, CellEngineConfigDataObject.ViewXMoveLongStep, -3000, 3000, "View X Change Using Long Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewPositionY, CellEngineConfigDataObject.ViewYMoveLongStep, -3000, 3000, "View Y Change Using Long Step", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.ViewPositionZ, CellEngineConfigDataObject.ViewZMoveLongStep, -3000, 3000, "View Z Change Using Long Step", IDButton);

                ImGui::Checkbox("View Change Using Long Step", &CellEngineConfigDataObject.ViewChangeUsingLongStep);
            }
        }
        CATCH("executing view move menu");
    }

    static void CameraMoveMenu(int& IDButton)
    {
        try
        {
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
        }
        CATCH("executing camera move menu");
    }

    static void RotationAngleMenu(int& IDButton)
    {
        try
        {
            if (ImGui::CollapsingHeader("Rotation Angle"))
            {
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.RotationAngle1, 1, -360, 360, "Rotation Angle 1", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.RotationAngle2, 1, -360, 360, "Rotation Angle 2", IDButton);
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.RotationAngle3, 1, -360, 360, "Rotation Angle 3", IDButton);
            }
        }
        CATCH("executing rotation angle menu");
    }

    static void SizeOfAtomsMenu(int& IDButton)
    {
        try
        {
            if (ImGui::CollapsingHeader("Size Of Atoms"))
            {
                DrawPlusMinusScalarButton(CellEngineConfigDataObject.SizeOfAtomChangeStep, 0.01, 0, 10, "Size Of Atoms Change Step", IDButton);

                static float SizeOfAtom3Axis = CellEngineConfigDataObject.SizeOfAtomX;
                DrawPlusMinusScalarButton(SizeOfAtom3Axis, CellEngineConfigDataObject.SizeOfAtomChangeStep, 0, 10, "Size Of Atoms 3 Axis", IDButton);
                CellEngineConfigDataObject.SizeOfAtomX = SizeOfAtom3Axis;
                CellEngineConfigDataObject.SizeOfAtomY = SizeOfAtom3Axis;
                CellEngineConfigDataObject.SizeOfAtomZ = SizeOfAtom3Axis;
            }
        }
        CATCH("executing size of atoms menu");
    }

    static void FilmMenu(int& IDButton)
    {
        try
        {
            if (ImGui::CollapsingHeader("Film"))
            {
                float Nothing;
                ColorButton(" START ", Nothing, 0, 0, 0, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue){ CellEngineDataFileObjectPointer->StartFilmOfStructures(); });

                ImGui::SameLine();
                if (ImGui::Button(" NEXT "))
                    CellEngineDataFileObjectPointer->ShowNextStructure();

                ImGui::SameLine();
                if (ImGui::Button(" PREV "))
                    CellEngineDataFileObjectPointer->ShowPrevStructure();

                ImGui::SameLine();

                ColorButton(" STOP ", Nothing, 0, 0, 0, 0, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue){ CellEngineDataFileObjectPointer->StopFilmOfStructures(); });
            }
        }
        CATCH("executing film menu");
    }

    static void ShapeOfAtomsMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Shape of particles - atoms"))
            {
                ImGui::RadioButton("Sphere", &CellEngineConfigDataObject.ChosenShapeOfAtoms, 1);
                ImGui::SameLine();
                ImGui::RadioButton("Cube", &CellEngineConfigDataObject.ChosenShapeOfAtoms, 2);
                ImGui::SameLine();
                ImGui::RadioButton("Torus", &CellEngineConfigDataObject.ChosenShapeOfAtoms, 3);
            }
        }
        CATCH("executing shape of atoms men");
    }

    static void PickedAtomDetailedInformationMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Picked Atom Detailed Information", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::Checkbox("Show details of picked atom", &CellEngineConfigDataObject.UseStencilBuffer);

                ImGui::Checkbox("Print atom description on screen", &CellEngineConfigDataObject.PrintAtomDescriptionOnScreen);

                if (CellEngineConfigDataObject.NumberOfStencilBufferLoops == 3)
                {
                    ImGui::Text("ATOM DATA:");
                    ImGui::Text("%s", string(CellEngineConfigDataObject.AtomDescriptionStr1 + " " + CellEngineConfigDataObject.AtomDescriptionStr2).c_str());
                    ImGui::Text("%s", string(CellEngineConfigDataObject.AtomDescriptionStr3 + " " + CellEngineConfigDataObject.AtomDescriptionStr4).c_str());
                }
            }
        }
        CATCH("executing picket atom detailed information");
    }

    static void RenderingInformationMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Rendering Information", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::Text("%s", CellEngineConfigDataObject.TimeParametersOfRenderingStr.c_str());
                ImGui::Text("%s", CellEngineConfigDataObject.NumberOfRenderedAtomsParametersOfRenderingStr.c_str());
                ImGui::Checkbox("Log parameters of rendering to file", &CellEngineConfigDataObject.LogParametersOfRenderingToFile);
            }
        }
        CATCH("executing rendering information menu");
    }

    static void BackgroundMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Background"))
            {
                for (UnsignedInt BackgroundColorIndex = 1; BackgroundColorIndex <= 3; BackgroundColorIndex++)
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
        }
        CATCH("executing background menu");
    }

    static void MenuWindow1(ImGuiWindowFlags WindowFlags, bool& ModifiableWindow)
    {
        try
        {
            if (ModifiableWindow == false)
                ImGui::Begin("Cell Engine Visualiser", nullptr, WindowFlags);
            else
                ImGui::Begin("Cell Engine Visualiser");

            ImGui::Text("STATUS");

            MenuParametersMenu(ModifiableWindow);

            int IDButton = 1;

            ViewMoveMenu(IDButton);

            CameraMoveMenu(IDButton);

            RotationAngleMenu(IDButton);

            SizeOfAtomsMenu(IDButton);

            FilmMenu(IDButton);

            ShapeOfAtomsMenu();

            PickedAtomDetailedInformationMenu();

            RenderingInformationMenu();

            BackgroundMenu();

            ImGui::End();
        }
        CATCH("executing menu window 1");
    }

public:
    static void DensityOfDrawnAtomsMenu()
    {
        try
        {
            const char* DensityOfDrawnAtomsComboBoxItems[] = {"1", "10", "100", "AUTOMATIC" };
            int DensityOfDrawnAtomsItemIndex;
            switch (CellEngineConfigDataObject.LoadOfAtomsStep)
            {
                case 1 : DensityOfDrawnAtomsItemIndex = 0; break;
                case 10 : DensityOfDrawnAtomsItemIndex = 1; break;
                case 100 : DensityOfDrawnAtomsItemIndex = 2; break;
                default : break;
            }
            if (CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep == true)
                DensityOfDrawnAtomsItemIndex = 3;
            ImGui::Combo(" Density Of Drawn Atoms", &DensityOfDrawnAtomsItemIndex, DensityOfDrawnAtomsComboBoxItems, IM_ARRAYSIZE(DensityOfDrawnAtomsComboBoxItems));
            switch (DensityOfDrawnAtomsItemIndex)
            {
                case 0 : CellEngineConfigDataObject.LoadOfAtomsStep = 1; CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = false; break;
                case 1 : CellEngineConfigDataObject.LoadOfAtomsStep = 10; CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = false; break;
                case 2 : CellEngineConfigDataObject.LoadOfAtomsStep = 100; CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = false; break;
                case 3 : CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = true; break;
                default : break;
            }
        }
        CATCH("executing of density of drawn atoms menu");
    }

    void TypesOfVisibilityMenu() const
    {
        try
        {
            const char* TypesOfVisibilityComboBoxItems[] = { "ALL", "ONLY DNA", "ONLY RNA", "SELECTED" };
            static int TypesOfVisibilityComboBoxCurrentItemIndex = 0;
            ImGui::Combo( " Types of Visibility", &TypesOfVisibilityComboBoxCurrentItemIndex, TypesOfVisibilityComboBoxItems, IM_ARRAYSIZE(TypesOfVisibilityComboBoxItems));
            switch (TypesOfVisibilityComboBoxCurrentItemIndex)
            {
                case 0 : CellEngineOpenGLVisualiserPointer->SetVisibilityOfAllParticles(true); break;
                case 1 : CellEngineOpenGLVisualiserPointer->SetVisibilityOfParticlesExcept(CellEngineConfigDataObject.DNAIdentifier, false); break;
                case 2 : CellEngineOpenGLVisualiserPointer->SetVisibilityOfParticlesExcept(CellEngineConfigDataObject.RNAIdentifier, false); break;
                default : break;
            }
            if (TypesOfVisibilityComboBoxCurrentItemIndex == 3)
                if (ImGui::CollapsingHeader("Particles Kinds", ImGuiTreeNodeFlags_DefaultOpen))
                {
                    for (auto& ParticlesKind : CellEngineConfigDataObject.ParticlesKinds)
                        ImGui::Checkbox(string(to_string(ParticlesKind.Identifier) + " " + ParticlesKind.NameFromDataFile).c_str(), &ParticlesKind.Visible);
                }
        }
        CATCH("executing types of visibility menu");
    }

    static void AtomKindsMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Atoms Kinds", ImGuiTreeNodeFlags_DefaultOpen))
            {
                bool ChangeColor = false;
                for (auto& AtomsKind : CellEngineConfigDataObject.AtomsKinds)
                {
                    if (ImGui::ColorEdit3(string(AtomsKind.Name + " Atom Color").c_str(), (float*)&AtomsKind.ColorVmathVec3) == true)
                    {
                        AtomsKind.Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(AtomsKind.ColorVmathVec3);
                        ChangeColor = true;
                    }
                }

                if (ChangeColor == true)
                    for (auto& ParticleCenter : CellEngineDataFileObjectPointer->GetParticlesCenters())
                    {
                        ParticleCenter.AtomColor = CellEngineConfigDataObject.GetAtomKindDataForAtom(ParticleCenter.Name[0])->Color;
                        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                            for (auto& AtomObject : CellEngineDataFileObjectPointer->GetAllAtoms()[ParticleCenter.AtomIndex])
                                AtomObject.AtomColor = CellEngineConfigDataObject.GetAtomKindDataForAtom(AtomObject.Name[0])->Color;
                    }
            }
        }
        CATCH("executing atoms kinds menu");
    }

    static void TypeOfGeneratedColorsMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Type of generated colors", ImGuiTreeNodeFlags_DefaultOpen))
            {
                static int MakeColorsTypeData = static_cast<int>(CellEngineConfigDataObject.MakeColorsTypeObject);
                ImGui::RadioButton("Draw Color For Every Atom", &MakeColorsTypeData, 1);
                ImGui::RadioButton("Draw Color For Every Particle Kind", &MakeColorsTypeData, 2);
                ImGui::RadioButton("Draw Random Color For Every Particle Kind", &MakeColorsTypeData, 3);
                ImGui::RadioButton("Draw Random Color For Every Unique Particle", &MakeColorsTypeData, 4);
                CellEngineConfigDataObject.MakeColorsTypeObject = static_cast<CellEngineConfigData::MakeColorsType>(MakeColorsTypeData);
            }
        }
        CATCH("executing types of generated colors menu");
    }

    void VoxelSimulationSpaceVisibility() const
    {
        try
        {
            if (CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer != nullptr)
            {
                std::lock_guard<std::mutex> LockGuardObject{CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderMenuAndVoxelSimulationSpaceMutexObject};

                auto  CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer = dynamic_cast<CellEngineOpenGLVisualiserOfVoxelSimulationSpace*>(CellEngineOpenGLVisualiserPointer.get());

                const auto StartPos = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetStartPositions();
                static int DrawSpaceStartXYZ[3] = { static_cast<int>(get<0>(StartPos)), static_cast<int>(get<1>(StartPos)), static_cast<int>(get<2>(StartPos)) };
                ImGui::DragInt3("StartX StartY StartZ", DrawSpaceStartXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                const auto Steps = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetSteps();
                static int DrawSpaceStepsXYZ[3] = { static_cast<int>(get<0>(Steps)), static_cast<int>(get<1>(Steps)), static_cast<int>(get<2>(Steps)) };
                ImGui::DragInt3("StepX  StepY  StepZ", DrawSpaceStepsXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                const auto Sizes = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetSizes();
                static int DrawSpaceSizesXYZ[3] = { static_cast<int>(get<0>(Sizes)), static_cast<int>(get<1>(Sizes)), static_cast<int>(get<2>(Sizes)) };
                ImGui::DragInt3("SizeX  SizeY  SizeZ", DrawSpaceSizesXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SetVoxelSpaceSelection(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

                ImGui::Checkbox("Draw empty voxels", &CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->DrawEmptyVoxels);

                static int TypeOfDrawingVoxelSpace = static_cast<int>(CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType);
                ImGui::RadioButton("Draw Voxel Space FULL", &TypeOfDrawingVoxelSpace, 1);
                ImGui::RadioButton("Draw Voxel Space SELECTED", &TypeOfDrawingVoxelSpace, 2);
                CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType = static_cast<CellEngineOpenGLVisualiserOfVoxelSimulationSpace::VoxelSpaceDrawingTypes>(TypeOfDrawingVoxelSpace);

                if (CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType == CellEngineOpenGLVisualiserOfVoxelSimulationSpace::VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
                    CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SetVoxelSpaceSelection(0, 0, 0, 64, 64, 64, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension);
                else
                if (CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType == CellEngineOpenGLVisualiserOfVoxelSimulationSpace::VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
                    CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SetVoxelSpaceSelection(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

                if (ImGui::Button("  SAVE MOUSE POSITION  ") == true)
                {
                    CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SaveVoxelPositionChosenByMouse();

                    const auto TempStartPos = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetStartPositions();
                    DrawSpaceStartXYZ[0] = static_cast<int>(get<0>(TempStartPos));
                    DrawSpaceStartXYZ[1] = static_cast<int>(get<1>(TempStartPos));
                    DrawSpaceStartXYZ[2] = static_cast<int>(get<2>(TempStartPos));
                }

                if (ImGui::Button(" DRAW RANDOM PARTICLES ") == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomParticlesInSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

                int IDButton = 1;
                float Nothing;
                ColorButton("    START DIFFUSION    ", Nothing, 0, 0, 0, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue){ CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfDiffusion(10, 0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]); });

                ColorButton("    STOP DIFFUSION     ", Nothing, 0, 0, 0, 0, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue){ CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetDNASequenceFromNucleotides(); });
            }
        }
        CATCH("executing voxel simulation space visibility menu")
    }

    void MenuWindow2(ImGuiWindowFlags WindowFlags, const bool ModifiableWindow) const
    {
        try
        {
            if (ModifiableWindow == false)
                ImGui::Begin("Details Of Cell Drawing", nullptr, WindowFlags);
            else
                ImGui::Begin("Details Of Cell Drawing");

            if (ImGui::CollapsingHeader("Visibility parameters of particles", ImGuiTreeNodeFlags_DefaultOpen))
            {
                DensityOfDrawnAtomsMenu();

                ImGui::Checkbox("Automatic Change Of Size Of Atom", &CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom);
                ImGui::Checkbox("Show Details In Atom Scale", &CellEngineConfigDataObject.ShowDetailsInAtomScale);
                ImGui::Checkbox("Show Atoms In Each Part Of the Cell", &CellEngineConfigDataObject.ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside);
                ImGui::Checkbox("Draw Bonds Between Atoms", &CellEngineConfigDataObject.DrawBondsBetweenAtoms);
                ImGui::Checkbox("Draw Bonds Between Particles Centers", &CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters);
                ImGui::Checkbox("Render Objects", &CellEngineOpenGLVisualiserPointer->RenderObjectsBool);

                TypesOfVisibilityMenu();

                AtomKindsMenu();

                TypeOfGeneratedColorsMenu();
            }

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                if (ImGui::CollapsingHeader("Visibility of Voxel Simulation Space", ImGuiTreeNodeFlags_DefaultOpen))
                    VoxelSimulationSpaceVisibility();

            ImGui::End();
        }
        CATCH("executing menu window 2");
    }

    void ImGuiMenuGLFWMainLoop(GLFWwindow* ImGuiMenuWindow) const
    {
        try
        {
            bool ModifiableWindow = false;

            while (!glfwWindowShouldClose(ImGuiMenuWindow))
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

                ImGuiWindowFlags WindowFlags = 0;
                WindowFlags |= ImGuiWindowFlags_NoMove;
                WindowFlags |= ImGuiWindowFlags_NoResize;

                MenuWindow1(WindowFlags, ModifiableWindow);

                MenuWindow2(WindowFlags, ModifiableWindow);

                ImGui::Render();

                int ImGuiMenuWindowWidth, ImGuiMenuWindowHeight;
                glfwGetFramebufferSize(ImGuiMenuWindow, &ImGuiMenuWindowWidth, &ImGuiMenuWindowHeight);
                glViewport(0, 0, ImGuiMenuWindowWidth, ImGuiMenuWindowHeight);
                ImVec4 BackgroundColor = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
                glClearColor(BackgroundColor.x * BackgroundColor.w, BackgroundColor.y * BackgroundColor.w, BackgroundColor.z * BackgroundColor.w, BackgroundColor.w);
                glClear(GL_COLOR_BUFFER_BIT);
                ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

                glfwSwapBuffers(ImGuiMenuWindow);
            }
        }
        CATCH("executing imgui menu glfw main loop");
    }

    static void ImGuiMenuGLFShutdown(GLFWwindow* ImGuiMenuWindow)
    {
        try
        {
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext();

            glfwDestroyWindow(ImGuiMenuWindow);
            glfwTerminate();
        }
        CATCH("shutting down imgui menu")
    }

public:
    unique_ptr<CellEngineOpenGLVisualiser> CellEngineOpenGLVisualiserPointer;

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wreturn-type"
    static unique_ptr<CellEngineOpenGLVisualiser> CreateCellEngineOpenGLVisualiserObject(const CellEngineConfigData::TypesOfSpace TypeOfSpace)
    {
        switch (TypeOfSpace)
        {
            case CellEngineConfigData::TypesOfSpace::FullAtomSpace : return make_unique<CellEngineOpenGLVisualiserOfFullAtomSimulationSpace>();
            case CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace : return make_unique<CellEngineOpenGLVisualiserOfVoxelSimulationSpace>();
            default : break;
        }
    }
    #pragma GCC diagnostic pop

    void CellEngineOpenGLVisualiserThreadFunction(int XPosWindow, int YPosWindow, int WidthWindow, int HeightWindow)
    {
        try
        {
            CellEngineOpenGLVisualiserPointer = CreateCellEngineOpenGLVisualiserObject(CellEngineConfigDataObject.TypeOfSpace);
            CellEngineOpenGLVisualiserPointer->Run(XPosWindow, YPosWindow, WidthWindow, HeightWindow);
        }
        CATCH("running cell engine opengl visualiser thread function");
    }

public:
    CellEngineImGuiMenu(int argc, const char** argv)
    {
        try
        {
            ReadInitConfiguration(argc, argv);

            //TestWellStirredChemicalReactionsSimulation(); exit(0);

            GLFWwindow* ImGuiMenuWindow = PrepareImGuiMenuGLFWData();

            thread CellEngineOpenGLVisualiserThreadObject(&CellEngineImGuiMenu::CellEngineOpenGLVisualiserThreadFunction, this, CellEngineConfigDataObject.XTopMainWindow, CellEngineConfigDataObject.YTopMainWindow, CellEngineConfigDataObject.WidthMainWindow, CellEngineConfigDataObject.HeightMainWindow);

            ImGuiMenuGLFWMainLoop(ImGuiMenuWindow);

            CellEngineOpenGLVisualiserThreadObject.detach();

            ImGuiMenuGLFShutdown(ImGuiMenuWindow);
        }
        CATCH("starting imgui menu and whole cell opengl visualization")
    }
};

int main(int argc, const char ** argv)
{
    CellEngineImGuiMenu CellEngineImGuiMenuObject(argc, argv);
    return 0;
}