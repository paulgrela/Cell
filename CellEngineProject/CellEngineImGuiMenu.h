
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "DestinationPlatform.h"

#include <string>
#include <memory>

#include "Logger.h"
#include "Combinatorics.h"
#include "DateTimeUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineUseful.h"
#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"
#include "CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h"
#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineParticlesKindsManager.h"
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
                unsigned int FullScreen: 1;
                unsigned int VSync: 1;
                unsigned int Cursor: 1;
                unsigned int Stereo: 1;
                unsigned int Debug: 1;
                unsigned int Robust: 1;
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
            const char *LogDirectory = "." OS_DIR_SEP;
            LoggersManagerObject.InitializeSpecialLogFiles(false, true, false, false, false, true, false, true, false);
            LoggersManagerObject.InitializeFilesNames({"AllMessages" });
            LoggersManagerObject.InitializeSelectiveWordsFunctions({[](const string &s){ return true; } });
            LoggersManagerObject.InitializeLoggerManagerDataForTask("CELL_RESULTS", LogDirectory, string("Logs." + GetActualDateTimeStandardCPP(".", ".", ".", ".", ".")), true, 0, function<void(const UnsignedInt &CurrentThreadId, const UnsignedInt FileNumber, const string &MessageStr)>());
            LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile, CellEngineConfigDataObject.PrintLogToCommonFileWhenPrintLogToSpecialFile);
        }
        CATCH("initializing logger manager parameters")
    }

    static void ReadInitConfiguration(const int argc, const char **argv)
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

            LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile, CellEngineConfigDataObject.PrintLogToCommonFileWhenPrintLogToSpecialFile);
        }
        CATCH("reading of data file")
    }

    template<class T>
    static void ColorButton(const char *Text, float &VariableToChange, const float Step, const float MinValue, const float MaxValue, const float ColorParam, int &IDButton, T FunctionToExecute)
    {
        try
        {
            ImGui::PushID(IDButton);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4) ImColor::HSV(ColorParam / 7.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4) ImColor::HSV(ColorParam / 7.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4) ImColor::HSV(ColorParam / 7.0f, 0.8f, 0.8f));
            if (ImGui::Button(Text))
                FunctionToExecute(VariableToChange, Step, MinValue, MaxValue);
            ImGui::PopStyleColor(3);
            ImGui::PopID();
            IDButton++;
        }
        CATCH("drawing color button");
    }

    static void DrawPlusMinusScalarButton(float &VariableToChange, const float Step, const float MinValue, const float MaxValue, const string &Description, int &IDButton)
    {
        try
        {
            ColorButton(" - ", VariableToChange, Step, MinValue, MaxValue, 0, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){ if (VariableToChange - Step >= MinValue) VariableToChange -= Step; });
            ImGui::SameLine();
            ColorButton(" + ", VariableToChange, Step, MinValue, MaxValue, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){ if (VariableToChange + Step <= MaxValue) VariableToChange += Step; });
            ImGui::SameLine();
            ImGui::Text("%s", string(to_string(VariableToChange) + " [" + Description + "]").c_str());
        }
        CATCH("drawing plus minus scalar button");
    }

    static void glfw_error_callback(int Error, const char *Description)
    {
        LoggersManagerObject.Log(STREAM("Glfw Error nr " << Error << " : " << Description << endl));
    }

    GLFWwindow *PrepareImGuiMenuGLFWData()
    {
        GLFWwindow *ImGuiMenuWindow;

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
            ImGuiIO &io = ImGui::GetIO();

            ImGui_ImplGlfw_InitForOpenGL(ImGuiMenuWindow, true);
            const char *glsl_version = "#version 130";
            ImGui_ImplOpenGL3_Init(glsl_version);
        }
        CATCH("preparing imgui menu glfw data");

        return ImGuiMenuWindow;
    };

public:
    static void MenuParametersMenu(bool &ModifiableWindow)
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

    static void ViewMoveMenu(int &IDButton)
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

    static void CameraMoveMenu(int &IDButton)
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

    static void RotationAngleMenu(int &IDButton)
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

    static void SizeOfAtomsMenu(int &IDButton)
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

    static void FilmMenu(int &IDButton)
    {
        try
        {
            if (ImGui::CollapsingHeader("Film"))
            {
                float Nothing;
                ColorButton(" START ", Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){ CellEngineDataFileObjectPointer->StartFilmOfStructures(); });

                ImGui::SameLine();
                if (ImGui::Button(" NEXT "))
                    CellEngineDataFileObjectPointer->ShowNextStructure();

                ImGui::SameLine();
                if (ImGui::Button(" PREV "))
                    CellEngineDataFileObjectPointer->ShowPrevStructure();

                ImGui::SameLine();

                ColorButton(" STOP ", Nothing, 0, 0, 0, 0, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){ CellEngineDataFileObjectPointer->StopFilmOfStructures(); });
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
                    ImGui::Text("%s", string(CellEngineUseful::AtomDescriptionTextsObject.Texts[0] + " " + CellEngineUseful::AtomDescriptionTextsObject.Texts[1]).c_str());
                    ImGui::Text("%s", string(CellEngineUseful::AtomDescriptionTextsObject.Texts[2] + " " + CellEngineUseful::AtomDescriptionTextsObject.Texts[3]).c_str());
                    ImGui::Text("%s", string(CellEngineUseful::AtomDescriptionTextsObject.Texts[4] + " " + CellEngineUseful::AtomDescriptionTextsObject.Texts[5]).c_str());
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
                    ImGui::ColorEdit3(string("Background Color " + to_string(BackgroundColorIndex)).c_str(), (float *) &BackgroundColor);
                    CellEngineConfigDataObject.BackgroundColors[BackgroundColorIndex] = vmath::vec3(BackgroundColor.x, BackgroundColor.y, BackgroundColor.z);
                }
                const char *BackgroundColorComboBoxItems[] = { "Background Color 1", "Background Color 2", "Background Color 3" };
                static int BackgroundColorComboBoxItemsIndex = static_cast<int>(CellEngineConfigDataObject.ChosenBackgroundColor - 1);
                ImGui::Combo(" Chosen Background Color", &BackgroundColorComboBoxItemsIndex, BackgroundColorComboBoxItems, IM_ARRAYSIZE(BackgroundColorComboBoxItems));
                CellEngineConfigDataObject.ChosenBackgroundColor = BackgroundColorComboBoxItemsIndex + 1;
            }
        }
        CATCH("executing background menu");
    }

    static void MenuWindow1(const ImGuiWindowFlags WindowFlags, bool& ModifiableWindow)
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
            const char *DensityOfDrawnAtomsComboBoxItems[] = { "1", "10", "100", "AUTOMATIC" };
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

    static void DetailedVisibilityParametersOfParticles()
    {
        try
        {
            ImGui::Checkbox("Automatic Change Of Size Of Atom", &CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom);
            ImGui::Checkbox("Show Details In Atom Scale", &CellEngineConfigDataObject.ShowDetailsInAtomScale);
            ImGui::Checkbox("Show Atoms In Each Part Of the Cell", &CellEngineConfigDataObject.ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside);
            ImGui::Checkbox("Draw Bonds Between Atoms", &CellEngineConfigDataObject.DrawBondsBetweenAtoms);
            ImGui::Checkbox("Render Objects", &CellEngineOpenGLVisualiserPointer->RenderObjectsBool);
        }
        CATCH("drawing visibility parameters of particles")
    }

    static bool ComparisonOfParticle(const ParticleKind& P1, const ParticleKind& P2)
    {
        if (P1.ParticleKindSpecialDataSector.empty() == false && P2.ParticleKindSpecialDataSector.empty() == false)
        {
            if (ParticlesKindsManager::ConvertParticleTypeToString(P1.ParticleKindSpecialDataSector.back().ParticleType) < ParticlesKindsManager::ConvertParticleTypeToString(P2.ParticleKindSpecialDataSector.back().ParticleType))
                return true;
            if (ParticlesKindsManager::ConvertParticleTypeToString(P1.ParticleKindSpecialDataSector.back().ParticleType) == ParticlesKindsManager::ConvertParticleTypeToString(P2.ParticleKindSpecialDataSector.back().ParticleType))
                return P1.IdStr < P2.IdStr;
        }
        return false;
    }

    static void TypesOfVisibilityMenu()
    {
        try
        {
            const char *TypesOfVisibilityComboBoxItems[] = { "ALL", "ONLY DNA", "ONLY RNA", "SELECTED" };

            static int TypesOfVisibilityComboBoxCurrentItemIndex = 0;
            static int PrevTypesOfVisibilityComboBoxCurrentItemIndex = 0;

            ImGui::Combo(" Types of Visibility", &TypesOfVisibilityComboBoxCurrentItemIndex, TypesOfVisibilityComboBoxItems, IM_ARRAYSIZE(TypesOfVisibilityComboBoxItems));
            if (PrevTypesOfVisibilityComboBoxCurrentItemIndex != TypesOfVisibilityComboBoxCurrentItemIndex)
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
                    vector<ParticleKind> LocalParticlesKindsForImGuiSelection;

                        for (const auto& ParticlesKindsMapElement : ParticlesKindsManagerObject.ParticlesKinds)
                            if (ParticlesKindsMapElement.second.IdStr.starts_with("JCVISYN3A_") || ParticlesKindsMapElement.second.IdStr.starts_with("particle_") || ParticlesKindsMapElement.second.IdStr.starts_with("trna_") || ParticlesKindsMapElement.second.IdStr.starts_with("mrna_") || ParticlesKindsMapElement.second.IdStr.starts_with("rrna_") || ParticlesKindsMapElement.second.IdStr.starts_with("M_"))
                                LocalParticlesKindsForImGuiSelection.emplace_back(ParticlesKindsMapElement.second);

                        sort(LocalParticlesKindsForImGuiSelection.begin(), LocalParticlesKindsForImGuiSelection.end(), ComparisonOfParticle);

                        for (auto &ParticlesKindForImGuiSelectionObject : LocalParticlesKindsForImGuiSelection)
                        {
                            string TypeStr = "Unknown Type ";
                            if (ParticlesKindForImGuiSelectionObject.ParticleKindSpecialDataSector.empty() == false)
                                TypeStr = ParticlesKindsManager::ConvertParticleTypeToString(ParticlesKindForImGuiSelectionObject.ParticleKindSpecialDataSector.back().ParticleType);
                            ImGui::Checkbox(string(to_string(ParticlesKindForImGuiSelectionObject.EntityId) + " " + "T = " + TypeStr + " " + ParticlesKindForImGuiSelectionObject.IdStr + " " + (ParticlesKindForImGuiSelectionObject.Name.starts_with("M_") ? ParticlesKindForImGuiSelectionObject.Name : "") + (ParticlesKindForImGuiSelectionObject.Formula != ParticlesKindForImGuiSelectionObject.IdStr && ParticlesKindForImGuiSelectionObject.Formula.empty() == false ? ("F = " + ParticlesKindForImGuiSelectionObject.Formula + " ") : "")).c_str(), &ParticlesKindForImGuiSelectionObject.GraphicData.Visible);
                        }

                    for (const auto &ParticleKindForImGuiSelectionObject: LocalParticlesKindsForImGuiSelection)
                        ParticlesKindsManagerObject.GetParticleKind(ParticleKindForImGuiSelectionObject.EntityId).GraphicData.Visible = ParticleKindForImGuiSelectionObject.GraphicData.Visible;
                }
            PrevTypesOfVisibilityComboBoxCurrentItemIndex = TypesOfVisibilityComboBoxCurrentItemIndex;
        }
        CATCH("executing types of visibility menu");
    }

    static void ChangeNumberOfParticlesMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("Particles Kinds Add / Remove", ImGuiTreeNodeFlags_None))
            {
                if (CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer != nullptr)
                {
                    vector<ParticleKind> LocalParticlesKindsForImGuiSelection;
                    for (const auto& ParticlesKindsMapElement : ParticlesKindsManagerObject.ParticlesKinds)
                        if (ParticlesKindsMapElement.second.IdStr.starts_with("JCVISYN3A_") || ParticlesKindsMapElement.second.IdStr.starts_with("particle_") || ParticlesKindsMapElement.second.IdStr.starts_with("trna_") || ParticlesKindsMapElement.second.IdStr.starts_with("mrna_") || ParticlesKindsMapElement.second.IdStr.starts_with("rrna_") || ParticlesKindsMapElement.second.IdStr.starts_with("M_"))
                            LocalParticlesKindsForImGuiSelection.emplace_back(ParticlesKindsMapElement.second);

                    sort(LocalParticlesKindsForImGuiSelection.begin(), LocalParticlesKindsForImGuiSelection.end(), ComparisonOfParticle);

                    vector<bool> LocalParticleKindSelection;
                    for (auto &ParticlesKindForImGuiSelectionObject : LocalParticlesKindsForImGuiSelection)
                        if (ParticlesKindForImGuiSelectionObject.ParticleKindSpecialDataSector.empty() == false)
                            ImGui::Checkbox(string(to_string(ParticlesKindForImGuiSelectionObject.EntityId) + " " + "T = " + ParticlesKindsManager::ConvertParticleTypeToString(ParticlesKindForImGuiSelectionObject.ParticleKindSpecialDataSector.back().ParticleType) + " " + ParticlesKindForImGuiSelectionObject.IdStr + " " + (ParticlesKindForImGuiSelectionObject.Name.starts_with("M_") ? ParticlesKindForImGuiSelectionObject.Name : "") + (ParticlesKindForImGuiSelectionObject.Formula != ParticlesKindForImGuiSelectionObject.IdStr && ParticlesKindForImGuiSelectionObject.Formula.empty() == false ? ("F = " + ParticlesKindForImGuiSelectionObject.Formula + " ") : "")).c_str(), &ParticlesKindForImGuiSelectionObject.GraphicData.Selected);

                    for (const auto &ParticleKindForImGuiSelectionObject : LocalParticlesKindsForImGuiSelection)
                        ParticlesKindsManagerObject.GetParticleKind(ParticleKindForImGuiSelectionObject.EntityId).GraphicData.Selected = ParticleKindForImGuiSelectionObject.GraphicData.Selected;

                    static int NumberOfParticlesToSub[1] = { 100 };
                    ImGui::DragInt("SUB", NumberOfParticlesToSub, 1, 0, 1000, "%d", ImGuiSliderFlags_AlwaysClamp);
                    if (ImGui::Button("SUB PARTICLES") == true)
                        for (const auto &ParticleKindForImGuiSelectionObject: LocalParticlesKindsForImGuiSelection)
                            if (ParticleKindForImGuiSelectionObject.GraphicData.Selected == true)
                                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RemoveParticlesWithChosenEntityId(ParticleKindForImGuiSelectionObject.EntityId, NumberOfParticlesToSub[0]);

                    if (ImGui::Button("REMOVE ALL PARTICLES") == true)
                        for (const auto &ParticleKindForImGuiSelectionObject: LocalParticlesKindsForImGuiSelection)
                            if (ParticleKindForImGuiSelectionObject.GraphicData.Selected == true)
                                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RemoveAllParticlesWithChosenEntityId(ParticleKindForImGuiSelectionObject.EntityId);

                    if (ImGui::Button("REMOVE ALL RNA PARTICLES ") == true)
                        CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RemoveAllRNAParticles();
                    if (ImGui::Button("REMOVE ALL tRNA PARTICLES ") == true)
                        CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RemoveAlltRNAParticles();
                    if (ImGui::Button("REMOVE ALL mRNA PARTICLES ") == true)
                        CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RemoveAllmRNAParticles();
                    if (ImGui::Button("REMOVE ALL rRNA PARTICLES ") == true)
                        CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RemoveAllrRNAParticles();
                }
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
                for (auto &AtomsKind: ParticlesKindsManagerObject.AtomsKindsGraphicData)
                {
                    if (ImGui::ColorEdit3(string(AtomsKind.Name + " Atom Color").c_str(), (float *) &AtomsKind.ColorVmathVec3) == true)
                    {
                        AtomsKind.Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(AtomsKind.ColorVmathVec3);
                        ChangeColor = true;
                    }
                }
                if (ChangeColor == true)
                    FOR_EACH_PARTICLE_IN_XYZ_ONLY
                        for (auto& ParticleObject: CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                            for (auto& AtomObject : ParticleObject.second.ListOfAtoms)
                                AtomObject.AtomColor = ParticlesKindsManagerObject.GetGraphicAtomKindDataFromAtomName(AtomObject.Name[0])->Color;
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

    static string AlignString(const string& InputStr, const UnsignedInt Size, UnsignedInt PrefixSize = 3)
    {
        return string(PrefixSize, ' ') + InputStr + string(Size - InputStr.length(), ' ');
    }

    static void VoxelSimulationSpaceParametersMenu(CellEngineOpenGLVisualiserOfVoxelSimulationSpace* CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer, int DrawSpaceStartXYZ[3], int DrawSpaceStepsXYZ[3], int DrawSpaceSizesXYZ[3], const UnsignedInt StringLength)
    {
        try
        {
            CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SetVoxelSpaceSelection(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            ImGui::Checkbox("Draw empty voxels", &CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->DrawEmptyVoxels);

            static int TypeOfDrawingVoxelSpace = static_cast<int>(CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType);
            ImGui::RadioButton("Draw Voxel Space FULL", &TypeOfDrawingVoxelSpace, 1);
            ImGui::RadioButton("Draw Voxel Space SELECTED", &TypeOfDrawingVoxelSpace, 2);
            CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType = static_cast<CellEngineOpenGLVisualiserOfVoxelSimulationSpace::VoxelSpaceDrawingTypes>(TypeOfDrawingVoxelSpace);

            if (CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType == CellEngineOpenGLVisualiserOfVoxelSimulationSpace::VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
                CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SetVoxelSpaceSelection(0, 0, 0, 64, 64, 64, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension);
            else
            if (CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SpaceDrawingType == CellEngineOpenGLVisualiserOfVoxelSimulationSpace::VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
                CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SetVoxelSpaceSelection(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            ImGui::Text("");

            static int SelectedSpaceStartParametersDrawTypesIndex = static_cast<int>(CellEngineConfigDataObject.SelectedSpaceStartParametersDrawTypesObject);
            ImGui::RadioButton("Draw Selected Space From Center", &SelectedSpaceStartParametersDrawTypesIndex, 1);
            ImGui::RadioButton("Draw Selected Space From Corner", &SelectedSpaceStartParametersDrawTypesIndex, 2);
            CellEngineConfigDataObject.SelectedSpaceStartParametersDrawTypesObject = static_cast<CellEngineConfigData::SelectedSpaceStartParametersDrawTypes>(SelectedSpaceStartParametersDrawTypesIndex);

            ImGui::Text("%s", string("Number of free indexes for particles = " + to_string(CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetFreeIndexesOfParticleSize())).c_str());

            if (ImGui::Button(AlignString("SAVE MOUSE POSITION", StringLength).c_str()) == true)
            {
                CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->SaveVoxelPositionChosenByMouse();

                const auto TempStartPos = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetStartPositions();
                DrawSpaceStartXYZ[0] = static_cast<int>(get<0>(TempStartPos));
                DrawSpaceStartXYZ[1] = static_cast<int>(get<1>(TempStartPos));
                DrawSpaceStartXYZ[2] = static_cast<int>(get<2>(TempStartPos));
            }
        }
        CATCH("modification of voxel simulation space parameters ")
    }

    static void OftenOperationsVoxelSimulationSpaceParametersMenu(int DrawSpaceStartXYZ[3], int DrawSpaceStepsXYZ[3], int DrawSpaceSizesXYZ[3], const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::Button(AlignString("READ GENES", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadGenes();
            if (ImGui::Button(AlignString("FIND INTER GENES", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->FindInterGenesSequences();

            if (ImGui::Button(AlignString("ADD ILLINOIS DATA", StringLength).c_str()) == true)
            {
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadAllIllinoisDataFromFiles();
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeSequenceFromFile(CellEngineConfigDataObject.DNAPaired);
            }
            if (ImGui::Button(AlignString("READ ILLINOIS CHEMICAL REACTIONS", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadChemicalReactionsFromFiles();
            if (ImGui::Button(AlignString("READ UPDATED ILLINOIS CHEMICAL REACTIONS", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadNewChemicalReactionsFromFiles();

            if (ImGui::Button(AlignString("ADD SPECIAL PARTICLE KINDS", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->AddSpecialParticlesKinds();
            if (ImGui::Button(AlignString("ADD PARTICLE KINDS", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->AddParticlesKinds();
            if (ImGui::Button(AlignString("ADD CHEMICAL REACTIONS", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->AddChemicalReactions();
            if (ImGui::Button(AlignString("ADD TEST CHEMICAL REACTIONS", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->AddTestChemicalReactions();
            if (ImGui::Button(AlignString("CLEAR SELECTED SPACE", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ClearSelectedSpace(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("ADD RANDOM PARTICLES", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomParticlesInSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("ADD PLANED CUBOID PARTICLES", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GeneratePlanedCuboidParticlesInSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("ADD PLANED ELLIPSOID PARTICLES", StringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GeneratePlanedEllipsoidParticlesInSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("ADD ALL FOR TEST DNA REACTIONS", StringLength).c_str()) == true)
            {
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->AddParticlesKinds();
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->AddTestChemicalReactions();
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ClearSelectedSpace(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GeneratePlanedCuboidParticlesInSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceStepsXYZ[0], DrawSpaceStepsXYZ[1], DrawSpaceStepsXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeDataFromFile(CellEngineConfigDataObject.DNAPaired);
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeSequenceFromFile(CellEngineConfigDataObject.DNAPaired);
            }
        }
        CATCH("modification of diffusion voxel simulation space parameters menu")
    }

    static void DiffusionVoxelSimulationSpaceParametersMenu(int DrawSpaceStartXYZ[3], int DrawSpaceStepsXYZ[3], int DrawSpaceSizesXYZ[3], const UnsignedInt StringLength)
    {
        try
        {
            int IDButton = 1;
            float Nothing;

            if (ImGui::CollapsingHeader("WORK", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ColorButton(AlignString("START WORK", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});
                ColorButton(AlignString("STOP WORK", StringLength).c_str(), Nothing, 0, 0, 0, 0, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});
            }

            if (ImGui::CollapsingHeader("DIFFUSION", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ColorButton(AlignString("START DIFFUSION", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});
                ColorButton(AlignString("STOP DIFFUSION", StringLength).c_str(), Nothing, 0, 0, 0, 0, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});

                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION FOR RANGE OF PARTICLES", StringLength).c_str(), Nothing, 0, 0, 0, 8, IDButton, [DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfDiffusionForSelectedRangeOfParticles(10, 0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                });

                constexpr UnsignedInt AdditionalSpaceBoundFactor = 20;
                constexpr double MultiplyElectricChargeFactor = 100;

                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION FOR SELECTED SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 8, IDButton, [MultiplyElectricChargeFactor, DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfDiffusionForSelectedSpace(true, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                });

                ColorButton(AlignString("MAKE ONE STEP OF ELECTRIC DIFFUSION FOR RANGE OF PARTICLES - FCP", StringLength).c_str(), Nothing, 0, 0, 0, 9, IDButton, [MultiplyElectricChargeFactor, DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity::FromChosenParticleAsCenter, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, 19, 0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                });
                ColorButton(AlignString("MAKE ONE STEP OF ELECTRIC DIFFUSION FOR RANGE OF PARTICLES - ISS", StringLength).c_str(), Nothing, 0, 0, 0, 9, IDButton, [MultiplyElectricChargeFactor, DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity::InChosenSectorOfSimulationSpace, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, 19, 0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                });

                ColorButton(AlignString("MAKE ONE STEP OF ELECTRIC DIFFUSION FOR SELECTED SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 9, IDButton, [MultiplyElectricChargeFactor, DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfElectricDiffusionForSelectedSpace(TypesOfLookingForParticlesInProximity::InChosenSectorOfSimulationSpace, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, 19, 0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
                });

                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION NOT IN BOUNDS FOR SELECTED BIG PART OF CELL", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfDiffusionForBigPartOfCellSpace(false, CellEngineConfigDataObject.SizeOfBigPartOfTheCellMultiplyFactor, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });
                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION NOT IN BOUNDS FOR WHOLE CELL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfDiffusionForWholeCellSpace(false, 0, 0, 0, 32, 32, 32, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });

                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION IN BOUNDS FOR SELECTED BIG PART OF CELL", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [DrawSpaceSizesXYZ, DrawSpaceStartXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfDiffusionForBigPartOfCellSpace(true, CellEngineConfigDataObject.SizeOfBigPartOfTheCellMultiplyFactor, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });
                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION IN BOUNDS FOR WHOLE CELL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfDiffusionForWholeCellSpace(true, 0, 0, 0, 32, 32, 32, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });
            }
        }
        CATCH("modification of diffusion operations voxel simulation space parameters menu")
    }

    static void MenuChemicalReactions(ImGuiWindowFlags WindowFlags, const bool& ModifiableWindow, const int DrawSpaceStartXYZ[], const int DrawSpaceSizesXYZ[], bool* OpenMenuChemicalReactionsWindow)
    {
        try
        {
            UnsignedInt StringLength = 90;
            UnsignedInt PrefixStringLength = 0;
            if (ModifiableWindow == false)
                ImGui::Begin("Chemical ChemicalReactions Menu", nullptr, WindowFlags);
            else
                ImGui::Begin("Chemical ChemicalReactions Menu");

            ImGui::Text("CHOOSE REACTION");

            if (ImGui::Button(AlignString("ONLY FIND PARTICLES NR = 0", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-STD ONLY WITH SEQ NR = 1101", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(1101, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 1 SEQ NR = 10", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 2 SEQ SHIFT 3 10 NR = 40", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(40, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 2 SEQ SHIFT 7 3 NR = 41", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(41, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 2 SEQ SHIFT 3 10 NR = 42", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(42, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 1 SEQ NR = 20", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(20, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 1 ANY NR = 30", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(30, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 ANY EQU SAME NR = 80", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(80, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 SEQ COMPLEMENT NR = 60", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(60, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 SEQ COMPLEMENT NR = 61", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(61, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 ANY COMPLEMENT NR = 70", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(70, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT CRISPER 1 NR = 100", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(100, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT CRISPER 2 NR = 110", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(110, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-POLYMERASE DNA START SEQ SPACE NR = 150", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(150, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-POLYMERASE DNA CONTINUE SPACE NR = 160", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(160, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-RIBOSOME RNA START SEQ SPACE NR = 170", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(170, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-RIBOSOME RNA CONTINUE SPACE NR = 180", StringLength, PrefixStringLength).c_str()) == true)
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(180, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("CLOSE", StringLength, PrefixStringLength).c_str()) == true)
                *OpenMenuChemicalReactionsWindow = false;

            ImGui::End();
        }
        CATCH("executing menu chemical reactions");
    }

    static void ReactionsVoxelSimulationSpaceParametersMenu(int DrawSpaceStartXYZ[3], int DrawSpaceStepsXYZ[3], int DrawSpaceSizesXYZ[3], const UnsignedInt StringLength, const ImGuiWindowFlags WindowFlags, const bool ModifiableWindow)
    {
        try
        {
            if (ImGui::CollapsingHeader("REACTIONS", ImGuiTreeNodeFlags_DefaultOpen))
            {
                int IDButton = 1;
                float Nothing;

                // ColorButton(AlignString("MAKE ONE STEP OF RANDOM REACTIONS FOR RANGE OF PARTICLES", StringLength).c_str(), Nothing, 0, 0, 0, 5, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                // {
                //     CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(10, 0);
                // });
                // ColorButton(AlignString("MAKE ONE RANDOM REACTION FOR ONE CHOSEN PARTICLE", StringLength).c_str(), Nothing, 0, 0, 0, 5, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                // {
                //     CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneStepOfRandomReactionsForOneParticleFromRangeOfParticles(10, 0, 4);
                // });

                ColorButton(AlignString("MAKE ONE RANDOM REACTION FOR SELECTED SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateOneRandomReactionForSelectedSpace(DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], true);
                });

                ColorButton(AlignString("MAKE ONE STEP OF RANDOM REACTIONS FOR SELECTED BIG PART OF CELL", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [DrawSpaceStartXYZ, DrawSpaceSizesXYZ](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(CellEngineConfigDataObject.SizeOfBigPartOfTheCellMultiplyFactor, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2], CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });

                ColorButton(AlignString("MAKE ONE STEP OF RANDOM REACTIONS FOR WHOLE CELL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfOneRandomReactionForWholeCellSpace(0, 0, 0, 32, 32, 32, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });
                ColorButton(AlignString("MAKE ONE STEP OF CHOSEN REACTIONS FOR WHOLE CELL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfOneChosenReactionForWholeCellSpace(10, 0, 0, 0, 32, 32, 32, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });

                ImGui::Text("");

                if (ImGui::Button(AlignString("PRINT CHEMICAL REACTIONS", StringLength).c_str()) == true)
                    ChemicalReactionsManagerObject.PrintChemicalReactions();

                static bool OpenMenuChemicalReactionsWindow = false;
                if (ImGui::Button(AlignString("CHOOSE CHEMICAL REACTION", StringLength).c_str()) == true)
                    OpenMenuChemicalReactionsWindow = true;

                if (OpenMenuChemicalReactionsWindow == true)
                    MenuChemicalReactions(WindowFlags, ModifiableWindow, DrawSpaceStartXYZ, DrawSpaceSizesXYZ, &OpenMenuChemicalReactionsWindow);
            }
        }
        CATCH("modification of reactions operations voxel simulation space parameters menu")
    }

    static void SimulationsVoxelSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("SIMULATIONS VOXEL SPACE", ImGuiTreeNodeFlags_DefaultOpen))
            {
                int IDButton = 1;
                float Nothing;

                ImGui::Text("");
                ImGui::Text("Number Of Simulation Steps Outside");
                ImGui::DragInt("Num Of Steps Outside", &CellEngineConfigDataObject.NumberOfStepsInSimulationOutside, 1, 1, 1000, "%d", ImGuiSliderFlags_AlwaysClamp);
                ImGui::Text("");
                ImGui::Text("Number Of Simulation Steps Inside");
                ImGui::DragInt("Num Of Steps Inside", &CellEngineConfigDataObject.NumberOfStepsInSimulationInside, 1, 1, 1000, "%d", ImGuiSliderFlags_AlwaysClamp);
                ImGui::Text("");
                ImGui::Text("Type Of Simulation");
                int TypeOfSimulation = static_cast<int>(CellEngineConfigDataObject.TypeOfSimulation);
                ImGui::RadioButton("BothReactionsAndDiffusion", &TypeOfSimulation, 1);
                ImGui::RadioButton("Only Reactions ", &TypeOfSimulation, 2);
                ImGui::RadioButton("Only Diffusion", &TypeOfSimulation, 3);
                CellEngineConfigDataObject.TypeOfSimulation = static_cast<CellEngineConfigData::TypesOfSimulation>(TypeOfSimulation);
                ImGui::Text("");
                ImGui::Checkbox("Use Mutex Between Main Screen Thread and Menu Threads", &CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads);
                ImGui::Text("");

                ColorButton(AlignString("START N STEPS OF SIMULATION FOR WHOLE CELL SPACE IN THREADS", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->FirstSendParticlesForThreads(false, true);
                });
                ColorButton(AlignString("MAKE N STEPS OF SIMULATION FOR WHOLE CELL SPACE IN THREADS", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfSimulationForWholeCellSpaceInThreads(CellEngineConfigDataObject.NumberOfStepsInSimulationOutside, CellEngineConfigDataObject.NumberOfStepsInSimulationInside);
                });
                ColorButton(AlignString("GATHER PARTICLES FROM THREADS AFTER SIMULATION", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GatherParticlesFromThreadsToParticlesInMainThread();
                });

                ImGui::Text("");

                ColorButton(AlignString("MAKE FULL N STEPS OF SIMULATION FOR WHOLE CELL SPACE IN THREADS", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(CellEngineConfigDataObject.NumberOfStepsInSimulationOutside, CellEngineConfigDataObject.NumberOfStepsInSimulationInside, true);
                });

                ImGui::Text("");

                ColorButton(AlignString("CHECK PARTICLES CENTERS", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->CheckParticlesCenters(false);
                });
                ColorButton(AlignString("CHECK PARTICLES CANCELLED STILL IN VOXEL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->CheckCancelledParticlesIndexes();
                });
                ColorButton(AlignString("CHECK IF PARTICLES FORMER EXISTED STILL IN VOXEL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->CheckFormerExistedParticlesIndexes();
                });
            }
        }
        CATCH("modification of simulations operations voxel simulation space parameters menu")
    }

    static void DNARandomGeneratorVoxelSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("RANDOM DNA GENERATOR", ImGuiTreeNodeFlags_DefaultOpen))
            {
                int IDButton = 1;
                float Nothing;

                ColorButton(AlignString("GENERATE RANDOM DNA 1 RANDOM TURN", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomDNAInWholeCell1or2RandomTurn(0, CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos + 3, CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos, 1, 1, 1, 1, 1, 1, 1, 1, false);
                });
                ColorButton(AlignString("GENERATE RANDOM DNA 2 RANDOM TURN", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomDNAInWholeCell1or2RandomTurn(0, CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos + 3, CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos, 1, 1, 1, 1, 1, 1, 1, 1, true);
                });
                ColorButton(AlignString("GENERATE RANDOM DNA 1 VERTICAL", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomDNAInWholeCell1or2Vertical(0, CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos + 3, CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos, 1, 1, 1, 1, 1, 1, 1, 1, false);
                });
                ColorButton(AlignString("GENERATE RANDOM DNA 2 VERTICAL", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float& VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomDNAInWholeCell1or2Vertical(0, CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos + 3, CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos, 2, 2, 2, 2, 2, 2, 2, 2, true);
                });


                if (ImGui::Button(AlignString("TRUE RANDOM GENERATOR SEED DEVICE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RandomGeneratorSetSeedByRandomDevice();
                if (ImGui::Button(AlignString("TRUE RANDOM GENERATOR SEED TIME", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->RandomGeneratorSetSeedByTime();
                if (ImGui::Button(AlignString("SAVE GENOME TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SaveGenomeDataToFile(2);
                if (ImGui::Button(AlignString("READ GENOME DATA FROM FILE", StringLength).c_str()) == true)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeDataFromFile(CellEngineConfigDataObject.DNAPaired);
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ReadGenomeSequenceFromFile(CellEngineConfigDataObject.DNAPaired);
                }
                if (ImGui::Button(AlignString("TEST GENOME DATA FROM FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->TestGeneratedGenomeCorrectness(2);

                if (ImGui::Button(AlignString("TEST GENOME PROMOTERS FINDER ALGORITHMS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->TestDifferentKindsOfPromotersFindingsAndTerminatorsFindingsAlgorithms();
            }
        }
        CATCH("modification of random dna generator operations voxel simulation space parameters menu")
    }

    static void RandomParticlesGeneratorVoxelSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("RANDOM PARTICLES GENERATOR"))
            {
                if (ImGui::Button(AlignString("SHOW ALL PARTICLES KINDS", StringLength).c_str()) == true)
                    ParticlesKindsManagerObject.PrintAllParticleKinds();
                if (ImGui::Button(AlignString("SHOW NUMBER OF PARTICLES TYPES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->PrintNumberOfParticlesForAllMainTypesOfParticles();
                if (ImGui::Button(AlignString("UPDATE SEQUENCE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->UpdateSequence(ParticlesTypes::mRNA);
                if (ImGui::Button(AlignString("CLEAR VOXEL SPACE AND PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ClearVoxelSpaceAndParticles();
                if (ImGui::Button(AlignString("CLEAR DNA PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->EraseAllDNAParticles();

                if (ImGui::Button(AlignString("GENERATE ALL RANDOM PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateAllRealRandomParticles();

                const UnsignedInt Radius1 = CellEngineConfigDataObject.Radius1ForGenerationOfParticles;
                const UnsignedInt Radius1Size = CellEngineConfigDataObject.Radius1SizeForGenerationOfParticles;
                const UnsignedInt Radius2 = CellEngineConfigDataObject.Radius2ForGenerationOfParticles;
                const UnsignedInt Radius2Size = CellEngineConfigDataObject.Radius2SizeForGenerationOfParticles;

                if (ImGui::Button(AlignString("SHOW DATA PARAMETERS", StringLength).c_str()) == true)
                    LoggersManagerObject.Log(STREAM("R = " << Radius1 << " " << Radius1Size << " " << Radius2 << " " << Radius2Size));

                if (ImGui::Button(AlignString("GENERATE RANDOM RIBOSOMES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Ribosome, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM RNA POLYMERASE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymerase, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM DNA POLYMERASE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::DNAPolymerase, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM MEMBRANE PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::MembraneProtein, false, Radius2, Radius2Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM RIBOSOMES PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::RibosomeProtein, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM DNA POLYMERASE P PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::PolymeraseProtein, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM RNA POLYMERASE P PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymeraseProtein, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM PROTEIN FRAC PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::ProteinFrac, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM OTHER PROTEIN PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::OtherProtein, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_uncharged", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_charged", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM mRNA PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::mRNA, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM rRNA PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::rRNA, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_uncharged M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, true, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_charged M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, true, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM mRNA PARTICLES M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::mRNA, true, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM rRNA PARTICLES M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::rRNA, true, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM BASIC PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Basic, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM LIPID PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Lipid, false, Radius2, Radius2Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM OTHER PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Other, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("SHOW POLYMERASE PARTICLES TYPES DATA", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ShowParticlesKindsData(ParticlesTypes::RNAPolymerase);
            }
        }
        CATCH("modification of often operations voxel simulation space parameters menu")
    }

    static void SavingAndReadingParticlesFromDataFileVoxelSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("SAVING AND READING PARTICLES TO AND FROM FILE", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Button(AlignString("SAVE PARTICLES DATA TO BINARY FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->SaveDataToFile();
                if (ImGui::Button(AlignString("READ PARTICLES DATA FROM BINARY FILE", StringLength).c_str()) == true)
                {
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->ClearWholeVoxelSpace();
                    CellEngineDataFileObjectPointer->ReadDataFromFile(false, false, CellEngineConfigData::TypesOfFileToRead::BinaryFile);
                }
            }
        }
        CATCH("modification of saving and reading particles from data file operations voxel simulation space parameters menu")
    }

    static void StatisticsVoxelSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("STATISTICS OF SIMULATION", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Button(AlignString("ZERO STATISTICS CONTAINERS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SetMakeSimulationStepNumberZero();
                if (ImGui::Button(AlignString("INCREMENT STATISTICS LEVEL", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SetIncSimulationStepNumber();
                if (ImGui::Button(AlignString("SAVE PARTICLES STATISTICS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SaveParticlesStatisticsOnce();
                if (ImGui::Button(AlignString("JOIN REACTIONS STATISTICS FROM THREADS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->JoinStatisticsFromThreads(CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SavedReactionsMap, CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SimulationStepNumber);
                if (ImGui::Button(AlignString("SAVE REACTIONS STATISTICS TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SaveReactionsStatisticsToFile();
                if (ImGui::Button(AlignString("SAVE HISTOGRAM OF PARTICLES TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SaveHistogramOfParticlesStatisticsToFile();
                if (ImGui::Button(AlignString("SAVE NUMBER OF PARTICLES TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SaveNumberOfParticlesStatisticsToFile();
            }
        }
        CATCH("modification of statistics operations voxel simulation space parameters menu")
    }

    static void CombinatoricsMenu()
    {
        try
        {
            if (ImGui::CollapsingHeader("COMBINATORICS", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Button("DEMONSTRATION") == true)
                    ShowDemoTest();
            }
        }
        CATCH("modification of combinatorics operations voxel simulation space parameters menu")
    }

    static void VoxelSimulationSpaceVisibility(const ImGuiWindowFlags WindowFlags, const bool ModifiableWindow)
    {
        try
        {
            if (CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer != nullptr)
            {
                conditional_lock_guard<recursive_mutex> LockGuardCond(CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads, &CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderMenuAndVoxelSimulationSpaceMutexObject);

                const auto CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer = dynamic_cast<CellEngineOpenGLVisualiserOfVoxelSimulationSpace*>(CellEngineOpenGLVisualiserPointer.get());

                const auto StartPos = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetStartPositions();
                static int DrawSpaceStartXYZ[3] = { static_cast<int>(get<0>(StartPos)), static_cast<int>(get<1>(StartPos)), static_cast<int>(get<2>(StartPos)) };
                ImGui::DragInt3("StartX StartY StartZ", DrawSpaceStartXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                const auto Steps = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetSteps();
                static int DrawSpaceStepsXYZ[3] = { static_cast<int>(get<0>(Steps)), static_cast<int>(get<1>(Steps)), static_cast<int>(get<2>(Steps)) };
                ImGui::DragInt3("StepX  StepY  StepZ", DrawSpaceStepsXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                const auto Sizes = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetSizes();
                static int DrawSpaceSizesXYZ[3] = { static_cast<int>(get<0>(Sizes)), static_cast<int>(get<1>(Sizes)), static_cast<int>(get<2>(Sizes)) };
                ImGui::DragInt3("SizeX  SizeY  SizeZ", DrawSpaceSizesXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                constexpr UnsignedInt StringLength = 70;


                VoxelSimulationSpaceParametersMenu(CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer, DrawSpaceStartXYZ, DrawSpaceStepsXYZ, DrawSpaceSizesXYZ, StringLength);

                OftenOperationsVoxelSimulationSpaceParametersMenu(DrawSpaceStartXYZ, DrawSpaceStepsXYZ, DrawSpaceSizesXYZ, StringLength);

                DiffusionVoxelSimulationSpaceParametersMenu(DrawSpaceStartXYZ, DrawSpaceStepsXYZ, DrawSpaceSizesXYZ, StringLength);

                ReactionsVoxelSimulationSpaceParametersMenu(DrawSpaceStartXYZ, DrawSpaceStepsXYZ, DrawSpaceSizesXYZ, StringLength, WindowFlags, ModifiableWindow);

                SimulationsVoxelSimulationSpaceParametersMenu(StringLength);

                DNARandomGeneratorVoxelSimulationSpaceParametersMenu(StringLength);

                RandomParticlesGeneratorVoxelSimulationSpaceParametersMenu(StringLength);

                SavingAndReadingParticlesFromDataFileVoxelSimulationSpaceParametersMenu(StringLength);

                StatisticsVoxelSimulationSpaceParametersMenu(StringLength);

                CombinatoricsMenu();
            }
        }
        CATCH("executing voxel simulation space visibility menu")
    }

    static void DiffusionFullAtomSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            int IDButton = 1;
            float Nothing;

            if (ImGui::CollapsingHeader("WORK", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ColorButton(AlignString("START WORK", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});
                ColorButton(AlignString("STOP WORK", StringLength).c_str(), Nothing, 0, 0, 0, 0, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});
            }

            if (ImGui::CollapsingHeader("DIFFUSION", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ColorButton(AlignString("START DIFFUSION", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});
                ColorButton(AlignString("STOP DIFFUSION", StringLength).c_str(), Nothing, 0, 0, 0, 0, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue){});

                constexpr UnsignedInt AdditionalSpaceBoundFactor = 20;
                constexpr double MultiplyElectricChargeFactor = 100;

                ColorButton(AlignString("MAKE ONE STEP OF DIFFUSION IN BOUNDS FOR WHOLE CELL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 3, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateNStepsOfDiffusionForWholeCellSpace(true, 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.NumberOfParticlesSectorsInX, CellEngineConfigDataObject.NumberOfParticlesSectorsInY, CellEngineConfigDataObject.NumberOfParticlesSectorsInZ, CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });
            }
        }
        CATCH("modification of diffusion operations atom simulation space parameters menu")
    }

    static void MenuFullAtomChemicalReactions(ImGuiWindowFlags WindowFlags, const bool& ModifiableWindow, const int DrawSpaceStartXYZ[], const int DrawSpaceSizesXYZ[], bool* OpenMenuChemicalReactionsWindow)
    {
        try
        {
            UnsignedInt StringLength = 90;
            UnsignedInt PrefixStringLength = 0;
            if (ModifiableWindow == false)
                ImGui::Begin("Chemical ChemicalReactions Menu", nullptr, WindowFlags);
            else
                ImGui::Begin("Chemical ChemicalReactions Menu");

            ImGui::Text("CHOOSE REACTION");

            // if (ImGui::Button(AlignString("ONLY FIND PARTICLES NR = 0", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(0, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-STD ONLY WITH SEQ NR = 1101", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(1101, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 1 SEQ NR = 10", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(10, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 2 SEQ SHIFT 3 10 NR = 40", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(40, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 2 SEQ SHIFT 7 3 NR = 41", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(41, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT 2 SEQ SHIFT 3 10 NR = 42", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(42, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 1 SEQ NR = 20", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(20, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 1 ANY NR = 30", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(30, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 ANY EQU SAME NR = 80", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(80, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 SEQ COMPLEMENT NR = 60", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(60, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 SEQ COMPLEMENT NR = 61", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(61, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-LINK 2 ANY COMPLEMENT NR = 70", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(70, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT CRISPER 1 NR = 100", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(100, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-CUT CRISPER 2 NR = 110", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(110, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-POLYMERASE DNA START SEQ SPACE NR = 150", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(150, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-POLYMERASE DNA CONTINUE SPACE NR = 160", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(160, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            //
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-RIBOSOME RNA START SEQ SPACE NR = 170", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(170, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);
            // if (ImGui::Button(AlignString("DO CHOSEN REACTION IN SELECTED SPACE-RIBOSOME RNA CONTINUE SPACE NR = 180", StringLength, PrefixStringLength).c_str()) == true)
            //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateOneChosenReactionForSelectedSpace(180, DrawSpaceStartXYZ[0], DrawSpaceStartXYZ[1], DrawSpaceStartXYZ[2], DrawSpaceSizesXYZ[0], DrawSpaceSizesXYZ[1], DrawSpaceSizesXYZ[2]);

            if (ImGui::Button(AlignString("CLOSE", StringLength, PrefixStringLength).c_str()) == true)
                *OpenMenuChemicalReactionsWindow = false;

            ImGui::End();
        }
        CATCH("executing menu full atom chemical reactions");
    }

    static void ReactionsFullAtomSimulationSpaceParametersMenu(const UnsignedInt StringLength, const ImGuiWindowFlags WindowFlags, const bool ModifiableWindow)
    {
        try
        {
            if (ImGui::CollapsingHeader("REACTIONS", ImGuiTreeNodeFlags_DefaultOpen))
            {
                int IDButton = 1;
                float Nothing;

                ColorButton(AlignString("MAKE ONE STEP OF RANDOM REACTIONS FOR WHOLE CELL SPACE", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateNStepsOfOneRandomReactionForWholeCellSpace(0, 0, 0, 32, 32, 32, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfStepsInSimulationOutside);
                });

                ImGui::Text("");

                if (ImGui::Button(AlignString("PRINT CHEMICAL REACTIONS", StringLength).c_str()) == true)
                    ChemicalReactionsManagerObject.PrintChemicalReactions();

                static bool OpenMenuChemicalReactionsWindow = false;
                if (ImGui::Button(AlignString("CHOOSE CHEMICAL REACTION", StringLength).c_str()) == true)
                    OpenMenuChemicalReactionsWindow = true;

                // if (OpenMenuChemicalReactionsWindow == true)
                //     MenuFullAtomChemicalReactions(WindowFlags, ModifiableWindow, DrawSpaceStartXYZ, DrawSpaceSizesXYZ, &OpenMenuChemicalReactionsWindow);
            }
        }
        CATCH("modification of reactions operations full atom simulation space parameters menu")
    }

    static void SimulationsFullAtomSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("SIMULATIONS FULL ATOM", ImGuiTreeNodeFlags_DefaultOpen))
            {
                int IDButton = 1;
                float Nothing;

                ImGui::Text("");
                ImGui::Text("Number Of Simulation Steps Outside");
                ImGui::DragInt("Num Of Steps Outside", &CellEngineConfigDataObject.NumberOfStepsInSimulationOutside, 1, 1, 1000, "%d", ImGuiSliderFlags_AlwaysClamp);
                ImGui::Text("");
                ImGui::Text("Number Of Simulation Steps Inside");
                ImGui::DragInt("Num Of Steps Inside", &CellEngineConfigDataObject.NumberOfStepsInSimulationInside, 1, 1, 1000, "%d", ImGuiSliderFlags_AlwaysClamp);
                ImGui::Text("");
                ImGui::Text("Type Of Simulation");
                int TypeOfSimulation = static_cast<int>(CellEngineConfigDataObject.TypeOfSimulation);
                ImGui::RadioButton("BothReactionsAndDiffusion", &TypeOfSimulation, 1);
                ImGui::RadioButton("Only Reactions ", &TypeOfSimulation, 2);
                ImGui::RadioButton("Only Diffusion", &TypeOfSimulation, 3);
                CellEngineConfigDataObject.TypeOfSimulation = static_cast<CellEngineConfigData::TypesOfSimulation>(TypeOfSimulation);
                ImGui::Text("");
                ImGui::Checkbox("Use Mutex Between Main Screen Thread and Menu Threads", &CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads);
                ImGui::Text("");

                ColorButton(AlignString("MAKE N STEPS OF SIMULATION FOR WHOLE CELL SPACE IN THREADS", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateNStepsOfSimulationForWholeCellSpaceInThreads(CellEngineConfigDataObject.NumberOfStepsInSimulationOutside, CellEngineConfigDataObject.NumberOfStepsInSimulationInside);
                });

                ColorButton(AlignString("CHECK PARTICLES CENTERS", StringLength).c_str(), Nothing, 0, 0, 0, 6, IDButton, [](float &VariableToChange, const float Step, const float MinValue, const float MaxValue)
                {
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->CheckParticlesCenters(false);
                });
            }
        }
        CATCH("modification of simulations operations full atom simulation space parameters menu")
    }

static void RandomParticlesGeneratorFullAtomSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("RANDOM PARTICLES GENERATOR"))
            {
                if (ImGui::Button(AlignString("SHOW ALL PARTICLES KINDS", StringLength).c_str()) == true)
                    ParticlesKindsManagerObject.PrintAllParticleKinds();
                if (ImGui::Button(AlignString("SHOW NUMBER OF PARTICLES TYPES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->PrintNumberOfParticlesForAllMainTypesOfParticles();
                if (ImGui::Button(AlignString("UPDATE SEQUENCE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->UpdateSequence(ParticlesTypes::mRNA);
                // if (ImGui::Button(AlignString("CLEAR VOXEL SPACE AND PARTICLES", StringLength).c_str()) == true)
                //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->ClearVoxelSpaceAndParticles();
                // if (ImGui::Button(AlignString("CLEAR DNA PARTICLES", StringLength).c_str()) == true)
                //     CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->EraseAllDNAParticles();

                if (ImGui::Button(AlignString("GENERATE ALL RANDOM PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->GenerateAllRealRandomParticles();

                const UnsignedInt Radius1 = CellEngineConfigDataObject.Radius1ForGenerationOfParticles;
                const UnsignedInt Radius1Size = CellEngineConfigDataObject.Radius1SizeForGenerationOfParticles;
                const UnsignedInt Radius2 = CellEngineConfigDataObject.Radius2ForGenerationOfParticles;
                const UnsignedInt Radius2Size = CellEngineConfigDataObject.Radius2SizeForGenerationOfParticles;

                if (ImGui::Button(AlignString("SHOW DATA PARAMETERS", StringLength).c_str()) == true)
                    LoggersManagerObject.Log(STREAM("R = " << Radius1 << " " << Radius1Size << " " << Radius2 << " " << Radius2Size));

                if (ImGui::Button(AlignString("GENERATE RANDOM RIBOSOMES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Ribosome, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM RNA POLYMERASE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymerase, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM DNA POLYMERASE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::DNAPolymerase, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM MEMBRANE PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::MembraneProtein, false, Radius2, Radius2Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM RIBOSOMES PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::RibosomeProtein, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM DNA POLYMERASE P PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::PolymeraseProtein, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM RNA POLYMERASE P PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymeraseProtein, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM PROTEIN FRAC PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::ProteinFrac, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM OTHER PROTEIN PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::OtherProtein, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_uncharged", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_charged", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM mRNA PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::mRNA, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM rRNA PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::rRNA, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_uncharged M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, true, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM tRNA PARTICLES_charged M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, true, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM mRNA PARTICLES M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::mRNA, true, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM rRNA PARTICLES M", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::rRNA, true, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("GENERATE RANDOM BASIC PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Basic, false, Radius1, Radius1Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM LIPID PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Lipid, false, Radius2, Radius2Size);
                if (ImGui::Button(AlignString("GENERATE RANDOM OTHER PARTICLES", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->InsertNewRandomParticlesForType(ParticlesTypes::Other, false, Radius1, Radius1Size);

                if (ImGui::Button(AlignString("SHOW POLYMERASE PARTICLES TYPES DATA", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->ShowParticlesKindsData(ParticlesTypes::RNAPolymerase);
            }
        }
        CATCH("modification of often operations full atom simulation space parameters menu")
    }

    static void SavingAndReadingParticlesFromDataFileFullAtomSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("SAVING AND READING PARTICLES TO AND FROM FILE", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Button(AlignString("SAVE PARTICLES DATA TO BINARY FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->SaveDataToFile();
                if (ImGui::Button(AlignString("READ PARTICLES DATA FROM BINARY FILE", StringLength).c_str()) == true)
                {
                    //CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->ClearWholeVoxelSpace();
                    CellEngineDataFileObjectPointer->ReadDataFromFile(false, false, CellEngineConfigData::TypesOfFileToRead::BinaryFile);
                }
            }
        }
        CATCH("modification of saving and reading particles from data file operations voxel simulation space parameters menu")
    }

    static void StatisticsFullAtomSimulationSpaceParametersMenu(const UnsignedInt StringLength)
    {
        try
        {
            if (ImGui::CollapsingHeader("STATISTICS OF SIMULATION", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Button(AlignString("ZERO STATISTICS CONTAINERS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SetMakeSimulationStepNumberZero();
                if (ImGui::Button(AlignString("INCREMENT STATISTICS LEVEL", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SetIncSimulationStepNumber();
                if (ImGui::Button(AlignString("SAVE PARTICLES STATISTICS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SaveParticlesStatisticsOnce();
                if (ImGui::Button(AlignString("JOIN REACTIONS STATISTICS FROM THREADS", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->JoinStatisticsFromThreads(CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMap, CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SimulationStepNumber);
                if (ImGui::Button(AlignString("SAVE REACTIONS STATISTICS TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SaveReactionsStatisticsToFile();
                if (ImGui::Button(AlignString("SAVE HISTOGRAM OF PARTICLES TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SaveHistogramOfParticlesStatisticsToFile();
                if (ImGui::Button(AlignString("SAVE NUMBER OF PARTICLES TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SaveNumberOfParticlesStatisticsToFile();

                if (ImGui::Button(AlignString("SAVE NUMBER OF PARTICLES IN SECTORS TO FILE", StringLength).c_str()) == true)
                    CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->WriteNumberOfParticlesInEachSectorToFile();

                if (ImGui::Button(AlignString("SAVE NUMBER OF PARTICLES KINDS WITHOUT ATOMS", StringLength).c_str()) == true)
                    CellEngineFullAtomSimulationSpace::WriteNumberOfParticlesKindsWithoutAtoms();
            }
        }
        CATCH("modification of statistics operations full atom simulation space parameters menu")
    }

    static void FullAtomSimulationSpaceVisibility(const ImGuiWindowFlags WindowFlags, const bool ModifiableWindow)
    {
        try
        {
            if (CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer != nullptr)
            {
                lock_guard<mutex> LockGuard{ CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderMenuAndFullAtomSimulationSpaceMutexObject };

                const auto CellEngineOpenGLFullAtomSimulationSpaceVisualiserObjectPointer = dynamic_cast<CellEngineOpenGLVisualiserOfVoxelSimulationSpace*>(CellEngineOpenGLVisualiserPointer.get());

                // const auto StartPos = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetStartPositions();
                // static int DrawSpaceStartXYZ[3] = { static_cast<int>(get<0>(StartPos)), static_cast<int>(get<1>(StartPos)), static_cast<int>(get<2>(StartPos)) };
                // ImGui::DragInt3("StartX StartY StartZ", DrawSpaceStartXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);
                //
                // const auto Steps = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetSteps();
                // static int DrawSpaceStepsXYZ[3] = { static_cast<int>(get<0>(Steps)), static_cast<int>(get<1>(Steps)), static_cast<int>(get<2>(Steps)) };
                // ImGui::DragInt3("StepX  StepY  StepZ", DrawSpaceStepsXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);
                //
                // const auto Sizes = CellEngineOpenGLVoxelSimulationSpaceVisualiserObjectPointer->GetSizes();
                // static int DrawSpaceSizesXYZ[3] = { static_cast<int>(get<0>(Sizes)), static_cast<int>(get<1>(Sizes)), static_cast<int>(get<2>(Sizes)) };
                // ImGui::DragInt3("SizeX  SizeY  SizeZ", DrawSpaceSizesXYZ, 1, 0, static_cast<int>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension), "%d", ImGuiSliderFlags_AlwaysClamp);

                constexpr UnsignedInt StringLength = 70;


                DiffusionFullAtomSimulationSpaceParametersMenu(StringLength);

                ReactionsFullAtomSimulationSpaceParametersMenu(StringLength, WindowFlags, ModifiableWindow);

                SimulationsFullAtomSimulationSpaceParametersMenu(StringLength);

                // DNARandomGeneratorVoxelSimulationSpaceParametersMenu(StringLength);

                RandomParticlesGeneratorFullAtomSimulationSpaceParametersMenu(StringLength);

                SavingAndReadingParticlesFromDataFileFullAtomSimulationSpaceParametersMenu(StringLength);

                StatisticsFullAtomSimulationSpaceParametersMenu(StringLength);
            }
        }
        CATCH("executing full atom simulation space visibility menu")
    }

    static void VisibilityParametersOfParticles()
    {
        try
        {
            DensityOfDrawnAtomsMenu();

            DetailedVisibilityParametersOfParticles();

            TypesOfVisibilityMenu();

            AtomKindsMenu();

            TypeOfGeneratedColorsMenu();

            ChangeNumberOfParticlesMenu();
        }
        CATCH("drawing visibility parameters of particles")
    }

    static void MenuWindow2(const ImGuiWindowFlags WindowFlags, const bool ModifiableWindow)
    {
        try
        {
            if (ModifiableWindow == false)
                ImGui::Begin("Details Of Cell Drawing", nullptr, WindowFlags);
            else
                ImGui::Begin("Details Of Cell Drawing");

            if (ImGui::CollapsingHeader("Visibility parameters of particles", ImGuiTreeNodeFlags_DefaultOpen))
                VisibilityParametersOfParticles();

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                if (ImGui::CollapsingHeader("Visibility of Voxel Simulation Space", ImGuiTreeNodeFlags_DefaultOpen))
                    VoxelSimulationSpaceVisibility(WindowFlags, ModifiableWindow);

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
                if (ImGui::CollapsingHeader("Visibility of Voxel Simulation Space", ImGuiTreeNodeFlags_DefaultOpen))
                    FullAtomSimulationSpaceVisibility(WindowFlags, ModifiableWindow);

            ImGui::End();
        }
        CATCH("executing menu window 2");
    }

    static void ImGuiMenuGLFWMainLoop(GLFWwindow* ImGuiMenuWindow)
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
    static inline unique_ptr<CellEngineOpenGLVisualiser> CellEngineOpenGLVisualiserPointer;

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wreturn-type"
    static unique_ptr<CellEngineOpenGLVisualiser> CreateCellEngineOpenGLVisualiserObject(const CellEngineConfigData::TypesOfSpace TypeOfSpace)
    {
        switch (TypeOfSpace)
        {
            case CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace : return make_unique<CellEngineOpenGLVisualiserOfFullAtomSimulationSpace>();
            case CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace : return make_unique<CellEngineOpenGLVisualiserOfVoxelSimulationSpace>();
            default : break;
        }
        return nullptr;
    }
    #pragma GCC diagnostic pop

    void CellEngineOpenGLVisualiserThreadFunction(const int XPosWindow, const int YPosWindow, const int WidthWindow, const int HeightWindow)
    {
        try
        {
            CellEngineOpenGLVisualiserPointer = CreateCellEngineOpenGLVisualiserObject(CellEngineConfigDataObject.TypeOfSpace);
            CellEngineOpenGLVisualiserPointer->Run(XPosWindow, YPosWindow, WidthWindow, HeightWindow);
        }
        CATCH("running cell engine opengl visualiser thread function");
    }

public:
    CellEngineImGuiMenu(const int argc, const char** argv)
    {
        try
        {
            ReadInitConfiguration(argc, argv);

            GLFWwindow* ImGuiMenuWindow = PrepareImGuiMenuGLFWData();

            thread CellEngineOpenGLVisualiserThreadObject(&CellEngineImGuiMenu::CellEngineOpenGLVisualiserThreadFunction, this, CellEngineConfigDataObject.XTopMainWindow, CellEngineConfigDataObject.YTopMainWindow, CellEngineConfigDataObject.WidthMainWindow, CellEngineConfigDataObject.HeightMainWindow);

            ImGuiMenuGLFWMainLoop(ImGuiMenuWindow);

            CellEngineOpenGLVisualiserThreadObject.detach();

            ImGuiMenuGLFShutdown(ImGuiMenuWindow);
        }
        CATCH("starting imgui menu and whole cell opengl visualization")
    }
};