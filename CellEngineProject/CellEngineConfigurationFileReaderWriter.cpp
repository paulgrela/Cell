
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <sb7color.h>

#include "StringUtils.h"
#include "ExceptionsMacro.h"
#include "DestinationPlatform.h"

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"

#include "CellEnginePDBDataFile.h"
#include "CellEngineConfigurationFileReaderWriter.h"
#include "CellEngineCIFDataFileReaderOfVoxelSimulationSpace.h"
#include "CellEngineCIFDataFileReaderOfFullAtomSimulationSpace.h"

using namespace std;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
std::unique_ptr<CellEngineDataFile> CreateCellEngineDataFileObject(const string_view& CellStateFileName)
{
    if (string_utils::check_end_str(CellStateFileName, ".pdb") == true)
        return make_unique<CellEnginePDBDataFile>();
    else
        switch (CellEngineConfigDataObject.TypeOfSpace)
        {
            case CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace : return make_unique<CellEngineCIFDataFileReaderOfVoxelSimulationSpace>();
            case CellEngineConfigData::TypesOfSpace::FullAtomSpace : return make_unique<CellEngineCIFDataFileReaderOfFullAtomSimulationSpace>();
            default : break;
        }
}
#pragma GCC diagnostic pop

void CellEngineConfigurationFileReaderWriter::ReadCellConfigurationFile(const char* ConfigFileNameParameter, const UnsignedInt ExecuteCellStateId)
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        ConfigFileName = ConfigFileNameParameter;

        read_xml(ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);
        LoggersManagerObject.Log(STREAM("Reading xml config file finished"));

        LoggersManagerObject.Log(STREAM("ExecuteCellStateId = " << ExecuteCellStateId));

        for (const ptree::value_type& MainConfigPropertyTreeElement : MainConfigPropertyTree.get_child("Settings"))
        {
            if (MainConfigPropertyTreeElement.first == "WindowsParameters")
            {
                #ifdef WINDOWS_PLATFORM
                CellEngineConfigDataObject.XTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMainWindows");
                CellEngineConfigDataObject.YTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMainWindows");
                CellEngineConfigDataObject.WidthMainWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMainWindows");
                CellEngineConfigDataObject.HeightMainWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMainWindows");
                CellEngineConfigDataObject.XTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMenuWindows");
                CellEngineConfigDataObject.YTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMenuWindows");
                CellEngineConfigDataObject.WidthMenuWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMenuWindows");
                CellEngineConfigDataObject.HeightMenuWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMenuWindows");
                CellEngineConfigDataObject.XTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondWindows");
                CellEngineConfigDataObject.YTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("YTopSecondWindows");
                CellEngineConfigDataObject.WidthSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondWindows");
                CellEngineConfigDataObject.HeightSecondWindow = MainConfigPropertyTreeElement.second.get<int>("HeightSecondWindows");
                #endif
                #ifdef UNIX_PLATFORM
                CellEngineConfigDataObject.XTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMainLinux");
                CellEngineConfigDataObject.YTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMainLinux");
                CellEngineConfigDataObject.WidthMainWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMainLinux");
                CellEngineConfigDataObject.HeightMainWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMainLinux");
                CellEngineConfigDataObject.XTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMenuLinux");
                CellEngineConfigDataObject.YTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMenuLinux");
                CellEngineConfigDataObject.WidthMenuWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMenuLinux");
                CellEngineConfigDataObject.HeightMenuWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMenuLinux");
                CellEngineConfigDataObject.XTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondLinux");
                CellEngineConfigDataObject.YTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("YTopSecondLinux");
                CellEngineConfigDataObject.WidthSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecondLinux");
                CellEngineConfigDataObject.HeightSecondWindow = MainConfigPropertyTreeElement.second.get<int>("HeightSecondLinux");
                #endif
            }
            else
            if (MainConfigPropertyTreeElement.first == "Algorithm")
            {
                CellEngineConfigDataObject.MultiThreaded = MainConfigPropertyTreeElement.second.get<bool>("MultiThreaded");
                CellEngineConfigDataObject.SetProcessPriorityHighest = MainConfigPropertyTreeElement.second.get<bool>("SetProcessPriorityHighest");
                CellEngineConfigDataObject.PrintAtomDescriptionOnScreen = MainConfigPropertyTreeElement.second.get<bool>("PrintAtomDescriptionOnScreen");
                CellEngineConfigDataObject.LogParametersOfRenderingToFile = MainConfigPropertyTreeElement.second.get<bool>("LogParametersOfRenderingToFile");
                CellEngineConfigDataObject.RandomColorEngineObject = static_cast<CellEngineConfigData::RandomColorEngineTypes>(MainConfigPropertyTreeElement.second.get<UnsignedInt>("RandomColorEngineTypes"));
            }
            else
            if (MainConfigPropertyTreeElement.first == "Logger")
            {
                CellEngineConfigDataObject.PrintLogToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintToConsole");
                CellEngineConfigDataObject.PrintLogToFiles = MainConfigPropertyTreeElement.second.get<bool>("PrintToFiles");

                CellEngineConfigDataObject.PrintLogLineNumberToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintLineNumberToConsole");
                CellEngineConfigDataObject.PrintLogDateTimeToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintDateTimeToConsole");
                CellEngineConfigDataObject.PrintLogProcessIdToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessIdToConsole");
                CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessPriorityLevelToConsole");
                CellEngineConfigDataObject.PrintLogThreadIdToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintThreadIdToConsole");

                CellEngineConfigDataObject.PrintLogLineNumberToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintLineNumberToFile");
                CellEngineConfigDataObject.PrintLogDateTimeToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintDateTimeToFile");
                CellEngineConfigDataObject.PrintLogProcessIdToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessIdToFile");
                CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessPriorityLevelToFile");
                CellEngineConfigDataObject.PrintLogThreadIdToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintThreadIdToFile");

                CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile = MainConfigPropertyTreeElement.second.get<UnsignedInt>("MaximalNumberOfLinesInOneFile");
            }
            else
            if (MainConfigPropertyTreeElement.first == "CellsStates")
                for (const ptree::value_type& CellStatePropertyTreeElement : MainConfigPropertyTree.get_child("Settings.CellsStates"))
                    if (CellStatePropertyTreeElement.second.get<UnsignedInt>("<xmlattr>.id") == ExecuteCellStateId)
                    {
                        CellEngineConfigDataObject.TypeOfSpace = static_cast<CellEngineConfigData::TypesOfSpace>(CellStatePropertyTreeElement.second.get<UnsignedInt>("TypeOfSpace"));

                        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                        {
                            CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension = CellStatePropertyTreeElement.second.get<UnsignedInt>("NumberOfVoxelsInVoxelSimulationSpaceInEachDimension");
                            CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace = CellStatePropertyTreeElement.second.get<float>("DivisionFactorForVoxelSimulationSpace");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartXPos = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionStartXPos");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartYPos = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionStartYPos");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartZPos = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionStartZPos");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStepX = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionStepX");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStepY = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionStepY");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStepZ = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionStepZ");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionSizeX = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionSizeX");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionSizeY = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionSizeY");
                            CellEngineConfigDataObject.VoxelSimulationSpaceSelectionSizeZ = CellStatePropertyTreeElement.second.get<UnsignedInt>("VoxelSimulationSpaceSelectionSizeZ");
                        }

                        auto CellStateFileName = CellStatePropertyTreeElement.second.get<string>("CellStateFileName");

                        CellEngineDataFileObjectPointer = CreateCellEngineDataFileObject(CellStateFileName);

                        CellEngineConfigDataObject.CellStateFileName = CellStateFileName;

                        CellEngineConfigDataObject.ChosenStructureIndex = CellStatePropertyTreeElement.second.get<UnsignedInt>("ChosenStructureIndex");

                        CellEngineConfigDataObject.SpecularPower = CellStatePropertyTreeElement.second.get<float>("SpecularPower");
                        CellEngineConfigDataObject.SpecularAlbedo = CellStatePropertyTreeElement.second.get<float>("SpecularAlbedo");
                        CellEngineConfigDataObject.MakeColorsTypeObject = static_cast<CellEngineConfigData::MakeColorsType>(CellStatePropertyTreeElement.second.get<UnsignedInt>("MakeColorsType"));

                        CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject = static_cast<CellEngineConfigData::StencilForDrawingObjectsTypes>(CellStatePropertyTreeElement.second.get<UnsignedInt>("StencilForDrawingObjectsTypes"));
                        CellEngineConfigDataObject.NumberOfStencilBufferLoops = CellStatePropertyTreeElement.second.get<UnsignedInt>("NumberOfStencilBufferLoops");

                        CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters = CellStatePropertyTreeElement.second.get<bool>("DrawBondsBetweenParticlesCenters");
                        CellEngineConfigDataObject.DrawBondsBetweenAtoms = CellStatePropertyTreeElement.second.get<bool>("DrawBondsBetweenAtoms");

                        CellEngineConfigDataObject.ShowDetailsInAtomScale = CellStatePropertyTreeElement.second.get<bool>("ShowDetailsInAtomScale");
                        CellEngineConfigDataObject.CheckAtomVisibility = CellStatePropertyTreeElement.second.get<bool>("CheckAtomVisibility");
                        CellEngineConfigDataObject.CutZ = CellStatePropertyTreeElement.second.get<float>("CutZ");
                        CellEngineConfigDataObject.Distance = CellStatePropertyTreeElement.second.get<float>("Distance");
                        CellEngineConfigDataObject.LoadOfAtomsStep = CellStatePropertyTreeElement.second.get<UnsignedInt>("LoadOfAtomsStep");

                        CellEngineConfigDataObject.XLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("XLowToDrawInAtomScale");
                        CellEngineConfigDataObject.XHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("XHighToDrawInAtomScale");
                        CellEngineConfigDataObject.YLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("YLowToDrawInAtomScale");
                        CellEngineConfigDataObject.YHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("YHighToDrawInAtomScale");
                        CellEngineConfigDataObject.ZLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("ZLowToDrawInAtomScale");
                        CellEngineConfigDataObject.ZHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("ZHighToDrawInAtomScale");

                        CellEngineConfigDataObject.SizeOfAtomsDrawingTypesObject = static_cast<CellEngineConfigData::SizeOfAtomsDrawingTypes>(CellStatePropertyTreeElement.second.get<UnsignedInt>("SizeOfAtomsDrawingTypes"));

                        CellEngineConfigDataObject.SizeOfAtomChangeStep = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomChangeStep");
                        CellEngineConfigDataObject.SizeOfAtomX = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomX");
                        CellEngineConfigDataObject.SizeOfAtomY = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomY");
                        CellEngineConfigDataObject.SizeOfAtomZ = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomZ");
                        CellEngineConfigDataObject.RotationAngle1 = CellStatePropertyTreeElement.second.get<float>("RotationAngle1");
                        CellEngineConfigDataObject.RotationAngle2 = CellStatePropertyTreeElement.second.get<float>("RotationAngle2");
                        CellEngineConfigDataObject.RotationAngle3 = CellStatePropertyTreeElement.second.get<float>("RotationAngle3");
                        CellEngineConfigDataObject.ViewPositionX = CellStatePropertyTreeElement.second.get<float>("ViewPositionX");
                        CellEngineConfigDataObject.ViewPositionY = CellStatePropertyTreeElement.second.get<float>("ViewPositionY");
                        CellEngineConfigDataObject.ViewPositionZ = CellStatePropertyTreeElement.second.get<float>("ViewPositionZ");
                        CellEngineConfigDataObject.CameraXPosition = CellStatePropertyTreeElement.second.get<float>("CameraXPosition");
                        CellEngineConfigDataObject.CameraYPosition = CellStatePropertyTreeElement.second.get<float>("CameraYPosition");
                        CellEngineConfigDataObject.CameraZPosition = CellStatePropertyTreeElement.second.get<float>("CameraZPosition");
                        CellEngineConfigDataObject.CameraXMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraXMoveShortStep");
                        CellEngineConfigDataObject.CameraYMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraYMoveShortStep");
                        CellEngineConfigDataObject.CameraZMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraZMoveShortStep");
                        CellEngineConfigDataObject.CameraXMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraXMoveLongStep");
                        CellEngineConfigDataObject.CameraYMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraYMoveLongStep");
                        CellEngineConfigDataObject.CameraZMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraZMoveLongStep");
                        CellEngineConfigDataObject.ViewPositionX = CellStatePropertyTreeElement.second.get<float>("ViewPositionX");
                        CellEngineConfigDataObject.ViewPositionY = CellStatePropertyTreeElement.second.get<float>("ViewPositionY");
                        CellEngineConfigDataObject.ViewPositionZ = CellStatePropertyTreeElement.second.get<float>("ViewPositionZ");
                        CellEngineConfigDataObject.ViewXMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewXMoveShortStep");
                        CellEngineConfigDataObject.ViewYMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewYMoveShortStep");
                        CellEngineConfigDataObject.ViewZMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewZMoveShortStep");
                        CellEngineConfigDataObject.ViewXMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewXMoveLongStep");
                        CellEngineConfigDataObject.ViewYMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewYMoveLongStep");
                        CellEngineConfigDataObject.ViewZMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewZMoveLongStep");

                        CellEngineConfigDataObject.ViewChangeUsingLongStep = CellStatePropertyTreeElement.second.get<bool>("ViewChangeUsingLongStep");
                        CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom = CellStatePropertyTreeElement.second.get<bool>("AutomaticChangeOfSizeOfAtom");
                        CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = CellStatePropertyTreeElement.second.get<bool>("AutomaticChangeOfLoadAtomsStep");

                        CellEngineConfigDataObject.BackgroundColors[1] = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor1")));
                        CellEngineConfigDataObject.BackgroundColors[2] = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor2")));
                        CellEngineConfigDataObject.BackgroundColors[3] = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor3")));
                        CellEngineConfigDataObject.ChosenBackgroundColor = CellStatePropertyTreeElement.second.get<UnsignedInt>("ChosenBackgroundColor");

                        CellEngineConfigDataObject.AtomsKinds.clear();
                        for (const ptree::value_type& AtomKindPropertyTreeElement : CellStatePropertyTreeElement.second.get_child("Atoms"))
                        {
                            GraphicAtomKind AtomKindObject;

                            LoggersManagerObject.Log(STREAM("Atom Kind Name = " << AtomKindPropertyTreeElement.second.get<string>("Name")));

                            AtomKindObject.Color = GetVector3FormVMathVec3(vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(AtomKindPropertyTreeElement.second.get<string>("Color"))));
                            AtomKindObject.ColorVmathVec3 = vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(AtomKindPropertyTreeElement.second.get<string>("Color")));
                            AtomKindObject.SizeX = AtomKindPropertyTreeElement.second.get<float>("SizeX");
                            AtomKindObject.SizeY = AtomKindPropertyTreeElement.second.get<float>("SizeY");
                            AtomKindObject.SizeZ = AtomKindPropertyTreeElement.second.get<float>("SizeZ");
                            AtomKindObject.Name = AtomKindPropertyTreeElement.second.get<string>("Name");

                            CellEngineConfigDataObject.AtomsKinds.emplace_back(AtomKindObject);
                        }

                        CellEngineConfigDataObject.ParticlesKindsXML.clear();
                        for (const ptree::value_type& ParticleKindPropertyTreeElement : CellStatePropertyTreeElement.second.get_child("Particles"))
                        {
                            GraphicParticleKind ParticleKindObject;

                            LoggersManagerObject.Log(STREAM("ParticleKind Name = " << ParticleKindPropertyTreeElement.second.get<string>("Name") << " ID = " << ParticleKindPropertyTreeElement.second.get<UnsignedInt>("<xmlattr>.id") << " COLOR = " << ParticleKindPropertyTreeElement.second.get<string>("Color")));

                            ParticleKindObject.NameFromXML = ParticleKindPropertyTreeElement.second.get<string>("Name");
                            ParticleKindObject.NameFromXML == "DNA" ? CellEngineConfigDataObject.DNAIdentifier = ParticleKindPropertyTreeElement.second.get<UnsignedInt>("<xmlattr>.id") : 0;
                            ParticleKindObject.NameFromXML == "RNA" ? CellEngineConfigDataObject.RNAIdentifier = ParticleKindPropertyTreeElement.second.get<UnsignedInt>("<xmlattr>.id") : 0;
                            ParticleKindObject.ParticleColor = GetVector3FormVMathVec3(vmath::FromVec4ToVec3(sb7::GetColorVec4FromColorName(ParticleKindPropertyTreeElement.second.get<string>("Color"))));
                            ParticleKindObject.Visible = ParticleKindPropertyTreeElement.second.get<bool>("Visible");
                            ParticleKindObject.SizeX = ParticleKindPropertyTreeElement.second.get<float>("SizeX");
                            ParticleKindObject.SizeY = ParticleKindPropertyTreeElement.second.get<float>("SizeY");
                            ParticleKindObject.SizeZ = ParticleKindPropertyTreeElement.second.get<float>("SizeZ");
                            CellEngineConfigDataObject.ParticlesKindsXML[ParticleKindPropertyTreeElement.second.get<UnsignedInt>("<xmlattr>.id")] = ParticleKindObject;
                        }
                    }
        }
    }
    CATCH("cell print configuration constructor")
}

void CellEngineConfigurationFileReaderWriter::SaveTestStatisticsToFile(const UnsignedInt ExecuteCellStateId) const
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        read_xml(ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);

        for (ptree::value_type& TestPropertyTreeElement : MainConfigPropertyTree.get_child("Settings.Tests"))
            if (ExecuteCellStateId == TestPropertyTreeElement.second.get<UnsignedInt>("<xmlattr>.id"))
            {
                ptree& TestPropertyTreeElementToWriteInFile = TestPropertyTreeElement.second;

                TestPropertyTreeElementToWriteInFile.put("CellStateFileName", CellEngineConfigDataObject.CellStateFileName);

                TestPropertyTreeElementToWriteInFile.put("ChosenStructureIndex", CellEngineConfigDataObject.ChosenStructureIndex);
            }

        std::ofstream OutputConfigFile(ConfigFileName);
        write_xml(OutputConfigFile, MainConfigPropertyTree, boost::property_tree::xml_writer_make_settings<string>('\t', 1));
    }
    CATCH("chess print configuration constructor")
}
