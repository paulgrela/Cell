
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <sb7color.h>

#include "StringUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"

#include "CellEnginePDBDataFile.h"
#include "CellEngineCIFDataFile.h"

#include "CellEngineConfigurationFileReaderWriter.h"

using namespace std;

void CellEngineConfigurationFileReaderWriter::ReadChessConfigurationFile(const char* ConfigFileNameParameter, const uint64_t ExecuteCellStateId)
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
            if (MainConfigPropertyTreeElement.first == "WindowParameters")
            {
                CellEngineConfigDataObject.XTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMain");
                CellEngineConfigDataObject.YTopMainWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMain");
                CellEngineConfigDataObject.WidthMainWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMain");
                CellEngineConfigDataObject.HeightMainWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMain");
                CellEngineConfigDataObject.XTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("XTopMenu");
                CellEngineConfigDataObject.YTopMenuWindow = MainConfigPropertyTreeElement.second.get<int>("YTopMenu");
                CellEngineConfigDataObject.WidthMenuWindow = MainConfigPropertyTreeElement.second.get<int>("WidthMenu");
                CellEngineConfigDataObject.HeightMenuWindow = MainConfigPropertyTreeElement.second.get<int>("HeightMenu");
                CellEngineConfigDataObject.XTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecond");
                CellEngineConfigDataObject.YTopSecondWindow = MainConfigPropertyTreeElement.second.get<int>("YTopSecond");
                CellEngineConfigDataObject.WidthSecondWindow = MainConfigPropertyTreeElement.second.get<int>("XTopSecond");
                CellEngineConfigDataObject.HeightSecondWindow = MainConfigPropertyTreeElement.second.get<int>("HeightSecond");
            }
            else
            if (MainConfigPropertyTreeElement.first == "Algorithm")
            {
                CellEngineConfigDataObject.MultiThreaded = MainConfigPropertyTreeElement.second.get<bool>("MultiThreaded");
                CellEngineConfigDataObject.SetProcessPriorityHighest = MainConfigPropertyTreeElement.second.get<bool>("SetProcessPriorityHighest");
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

                CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile = MainConfigPropertyTreeElement.second.get<uint64_t>("MaximalNumberOfLinesInOneFile");
            }
            else
            if (MainConfigPropertyTreeElement.first == "CellsStates")
                for (const ptree::value_type& CellStatePropertyTreeElement : MainConfigPropertyTree.get_child("Settings.CellsStates"))
                    if (CellStatePropertyTreeElement.second.get<uint64_t>("<xmlattr>.id") == ExecuteCellStateId)
                    {
                        auto CellStateFileName = CellStatePropertyTreeElement.second.get<string>("CellStateFileName");

                        if (string_utils::check_end_str(CellStateFileName, ".pdb") == true)
                            CellEngineDataFileObjectPointer = make_unique<CellEnginePDBDataFile>();
                        else
                            CellEngineDataFileObjectPointer = make_unique<CellEngineCIFDataFile>();

                        CellEngineConfigDataObject.CellStateFileName = CellStateFileName;

                        CellEngineConfigDataObject.ChosenStructureIndex = CellStatePropertyTreeElement.second.get<uint64_t>("ChosenStructureIndex");
                        CellEngineConfigDataObject.FilmOfStructuresActive = CellStatePropertyTreeElement.second.get<bool>("FilmOfStructuresActive");

                        CellEngineConfigDataObject.SpecularPower = CellStatePropertyTreeElement.second.get<float>("SpecularPower");
                        CellEngineConfigDataObject.SpecularAlbedo = CellStatePropertyTreeElement.second.get<float>("SpecularAlbedo");
                        CellEngineConfigDataObject.MakeColorsTypeObject = static_cast<CellEngineConfigData::MakeColorsType>(CellStatePropertyTreeElement.second.get<uint64_t>("MakeColorsType"));

                        CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject = static_cast<CellEngineConfigData::StencilForDrawingObjectsTypes>(CellStatePropertyTreeElement.second.get<uint64_t>("StencilForDrawingObjectsTypes"));
                        CellEngineConfigDataObject.NumberOfStencilBufferLoops = CellStatePropertyTreeElement.second.get<uint64_t>("NumberOfStencilBufferLoops");

                        CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters = CellStatePropertyTreeElement.second.get<bool>("DrawBondsBetweenParticlesCenters");
                        CellEngineConfigDataObject.DrawBondsBetweenAtoms = CellStatePropertyTreeElement.second.get<bool>("DrawBondsBetweenAtoms");

                        CellEngineConfigDataObject.ShowDetailsInAtomScale = CellStatePropertyTreeElement.second.get<bool>("ShowDetailsInAtomScale");
                        CellEngineConfigDataObject.CheckAtomVisibility = CellStatePropertyTreeElement.second.get<bool>("CheckAtomVisibility");
                        CellEngineConfigDataObject.CutZ = CellStatePropertyTreeElement.second.get<float>("CutZ");
                        CellEngineConfigDataObject.Distance = CellStatePropertyTreeElement.second.get<float>("Distance");
                        CellEngineConfigDataObject.LoadOfAtomsStep = CellStatePropertyTreeElement.second.get<uint64_t>("LoadOfAtomsStep");

                        CellEngineConfigDataObject.XLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("XLowToDrawInAtomScale");
                        CellEngineConfigDataObject.XHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("XHighToDrawInAtomScale");
                        CellEngineConfigDataObject.YLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("YLowToDrawInAtomScale");
                        CellEngineConfigDataObject.YHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("YHighToDrawInAtomScale");
                        CellEngineConfigDataObject.ZLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("ZLowToDrawInAtomScale");
                        CellEngineConfigDataObject.ZHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("ZHighToDrawInAtomScale");

                        CellEngineConfigDataObject.SizeOfAtomsDrawingTypesObject = static_cast<CellEngineConfigData::SizeOfAtomsDrawingTypes>(CellStatePropertyTreeElement.second.get<uint64_t>("SizeOfAtomsDrawingTypes"));

                        CellEngineConfigDataObject.SizeStep = CellStatePropertyTreeElement.second.get<float>("SizeStep");
                        CellEngineConfigDataObject.SizeOfAtomX = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomX");
                        CellEngineConfigDataObject.SizeOfAtomY = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomY");
                        CellEngineConfigDataObject.SizeOfAtomZ = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomZ");
                        CellEngineConfigDataObject.RotationAngle1 = CellStatePropertyTreeElement.second.get<float>("RotationAngle1");
                        CellEngineConfigDataObject.RotationAngle2 = CellStatePropertyTreeElement.second.get<float>("RotationAngle2");
                        CellEngineConfigDataObject.RotationAngle3 = CellStatePropertyTreeElement.second.get<float>("RotationAngle3");
                        CellEngineConfigDataObject.ViewX = CellStatePropertyTreeElement.second.get<float>("ViewX");
                        CellEngineConfigDataObject.ViewY = CellStatePropertyTreeElement.second.get<float>("ViewY");
                        CellEngineConfigDataObject.ViewZ = CellStatePropertyTreeElement.second.get<float>("ViewZ");
                        CellEngineConfigDataObject.CameraXPosition = CellStatePropertyTreeElement.second.get<float>("CameraXPosition");
                        CellEngineConfigDataObject.CameraYPosition = CellStatePropertyTreeElement.second.get<float>("CameraYPosition");
                        CellEngineConfigDataObject.CameraZPosition = CellStatePropertyTreeElement.second.get<float>("CameraZPosition");
                        CellEngineConfigDataObject.CameraXMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraXMoveShortStep");
                        CellEngineConfigDataObject.CameraYMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraYMoveShortStep");
                        CellEngineConfigDataObject.CameraZMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraZMoveShortStep");
                        CellEngineConfigDataObject.CameraXMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraXMoveLongStep");
                        CellEngineConfigDataObject.CameraYMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraYMoveLongStep");
                        CellEngineConfigDataObject.CameraZMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraZMoveLongStep");
                        CellEngineConfigDataObject.ViewX = CellStatePropertyTreeElement.second.get<float>("ViewX");
                        CellEngineConfigDataObject.ViewY = CellStatePropertyTreeElement.second.get<float>("ViewY");
                        CellEngineConfigDataObject.ViewZ = CellStatePropertyTreeElement.second.get<float>("ViewZ");
                        CellEngineConfigDataObject.ViewXMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewXMoveShortStep");
                        CellEngineConfigDataObject.ViewYMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewYMoveShortStep");
                        CellEngineConfigDataObject.ViewZMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewZMoveShortStep");
                        CellEngineConfigDataObject.ViewXMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewXMoveLongStep");
                        CellEngineConfigDataObject.ViewYMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewYMoveLongStep");
                        CellEngineConfigDataObject.ViewZMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewZMoveLongStep");

                        CellEngineConfigDataObject.ViewChangeUsingLongStep = CellStatePropertyTreeElement.second.get<bool>("ViewChangeUsingLongStep");
                        CellEngineConfigDataObject.AutomaticChangeOfSizeOfAtom = CellStatePropertyTreeElement.second.get<bool>("AutomaticChangeOfSizeOfAtom");
                        CellEngineConfigDataObject.AutomaticChangeOfLoadAtomsStep = CellStatePropertyTreeElement.second.get<bool>("AutomaticChangeOfLoadAtomsStep");

                        CellEngineConfigDataObject.BackgroundColors[1] = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor1")));
                        CellEngineConfigDataObject.BackgroundColors[2] = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor2")));
                        CellEngineConfigDataObject.BackgroundColors[3] = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor3")));
                        CellEngineConfigDataObject.ChosenBackgroundColor = CellStatePropertyTreeElement.second.get<uint64_t>("ChosenBackgroundColor");

                        CellEngineConfigDataObject.AtomsKinds.clear();
                        for (const ptree::value_type& AtomKindPropertyTreeElement : CellStatePropertyTreeElement.second.get_child("Atoms"))
                        {
                            AtomKind AtomKindObject;

                            LoggersManagerObject.Log(STREAM("Atom Kind Name = " << AtomKindPropertyTreeElement.second.get<string>("Name")));

                            AtomKindObject.Color = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(AtomKindPropertyTreeElement.second.get<string>("Color")));
                            AtomKindObject.SizeX = AtomKindPropertyTreeElement.second.get<float>("SizeX");
                            AtomKindObject.SizeY = AtomKindPropertyTreeElement.second.get<float>("SizeY");
                            AtomKindObject.SizeZ = AtomKindPropertyTreeElement.second.get<float>("SizeZ");
                            AtomKindObject.Name = AtomKindPropertyTreeElement.second.get<string>("Name");
                            //CellEngineConfigDataObject.AtomsKinds[AtomKindPropertyTreeElement.second.get<string>("Name")] = AtomKindObject;
                            CellEngineConfigDataObject.AtomsKinds.emplace_back(AtomKindObject);
                        }

                        CellEngineConfigDataObject.ParticlesKindsXML.clear();
                        for (const ptree::value_type& ParticleKindPropertyTreeElement : CellStatePropertyTreeElement.second.get_child("Particles"))
                        {
                            ParticleKind ParticleKindObject;

                            LoggersManagerObject.Log(STREAM("Particle Kind Name = " << ParticleKindPropertyTreeElement.second.get<string>("Name") << " ID = " << ParticleKindPropertyTreeElement.second.get<UnsignedIntType>("<xmlattr>.id") << " COLOR = " << ParticleKindPropertyTreeElement.second.get<string>("Color")));

                            ParticleKindObject.NameFromXML = ParticleKindPropertyTreeElement.second.get<string>("Name");
                            ParticleKindObject.Color = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(ParticleKindPropertyTreeElement.second.get<string>("Color")));
                            ParticleKindObject.Visible = ParticleKindPropertyTreeElement.second.get<bool>("Visible");
                            ParticleKindObject.SizeX = ParticleKindPropertyTreeElement.second.get<float>("SizeX");
                            ParticleKindObject.SizeY = ParticleKindPropertyTreeElement.second.get<float>("SizeY");
                            ParticleKindObject.SizeZ = ParticleKindPropertyTreeElement.second.get<float>("SizeZ");
                            CellEngineConfigDataObject.ParticlesKindsXML[ParticleKindPropertyTreeElement.second.get<UnsignedIntType>("<xmlattr>.id")] = ParticleKindObject;
                        }
                    }
        }
    }
    CATCH("cell print configuration constructor")
}

void CellEngineConfigurationFileReaderWriter::SaveTestStatisticsToFile(const uint64_t ExecuteCellStateId) const
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        read_xml(ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);

        for (ptree::value_type& TestPropertyTreeElement : MainConfigPropertyTree.get_child("Settings.Tests"))
            if (ExecuteCellStateId == TestPropertyTreeElement.second.get<uint64_t>("<xmlattr>.id"))
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
