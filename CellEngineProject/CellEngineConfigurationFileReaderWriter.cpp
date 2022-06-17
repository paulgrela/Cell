
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <sb7color.h>

#include "StringUtils.h"
#include "ExceptionsMacro.h"

#include "CellEngineDataFile.h"

#include "CellEnginePDBDataFile.h"
#include "CellEngineCIFDataFile.h"

#include "CellEngineConfigurationFileReaderWriter.h"

using namespace std;

void CellEngineConfigurationFileReaderWriter::ReadChessConfigurationFile(const char* ConfigFileNameParameter, unique_ptr<CellEngineDataFile>& CellEngineDataFileObjectPointer, const uint64_t ExecuteCellStateId)
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
            if (MainConfigPropertyTreeElement.first == "Algorithm")
            {
                MultiThreaded = MainConfigPropertyTreeElement.second.get<bool>("MultiThreaded");
                SetProcessPriorityHighest = MainConfigPropertyTreeElement.second.get<bool>("SetProcessPriorityHighest");
            }
            else
            if (MainConfigPropertyTreeElement.first == "Logger")
            {
                PrintLogToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintToConsole");
                PrintLogToFiles = MainConfigPropertyTreeElement.second.get<bool>("PrintToFiles");

                PrintLogLineNumberToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintLineNumberToConsole");
                PrintLogDateTimeToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintDateTimeToConsole");
                PrintLogProcessIdToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessIdToConsole");
                PrintLogProcessPriorityLevelToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessPriorityLevelToConsole");
                PrintLogThreadIdToConsole = MainConfigPropertyTreeElement.second.get<bool>("PrintThreadIdToConsole");

                PrintLogLineNumberToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintLineNumberToFile");
                PrintLogDateTimeToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintDateTimeToFile");
                PrintLogProcessIdToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessIdToFile");
                PrintLogProcessPriorityLevelToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintProcessPriorityLevelToFile");
                PrintLogThreadIdToFile = MainConfigPropertyTreeElement.second.get<bool>("PrintThreadIdToFile");
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

                        CellEngineDataFileObjectPointer->CellStateFileName = CellStateFileName;

                        CellEngineDataFileObjectPointer->ChosenStructureIndex = CellStatePropertyTreeElement.second.get<uint64_t>("ChosenStructureIndex");
                        CellEngineDataFileObjectPointer->FilmOfStructuresActive = CellStatePropertyTreeElement.second.get<bool>("FilmOfStructuresActive");

                        CellEngineDataFileObjectPointer->SpecularPower = CellStatePropertyTreeElement.second.get<float>("SpecularPower");
                        CellEngineDataFileObjectPointer->SpecularAlbedo = CellStatePropertyTreeElement.second.get<float>("SpecularAlbedo");
                        CellEngineDataFileObjectPointer->MakeColorsTypeObject = static_cast<CellEngineDataFile::MakeColorsType>(CellStatePropertyTreeElement.second.get<uint64_t>("MakeColorsType"));

                        CellEngineDataFileObjectPointer->StencilForDrawingObjectsTypesObject = static_cast<CellEngineDataFile::StencilForDrawingObjectsTypes>(CellStatePropertyTreeElement.second.get<uint64_t>("StencilForDrawingObjectsTypes"));
                        CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops = CellStatePropertyTreeElement.second.get<uint64_t>("NumberOfStencilBufferLoops");

                        CellEngineDataFileObjectPointer->DrawBondsBetweenParticlesCenters = CellStatePropertyTreeElement.second.get<bool>("DrawBondsBetweenParticlesCenters");
                        CellEngineDataFileObjectPointer->DrawBondsBetweenAtoms = CellStatePropertyTreeElement.second.get<bool>("DrawBondsBetweenAtoms");

                        CellEngineDataFileObjectPointer->ShowDetailsInAtomScale = CellStatePropertyTreeElement.second.get<bool>("ShowDetailsInAtomScale");
                        CellEngineDataFileObjectPointer->CheckAtomVisibility = CellStatePropertyTreeElement.second.get<bool>("CheckAtomVisibility");
                        CellEngineDataFileObjectPointer->CutZ = CellStatePropertyTreeElement.second.get<float>("CutZ");
                        CellEngineDataFileObjectPointer->Distance = CellStatePropertyTreeElement.second.get<float>("Distance");
                        CellEngineDataFileObjectPointer->LoadOfAtomsStep = CellStatePropertyTreeElement.second.get<uint64_t>("LoadOfAtomsStep");

                        CellEngineDataFileObjectPointer->XLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("XLowToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->XHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("XHighToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->YLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("YLowToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->YHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("YHighToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->ZLowToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("ZLowToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->ZHighToDrawInAtomScale = CellStatePropertyTreeElement.second.get<float>("ZHighToDrawInAtomScale");

                        CellEngineDataFileObjectPointer->SizeOfAtomsDrawingTypesObject = static_cast<CellEngineDataFile::SizeOfAtomsDrawingTypes>(CellStatePropertyTreeElement.second.get<uint64_t>("SizeOfAtomsDrawingTypes"));

                        CellEngineDataFileObjectPointer->SizeStep = CellStatePropertyTreeElement.second.get<float>("SizeStep");
                        CellEngineDataFileObjectPointer->SizeOfAtomX = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomX");
                        CellEngineDataFileObjectPointer->SizeOfAtomY = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomY");
                        CellEngineDataFileObjectPointer->SizeOfAtomZ = CellStatePropertyTreeElement.second.get<float>("SizeOfAtomZ");
                        CellEngineDataFileObjectPointer->CameraXPosition = CellStatePropertyTreeElement.second.get<float>("CameraXPosition");
                        CellEngineDataFileObjectPointer->CameraYPosition = CellStatePropertyTreeElement.second.get<float>("CameraYPosition");
                        CellEngineDataFileObjectPointer->CameraZPosition = CellStatePropertyTreeElement.second.get<float>("CameraZPosition");
                        CellEngineDataFileObjectPointer->CameraXMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraXMoveShortStep");
                        CellEngineDataFileObjectPointer->CameraYMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraYMoveShortStep");
                        CellEngineDataFileObjectPointer->CameraZMoveShortStep = CellStatePropertyTreeElement.second.get<float>("CameraZMoveShortStep");
                        CellEngineDataFileObjectPointer->CameraXMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraXMoveLongStep");
                        CellEngineDataFileObjectPointer->CameraYMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraYMoveLongStep");
                        CellEngineDataFileObjectPointer->CameraZMoveLongStep = CellStatePropertyTreeElement.second.get<float>("CameraZMoveLongStep");
                        CellEngineDataFileObjectPointer->ViewXMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewXMoveShortStep");
                        CellEngineDataFileObjectPointer->ViewYMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewYMoveShortStep");
                        CellEngineDataFileObjectPointer->ViewZMoveShortStep = CellStatePropertyTreeElement.second.get<float>("ViewZMoveShortStep");
                        CellEngineDataFileObjectPointer->ViewXMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewXMoveLongStep");
                        CellEngineDataFileObjectPointer->ViewYMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewYMoveLongStep");
                        CellEngineDataFileObjectPointer->ViewZMoveLongStep = CellStatePropertyTreeElement.second.get<float>("ViewZMoveLongStep");

                        CellEngineDataFileObjectPointer->ViewChangeUsingLongStep = CellStatePropertyTreeElement.second.get<bool>("ViewChangeUsingLongStep");
                        CellEngineDataFileObjectPointer->AutomaticChangeOfSizeOfAtom = CellStatePropertyTreeElement.second.get<bool>("AutomaticChangeOfSizeOfAtom");
                        CellEngineDataFileObjectPointer->AutomaticChangeOfLoadAtomsStep = CellStatePropertyTreeElement.second.get<bool>("AutomaticChangeOfLoadAtomsStep");

                        CellEngineDataFileObjectPointer->BackgroundColors[1] = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor1")));
                        CellEngineDataFileObjectPointer->BackgroundColors[2] = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor2")));
                        CellEngineDataFileObjectPointer->BackgroundColors[3] = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(CellStatePropertyTreeElement.second.get<string>("BackgroundColor3")));
                        CellEngineDataFileObjectPointer->ChosenBackgroundColor = CellStatePropertyTreeElement.second.get<uint64_t>("ChosenBackgroundColor");

                        CellEngineDataFileObjectPointer->AtomsKinds.clear();
                        for (const ptree::value_type& AtomKindPropertyTreeElement : CellStatePropertyTreeElement.second.get_child("Atoms"))
                        {
                            AtomKind AtomKindObject;

                            LoggersManagerObject.Log(STREAM("Atom Kind Name = " << AtomKindPropertyTreeElement.second.get<string>("Name")));

                            AtomKindObject.Color = sb7::FromVec4ToVec3(sb7::GetColorVec4FromColorName(AtomKindPropertyTreeElement.second.get<string>("Color")));
                            AtomKindObject.SizeX = AtomKindPropertyTreeElement.second.get<float>("SizeX");
                            AtomKindObject.SizeY = AtomKindPropertyTreeElement.second.get<float>("SizeY");
                            AtomKindObject.SizeZ = AtomKindPropertyTreeElement.second.get<float>("SizeZ");
                            CellEngineDataFileObjectPointer->AtomsKinds[AtomKindPropertyTreeElement.second.get<string>("Name")] = AtomKindObject;
                        }

                        CellEngineDataFileObjectPointer->ParticlesKinds.clear();
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
                            CellEngineDataFileObjectPointer->ParticlesKinds[ParticleKindPropertyTreeElement.second.get<UnsignedIntType>("<xmlattr>.id")] = ParticleKindObject;
                        }
                    }
        }
    }
    CATCH("cell print configuration constructor")
}

void CellEngineConfigurationFileReaderWriter::SaveTestStatisticsToFile(unique_ptr<CellEngineDataFile>& CellEngineDataFileObjectPointer, const uint64_t ExecuteCellStateId) const
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

                TestPropertyTreeElementToWriteInFile.put("CellStateFileName", CellEngineDataFileObjectPointer->CellStateFileName);

                TestPropertyTreeElementToWriteInFile.put("ChosenStructureIndex", CellEngineDataFileObjectPointer->ChosenStructureIndex);
            }

        std::ofstream OutputConfigFile(ConfigFileName);
        write_xml(OutputConfigFile, MainConfigPropertyTree, boost::property_tree::xml_writer_make_settings<string>('\t', 1));
    }
    CATCH("chess print configuration constructor")
}
