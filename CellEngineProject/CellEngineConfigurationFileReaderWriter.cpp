
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

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
                for (const ptree::value_type& TestPropertyTreeElement : MainConfigPropertyTree.get_child("Settings.CellsStates"))
                    if (TestPropertyTreeElement.second.get<uint64_t>("<xmlattr>.id") == ExecuteCellStateId)
                    {
                        auto CellStateFileName = TestPropertyTreeElement.second.get<string>("CellStateFileName");

                        if (string_utils::check_end_str(CellStateFileName, ".pdb") == true)
                            CellEngineDataFileObjectPointer = make_unique<CellEnginePDBDataFile>();
                        else
                            CellEngineDataFileObjectPointer = make_unique<CellEngineCIFDataFile>();

                        CellEngineDataFileObjectPointer->CellStateFileName = CellStateFileName;

                        CellEngineDataFileObjectPointer->ChosenStructureIndex = TestPropertyTreeElement.second.get<uint64_t>("ChosenStructureIndex");
                        CellEngineDataFileObjectPointer->FilmOfStructuresActive = TestPropertyTreeElement.second.get<bool>("FilmOfStructuresActive");

                        CellEngineDataFileObjectPointer->SpecularPower = TestPropertyTreeElement.second.get<float>("SpecularPower");
                        CellEngineDataFileObjectPointer->SpecularAlbedo = TestPropertyTreeElement.second.get<float>("SpecularAlbedo");
                        CellEngineDataFileObjectPointer->MakeColorsTypeObject = static_cast<CellEngineDataFile::MakeColorsType>(TestPropertyTreeElement.second.get<uint64_t>("MakeColorsType"));

                        CellEngineDataFileObjectPointer->StencilForDrawingObjectsTypesObject = static_cast<CellEngineDataFile::StencilForDrawingObjectsTypes>(TestPropertyTreeElement.second.get<uint64_t>("StencilForDrawingObjectsTypes"));
                        CellEngineDataFileObjectPointer->NumberOfStencilBufferLoops = TestPropertyTreeElement.second.get<uint64_t>("NumberOfStencilBufferLoops");

                        CellEngineDataFileObjectPointer->DrawBondsBetweenParticlesCenters = TestPropertyTreeElement.second.get<bool>("DrawBondsBetweenParticlesCenters");
                        CellEngineDataFileObjectPointer->DrawBondsBetweenAtoms = TestPropertyTreeElement.second.get<bool>("DrawBondsBetweenAtoms");

                        CellEngineDataFileObjectPointer->ShowDetailsInAtomScale = TestPropertyTreeElement.second.get<bool>("ShowDetailsInAtomScale");
                        CellEngineDataFileObjectPointer->CheckAtomVisibility = TestPropertyTreeElement.second.get<bool>("CheckAtomVisibility");
                        CellEngineDataFileObjectPointer->CutZ = TestPropertyTreeElement.second.get<float>("CutZ");
                        CellEngineDataFileObjectPointer->Distance = TestPropertyTreeElement.second.get<float>("Distance");
                        CellEngineDataFileObjectPointer->LoadOfAtomsStep = TestPropertyTreeElement.second.get<uint64_t>("LoadOfAtomsStep");

                        CellEngineDataFileObjectPointer->XLowToDrawInAtomScale = TestPropertyTreeElement.second.get<float>("XLowToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->XHighToDrawInAtomScale = TestPropertyTreeElement.second.get<float>("XHighToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->YLowToDrawInAtomScale = TestPropertyTreeElement.second.get<float>("YLowToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->YHighToDrawInAtomScale = TestPropertyTreeElement.second.get<float>("YHighToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->ZLowToDrawInAtomScale = TestPropertyTreeElement.second.get<float>("ZLowToDrawInAtomScale");
                        CellEngineDataFileObjectPointer->ZHighToDrawInAtomScale = TestPropertyTreeElement.second.get<float>("ZHighToDrawInAtomScale");

                        CellEngineDataFileObjectPointer->SizeOfAtomsDrawingTypesObject = static_cast<CellEngineDataFile::SizeOfAtomsDrawingTypes>(TestPropertyTreeElement.second.get<uint64_t>("SizeOfAtomsDrawingTypes"));

                        CellEngineDataFileObjectPointer->SizeStep = TestPropertyTreeElement.second.get<float>("SizeStep");
                        CellEngineDataFileObjectPointer->SizeOfAtomX = TestPropertyTreeElement.second.get<float>("SizeOfAtomX");
                        CellEngineDataFileObjectPointer->SizeOfAtomY = TestPropertyTreeElement.second.get<float>("SizeOfAtomY");
                        CellEngineDataFileObjectPointer->SizeOfAtomZ = TestPropertyTreeElement.second.get<float>("SizeOfAtomZ");
                        CellEngineDataFileObjectPointer->CameraXPosition = TestPropertyTreeElement.second.get<float>("CameraXPosition");
                        CellEngineDataFileObjectPointer->CameraYPosition = TestPropertyTreeElement.second.get<float>("CameraYPosition");
                        CellEngineDataFileObjectPointer->CameraZPosition = TestPropertyTreeElement.second.get<float>("CameraZPosition");
                        CellEngineDataFileObjectPointer->CameraXMoveShortStep = TestPropertyTreeElement.second.get<float>("CameraXMoveShortStep");
                        CellEngineDataFileObjectPointer->CameraYMoveShortStep = TestPropertyTreeElement.second.get<float>("CameraYMoveShortStep");
                        CellEngineDataFileObjectPointer->CameraZMoveShortStep = TestPropertyTreeElement.second.get<float>("CameraZMoveShortStep");
                        CellEngineDataFileObjectPointer->CameraXMoveLongStep = TestPropertyTreeElement.second.get<float>("CameraXMoveLongStep");
                        CellEngineDataFileObjectPointer->CameraYMoveLongStep = TestPropertyTreeElement.second.get<float>("CameraYMoveLongStep");
                        CellEngineDataFileObjectPointer->CameraZMoveLongStep = TestPropertyTreeElement.second.get<float>("CameraZMoveLongStep");
                        CellEngineDataFileObjectPointer->ViewXMoveShortStep = TestPropertyTreeElement.second.get<float>("ViewXMoveShortStep");
                        CellEngineDataFileObjectPointer->ViewYMoveShortStep = TestPropertyTreeElement.second.get<float>("ViewYMoveShortStep");
                        CellEngineDataFileObjectPointer->ViewZMoveShortStep = TestPropertyTreeElement.second.get<float>("ViewZMoveShortStep");
                        CellEngineDataFileObjectPointer->ViewXMoveLongStep = TestPropertyTreeElement.second.get<float>("ViewXMoveLongStep");
                        CellEngineDataFileObjectPointer->ViewYMoveLongStep = TestPropertyTreeElement.second.get<float>("ViewYMoveLongStep");
                        CellEngineDataFileObjectPointer->ViewZMoveLongStep = TestPropertyTreeElement.second.get<float>("ViewZMoveLongStep");

                        CellEngineDataFileObjectPointer->ViewChangeUsingLongStep = TestPropertyTreeElement.second.get<bool>("ViewChangeUsingLongStep");
                        CellEngineDataFileObjectPointer->AutomaticChangeOfSizeOfAtom = TestPropertyTreeElement.second.get<bool>("AutomaticChangeOfSizeOfAtom");
                        CellEngineDataFileObjectPointer->AutomaticChangeOfLoadAtomsStep = TestPropertyTreeElement.second.get<bool>("AutomaticChangeOfLoadAtomsStep");

                        CellEngineDataFileObjectPointer->BackgroundColors[1] = FromVec4ToVec3(GetColorVec4FromColorName(TestPropertyTreeElement.second.get<string>("BackgroundColor1")));
                        CellEngineDataFileObjectPointer->BackgroundColors[2] = FromVec4ToVec3(GetColorVec4FromColorName(TestPropertyTreeElement.second.get<string>("BackgroundColor2")));
                        CellEngineDataFileObjectPointer->BackgroundColors[3] = FromVec4ToVec3(GetColorVec4FromColorName(TestPropertyTreeElement.second.get<string>("BackgroundColor3")));
                        CellEngineDataFileObjectPointer->ChosenBackgroundColor = TestPropertyTreeElement.second.get<uint64_t>("ChosenBackgroundColor");
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

        read_xml(this->ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);

        for (ptree::value_type& TestPropertyTreeElement : MainConfigPropertyTree.get_child("Settings.Tests"))
            if (ExecuteCellStateId == TestPropertyTreeElement.second.get<uint64_t>("<xmlattr>.id"))
            {
                ptree& TestPropertyTreeElementToWriteInFile = TestPropertyTreeElement.second;

                TestPropertyTreeElementToWriteInFile.put("CellStateFileName", CellEngineDataFileObjectPointer->CellStateFileName);

                TestPropertyTreeElementToWriteInFile.put("ChosenStructureIndex", CellEngineDataFileObjectPointer->ChosenStructureIndex);
            }

        std::ofstream OutputConfigFile(this->ConfigFileName);
        write_xml(OutputConfigFile, MainConfigPropertyTree, boost::property_tree::xml_writer_make_settings<string>('\t', 1));
    }
    CATCH("chess print configuration constructor")
}
