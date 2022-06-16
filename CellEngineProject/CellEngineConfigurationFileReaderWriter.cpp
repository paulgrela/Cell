
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "ExceptionsMacro.h"

#include "CellEngineDataFile.h"
#include "CellEngineConfigurationFileReaderWriter.h"

using namespace std;

void CellEngineConfigurationFileReaderWriter::ReadChessConfigurationFile(const char* ConfigFileNameParameter, const unique_ptr<CellEngineDataFile>& CellEngineDataFileObjectPointer)
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        ConfigFileName = ConfigFileNameParameter;

        read_xml(ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);
        LoggersManagerObject.Log(STREAM("Reading xml config file finished"));

        uint64_t ExecuteCellStateId = MainConfigPropertyTree.get_child("Settings").get<uint64_t>("ExecuteCellStateId");
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
                    CellEngineState CellEngineStateObject;

                    CellEngineStateObject.CellStateFileName = TestPropertyTreeElement.second.get<string>("CellStateFileName");
                    CellEngineDataFileObjectPointer->ChosenStructureIndex = TestPropertyTreeElement.second.get<uint64_t>("ChosenStructureIndex");

                }
        }

    }
    CATCH("cell print configuration constructor")
}

void CellEngineConfigurationFileReaderWriter::SaveTestStatisticsToFile(const CellEngineState& CellEngineStateObject) const
{
    try
    {
        using boost::property_tree::ptree;

        ptree MainConfigPropertyTree;

        read_xml(this->ConfigFileName, MainConfigPropertyTree, boost::property_tree::xml_parser::trim_whitespace);

        for (ptree::value_type& TestPropertyTreeElement : MainConfigPropertyTree.get_child("Settings.Tests"))
            if (CellEngineStateObject.CellStateId == TestPropertyTreeElement.second.get<uint64_t>("<xmlattr>.id"))
            {
                ptree& TestPropertyTreeElementToWriteInFile = TestPropertyTreeElement.second;
                TestPropertyTreeElementToWriteInFile.put("ChosenStructureIndex", CellEngineStateObject.ChosenStructureIndex);
            }

        std::ofstream OutputConfigFile(this->ConfigFileName);
        write_xml(OutputConfigFile, MainConfigPropertyTree, boost::property_tree::xml_writer_make_settings<string>('\t', 1));
    }
    CATCH("chess print configuration constructor")
}
