
#ifndef CELL_ENGINE_USEFUL_H
#define CELL_ENGINE_USEFUL_H

#include "CellEngineTypes.h"
#include "CellEngineConfigData.h"

#define SIMULATION_DETAILED_LOG

namespace CellEngineUseful
{
    static void SwitchOffLogs()
    {
        #ifdef SIMULATION_DETAILED_LOG
        LoggersManagerObject.InitializePrintingParameters(false, false, false, false, false, false, false, false, false, false, false, false, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile, false);
        #endif
    }

    static void SwitchOnLogs()
    {
        #ifdef SIMULATION_DETAILED_LOG
        LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile, CellEngineConfigDataObject.PrintLogToCommonFileWhenPrintLogToSpecialFile);
        #endif
    }

    inline bool IsDNAorRNA(const EntityIdInt EntityId)
    {
        return EntityId == CellEngineConfigDataObject.DNAIdentifier || EntityId == CellEngineConfigDataObject.RNAIdentifier;
    }

    inline bool IsDNA(const EntityIdInt EntityId)
    {
        return EntityId == CellEngineConfigDataObject.DNAIdentifier;
    }

    inline bool IsRNA(const EntityIdInt EntityId)
    {
        return EntityId == CellEngineConfigDataObject.RNAIdentifier;
    }

    inline bool IsNucleotide(const std::string_view ChainName)
    {
        return (ChainName == "NU01" || ChainName == "NU11" || ChainName == "NU02" || ChainName == "NU12" || ChainName == "NU03" || ChainName == "NU13" || ChainName == "NU04" || ChainName == "NU14");
    }

    inline ChainIdInt GetChainIdFromChainName(const std::string_view ChainName)
    {
        return ChainName.substr(0, 2) == "NU" || ChainName.substr(0, 2) == "NV" ? stoi(std::string(ChainName).substr(2, 2)) : 0;
    }

    inline char GetLetterFromChainIdForDNAorRNA(const ChainIdInt ChainId)
    {
        switch (ChainId)
        {
            case 1 : return 'A';
            case 2 : return 'C';
            case 3 : return 'G';
            case 4 : return 'T';
            case 5 : return 'U';

            case 11 : return 'T';
            case 12 : return 'G';
            case 13 : return 'C';
            case 14 : return 'A';
            case 15 : return 'U';

            default : break;
        }

        return 'I';
    }

    inline ChainIdInt GetChainIdFromLetterForDNAorRNA(const char Letter)
    {
        switch (Letter)
        {
            case 'A' : return 1;
            case 'C' : return 2;
            case 'G' : return 3;
            case 'T' : return 4;
            case 'U' : return 5;

            default : break;
        }

        return 0;
    }

    inline ChainIdInt GetPairedChainIdForDNAorRNA(ChainIdInt ChainId)
    {
        switch (ChainId)
        {
            case 1 : return 4;
            case 2 : return 3;
            case 3 : return 2;
            case 4 : return 1;
            case 5 : return 1;
            default : break;
        }

        return 0;
    }

    inline vector3_16 GetVector3FormVMathVec3ForColor(const vmath::vec3& Color)
    {
        return { static_cast<uint16_t>(Color.X() * 100.00f), static_cast<uint16_t>(Color.Y() * 100.00f), static_cast<uint16_t>(Color.Z() * 100.00f) };
    }

    inline vmath::vec3 GetVMathVec3FromVector3ForColor(vector3_16 Color)
    {
        return { float(Color.X) / 100.00f, float(Color.Y) / 100.00f, float(Color.Z) / 100.00f };
    }

    class AtomDescriptionTexts
    {
    public:
        std::string Texts[6];
    };

    inline AtomDescriptionTexts AtomDescriptionTextsObject;
}

#endif
