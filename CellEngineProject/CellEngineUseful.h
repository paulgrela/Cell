
#ifndef CELL_ENGINE_USEFUL_H
#define CELL_ENGINE_USEFUL_H

#include "CellEngineTypes.h"
#include "CellEngineConstants.h"
#include "CellEngineConfigData.h"

#define SIMULATION_DETAILED_LOG

namespace CellEngineUseful
{
    template<typename T>
    bool IsIn(T first, std::vector<T> values)
    {
        return find(values.begin(), values.end(), first) != values.end();
    }

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
            case ACode : return 'A';
            case CCode : return 'C';
            case GCode : return 'G';
            case TCode : return 'T';
            case UCode : return 'U';

            case TCodeAdd : return 'T';
            case GCodeAdd : return 'G';
            case CCodeAdd : return 'C';
            case ACodeAdd : return 'A';
            case UCodeAdd : return 'U';

            default : break;
        }

        return 'I';
    }

    inline ChainIdInt GetChainIdFromLetterForDNAorRNA(const char Letter)
    {
        switch (Letter)
        {
            case 'A' : return ACode;
            case 'C' : return CCode;
            case 'G' : return GCode;
            case 'T' : return TCode;
            case 'U' : return UCode;

            case 'R' : return RCode;
            case 'Y' : return YCode;
            case 'S' : return SCode;
            case 'W' : return WCode;
            case 'K' : return KCode;
            case 'M' : return MCode;
            case 'B' : return BCode;
            case 'D' : return DCode;
            case 'H' : return HCode;
            case 'V' : return VCode;

            case 'N' : return NCode;

            default : break;
        }

        return 0;
    }

    inline ChainIdInt GetPairedChainIdForDNAorRNA(const ChainIdInt ChainId)
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

    bool inline CompareIUPACNucleotideCode(const ChainIdInt LetterIUPAC, const ChainIdInt LetterToCheck)
    {
        if (LetterIUPAC == NCode)
            return true;

        if (LetterIUPAC == ACode && LetterToCheck == ACode)
            return true;

        if (LetterIUPAC == CCode && LetterToCheck == CCode)
            return true;

        if (LetterIUPAC == GCode && LetterToCheck == GCode)
            return true;

        if (LetterIUPAC == TCode && LetterToCheck == TCode)
            return true;

        if (LetterIUPAC == RCode && (LetterToCheck == ACode || LetterToCheck == GCode))
            return true;

        if (LetterIUPAC == YCode && (LetterToCheck == CCode || LetterToCheck == TCode))
            return true;

        if (LetterIUPAC == SCode && (LetterToCheck == GCode || LetterToCheck == CCode))
            return true;

        if (LetterIUPAC == WCode && (LetterToCheck == ACode || LetterToCheck == TCode))
            return true;

        if (LetterIUPAC == KCode && (LetterToCheck == GCode || LetterToCheck == TCode))
            return true;

        if (LetterIUPAC == MCode && (LetterToCheck == ACode || LetterToCheck == CCode))
            return true;

        if (LetterIUPAC == BCode && (LetterToCheck == CCode || LetterToCheck == GCode || LetterToCheck == TCode))
            return true;

        if (LetterIUPAC == DCode && (LetterToCheck == ACode || LetterToCheck == GCode || LetterToCheck == TCode))
            return true;

        if (LetterIUPAC == HCode && (LetterToCheck == ACode || LetterToCheck == CCode || LetterToCheck == TCode))
            return true;

        if (LetterIUPAC == VCode && (LetterToCheck == ACode || LetterToCheck == CCode || LetterToCheck == GCode))
            return true;

        return false;
    }

    bool inline CompareIUPACNucleotideLetterCode(const char LetterIUPAC, const char LetterToCheck)
    {
        if (LetterIUPAC == 'N')
            return true;

        if (LetterIUPAC == 'A' && LetterToCheck == 'A')
            return true;

        if (LetterIUPAC == 'C' && LetterToCheck == 'C')
            return true;

        if (LetterIUPAC == 'G' && LetterToCheck == 'G')
            return true;

        if (LetterIUPAC == 'T' && LetterToCheck == 'T')
            return true;

        if (LetterIUPAC == 'R' && (LetterToCheck == 'A' || LetterToCheck == 'G'))
            return true;

        if (LetterIUPAC == 'Y' && (LetterToCheck == 'C' || LetterToCheck == 'T'))
            return true;

        if (LetterIUPAC == 'S' && (LetterToCheck ==	'G' || LetterToCheck == 'C'))
            return true;

        if (LetterIUPAC == 'W' && (LetterToCheck ==	'A' || LetterToCheck == 'T'))
            return true;

        if (LetterIUPAC == 'K' && (LetterToCheck ==	'G' || LetterToCheck == 'T'))
            return true;

        if (LetterIUPAC == 'M' && (LetterToCheck ==	'A' || LetterToCheck == 'C'))
            return true;

        if (LetterIUPAC == 'B' && (LetterToCheck ==	'C' || LetterToCheck == 'G' || LetterToCheck == 'T'))
            return true;

        if (LetterIUPAC == 'D' && (LetterToCheck ==	'A' || LetterToCheck == 'G' || LetterToCheck == 'T'))
            return true;

        if (LetterIUPAC == 'H' && (LetterToCheck ==	'A' || LetterToCheck == 'C' || LetterToCheck == 'T'))
            return true;

        if (LetterIUPAC == 'V' && (LetterToCheck ==	'A' || LetterToCheck == 'C' || LetterToCheck == 'G'))
            return true;

        return false;
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
