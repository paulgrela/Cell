
#ifndef CELL_ENGINE_USEFUL_H
#define CELL_ENGINE_USEFUL_H

#include "CellEngineTypes.h"
#include "CellEngineConfigData.h"

namespace CellEngineUseful
{
    inline bool IsDNAorRNA(const EntityIdInt EntityId)
    {
        return EntityId == CellEngineConfigDataObject.DNAIdentifier || EntityId == CellEngineConfigDataObject.RNAIdentifier;
    }

    inline bool IsNucleotide(const std::string_view ChainName)
    {
        return (ChainName == "NU01" || ChainName == "NU11" || ChainName == "NU02" || ChainName == "NU12" || ChainName == "NU03" || ChainName == "NU13" || ChainName == "NU04" || ChainName == "NU14");
    }

    inline ChainIdInt GetChainIdFromChainName(const std::string_view ChainName)
    {
        return ChainName.substr(0, 2) == "NU" || ChainName.substr(0, 2) == "NV" ? stoi(std::string(ChainName).substr(2, 2)) : 0;
    }

    inline UnsignedInt GetParticleKindIndexFromChainId(const ChainIdInt ChainId)
    {
        return ChainId > 10 ? (ChainId - 10) : ChainId;
    }

    inline char GetLetterForDNAChainId0(const ChainIdInt ChainId)
    {
        switch (ChainId)
        {
            case 1 : return 'A';
            case 2 : return 'C';
            case 3 : return 'G';
            case 4 : return 'T';

            case 11 : return 'T';
            case 12 : return 'G';
            case 13 : return 'C';
            case 14 : return 'A';

            default : break;
        }

        return 'I';
    }

    inline char GetLetterForDNAChainId1(const ChainIdInt ChainId)
    {
        switch (ChainId)
        {
            case 1 : return 'T';
            case 2 : return 'G';
            case 3 : return 'C';
            case 4 : return 'A';

            case 11 : return 'A';
            case 12 : return 'C';
            case 13 : return 'G';
            case 14 : return 'T';

            default : break;
        }

        return 'I';
    }

    inline vector3_16 GetVector3FormVMathVec3ForColor(const vmath::vec3& Color)
    {
        return { static_cast<uint16_t>(Color.X() * 100.00f), static_cast<uint16_t>(Color.Y() * 100.00f), static_cast<uint16_t>(Color.Z() * 100.00f) };
    }

    inline vmath::vec3 GetVMathVec3FromVector3ForColor(vector3_16 Color)
    {
        return { float(Color.X) / 100.00f, float(Color.Y) / 100.00f, float(Color.Z) / 100.00f };
    }
}

#endif
