
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
        return stoi(std::string(ChainName).substr(2, 2));
    }

    inline UnsignedInt GetParticleKindIndexFromChainId(const ChainIdInt ChainId)
    {
        return ChainId > 10 ? (ChainId - 10) : ChainId;;
    }

    inline vector3_16 GetVector3FormVMathVec3(const vmath::vec3& Color)
    {
        return { static_cast<uint16_t>(Color.X() * 100.00f), static_cast<uint16_t>(Color.Y() * 100.00f), static_cast<uint16_t>(Color.Z() * 100.00f) };
    }

    inline vmath::vec3 GetVMathVec3FromVector3(vector3_16 Color)
    {
        return { float(Color.X) / 100.00f, float(Color.Y) / 100.00f, float(Color.Z) / 100.00f };
    }
}

#endif
