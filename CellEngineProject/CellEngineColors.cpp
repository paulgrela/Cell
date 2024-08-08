#define USING_MODULES_

#ifdef USING_MODULES
module;
#endif

#include <chrono>
#include <random>

#include "vmath.h"
#include "CellEngineTypes.h"
#include "CellEngineUseful.h"

#include "CellEngineConfigData.h"

#ifdef USING_MODULES
export module CellEngineColors;

export class CellEngineColors;

class CellEngineColors
{
private:
    std::mt19937_64 mt64{ std::random_device{}() };
    std::default_random_engine DefaultRandomEngineObject{ static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()) };
    std::uniform_real_distribution<float> UniformDistributionObject;
private:
    vector3_16 DNANR1Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR2Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR3Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR4Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR1Color1 = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR2Color2 = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR3Color3 = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 DNANR4Color4 = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());

    vector3_16 RNANR1Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 RNANR2Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 RNANR3Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
    vector3_16 RNANR4Color = CellEngineUseful::GetVector3FormVMathVec3ForColor(GetRandomColor());
public:
    void SelectRandomEngineForColors();
    vmath::vec3 GetRandomColor();
    vector3_16 GetDNAorRNAColor(EntityIdInt EntityId, ChainIdInt ChainId);
};
#else
#include "CellEngineColors.h"
#endif

void CellEngineColors::SelectRandomEngineForColors()
{
    try
    {
        switch (CellEngineConfigDataObject.RandomColorEngineObject)
        {
            case CellEngineConfigData::RandomColorEngineTypes::Rand : srand((unsigned int)time(nullptr)); break;
            case CellEngineConfigData::RandomColorEngineTypes::mt19937 : mt64.seed(std::random_device{}()); break;
            case CellEngineConfigData::RandomColorEngineTypes::DefaultRandomEngine: DefaultRandomEngineObject.seed(static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count())); break;
            default:  break;
        }
    }
    CATCH("selecting random engines for colors");
}

vmath::vec3 CellEngineColors::GetRandomColor()
{
    vmath::vec3 ReturnRandomColor;

    try
    {
        switch (CellEngineConfigDataObject.RandomColorEngineObject)
        {
            case CellEngineConfigData::RandomColorEngineTypes::Rand : ReturnRandomColor = vmath::vec3((float)rand() / static_cast<float>(RAND_MAX), (float)rand() / static_cast<float>(RAND_MAX), (float)rand() / static_cast<float>(RAND_MAX)); break;
            case CellEngineConfigData::RandomColorEngineTypes::mt19937 : ReturnRandomColor = vmath::vec3(UniformDistributionObject(mt64), UniformDistributionObject(mt64), UniformDistributionObject(mt64)); break;
            case CellEngineConfigData::RandomColorEngineTypes::DefaultRandomEngine : ReturnRandomColor = vmath::vec3(UniformDistributionObject(DefaultRandomEngineObject), UniformDistributionObject(DefaultRandomEngineObject), UniformDistributionObject(DefaultRandomEngineObject));break;
            default:  break;
        }
    }
    CATCH("getting random color")

    return ReturnRandomColor;
}

vector3_16 CellEngineColors::GetDNAorRNAColor(EntityIdInt EntityId, ChainIdInt ChainId)
{
    if (EntityId == CellEngineConfigDataObject.DNAIdentifier)
    {
        switch (ChainId)
        {
            case 01: return DNANR1Color1;
            case 11: return DNANR1Color;
            case 02: return DNANR2Color2;
            case 12: return DNANR2Color;
            case 03: return DNANR3Color3;
            case 13: return DNANR3Color;
            case 04: return DNANR4Color4;
            case 14: return DNANR4Color;
            default : break;
        }
    }
    else
    if (EntityId == CellEngineConfigDataObject.RNAIdentifier)
    {
        switch (ChainId)
        {
            case 01: return RNANR1Color;
            case 02: return RNANR2Color;
            case 03: return RNANR3Color;
            case 04: return RNANR4Color;
            default : break;
        }
    }

    return { 0, 0, 0 };
}

#ifdef USING_MODULES
export CellEngineColors CellEngineColorsObject;
#endif