
#ifndef CELL_ENGINE_COLORS_H
#define CELL_ENGINE_COLORS_H

#include <chrono>
#include <random>

#include "vmath.h"
#include "CellEngineTypes.h"
#include "CellEngineConfigData.h"

class CellEngineColors
{
private:
    std::mt19937_64 mt64{ std::random_device{}() };
    std::default_random_engine DefaultRandomEngineObject{ static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()) };
    std::uniform_real_distribution<float> UniformDistributionObject;
private:
    vector3_16 DNANR1Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR2Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR3Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR4Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR1Color1 = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR2Color2 = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR3Color3 = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 DNANR4Color4 = GetVector3FormVMathVec3(GetRandomColor());

    vector3_16 RNANR1Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 RNANR2Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 RNANR3Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3_16 RNANR4Color = GetVector3FormVMathVec3(GetRandomColor());
public:
    void SelectRandomEngineForColors();
    vmath::vec3 GetRandomColor();
    vector3_16 GetDNAorRNAColor(EntityIdInt EntityId, ChainIdInt ChainId);
};

inline CellEngineColors CellEngineColorsObject;

#endif
