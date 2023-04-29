
#ifndef CELL_ENGINE_COLORS_H
#define CELL_ENGINE_COLORS_H

#include <chrono>
#include <random>

#include "vmath.h"
#include "CellEngineTypes.h"
#include "CellEngineUseful.h"
#include "CellEngineConfigData.h"

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

inline CellEngineColors CellEngineColorsObject;

#endif
