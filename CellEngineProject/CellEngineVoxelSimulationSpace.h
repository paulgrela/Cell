#pragma once

#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_H

#include <stack>
#include <unordered_set>

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineConfigData.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineIllinoisDataCreator.h"
#include "CellEngineBasicVoxelsOperations.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineBasicParticlesOperations.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"
#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"
#include "CellEngineTestParticlesInVoxelSpaceGenerator.h"
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"
#include "CellEngineVoxelSimulationSpaceStatistics.h"
#include "CellEngineChemicalReactionsInVoxelSimulationSpace.h"
#include "CellEngineSimulationSpace.h"

class CellEngineVoxelSimulationSpace : public CellEngineSimulationSpace, public CellEngineTestParticlesInVoxelSpaceGenerator, public CellEngineChemicalReactionsInVoxelSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
{
protected:
    UnsignedInt XMinGlobal{}, XMaxGlobal{}, YMinGlobal{}, YMaxGlobal{}, ZMinGlobal{}, ZMaxGlobal{};
public:
    SimulationSpaceVoxel GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
    Particle& GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex);
private:
    ParticlesContainer<Particle>& Particles;
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    void SetAtomInVoxelSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
    Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) override;
public:
    void ClearVoxelSpaceAndParticles() override;
public:
    void ClearWholeVoxelSpace();
    void ClearSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_64 NewPointElement) override;
    bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ, float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam) override;
    // bool MoveParticleByVectorIfSpaceIsEmpty(Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ) override;
    bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) override;
protected:
    void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, UnsignedInt VectorX, UnsignedInt VectorY, UnsignedInt VectorZ) override;
public:
    template <class T>
    static void CheckParticlesIndexes(ParticlesContainerInternal<T>& FormerParticlesIndexes, const string& FormerState);
    void CheckCancelledParticlesIndexes();
    void CheckFormerExistedParticlesIndexes();
public:
    explicit CellEngineVoxelSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, bool GetMemoryForVoxelSpace, ThreadIdType ThreadIndexParam, const CurrentThreadPosType& CurrentThreadPos);
    ~CellEngineVoxelSimulationSpace() override;
};

#endif
