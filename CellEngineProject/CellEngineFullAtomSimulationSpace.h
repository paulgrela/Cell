
#ifndef CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_H
#define CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_H

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineIllinoisDataCreator.h"
#include "CellEngineBasicVoxelsOperations.h"
#include "CellEngineTestParticlesInVoxelSpaceGenerator.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"
#include "CellEngineVoxelSimulationSpaceStatistics.h"
#include "CellEngineChemicalReactionsInFullAtomSimulationSpace.h"
#include "CellEngineRealRandomParticlesInFullAtomSpaceGenerator.h"
#include "CellEngineSimulationSpace.h"

//class CellEngineFullAtomSimulationSpace : public CellEngineSimulationSpace, public CellEngineTestParticlesInVoxelSpaceGenerator, public CellEngineChemicalReactionsInFullAtomSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
class CellEngineFullAtomSimulationSpace : public CellEngineSimulationSpace, public CellEngineRealRandomParticlesInFullAtomSpaceGenerator, public CellEngineChemicalReactionsInFullAtomSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
{
protected:
    float XMinGlobal{}, XMaxGlobal{}, YMinGlobal{}, YMaxGlobal{}, ZMinGlobal{}, ZMaxGlobal{};
public:
    Particle& GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex);
private:
    ParticlesContainer<Particle>& Particles;
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) override;
public:
    void ClearVoxelSpaceAndParticles() override;
public:
    void ClearWholeVoxelSpace();
    void ClearSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_64 NewPointElement) override;
    bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ, float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam) override;
    bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) override;
protected:
    void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, UnsignedInt VectorX, UnsignedInt VectorY, UnsignedInt VectorZ) override;
public:
    explicit CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, bool GetMemoryForVoxelSpace, ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPos);
    ~CellEngineFullAtomSimulationSpace() override;
};

#endif
