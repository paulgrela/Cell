
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
#include "CellEngineChemicalReactionsInVoxelSimulationSpace.h"
#include "CellEngineSimulationSpace.h"

class CellEngineFullAtomSimulationSpace : public CellEngineSimulationSpace, public CellEngineTestParticlesInVoxelSpaceGenerator, public CellEngineChemicalReactionsInVoxelSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
{
protected:
    float XMinGlobal{}, XMaxGlobal{}, YMinGlobal{}, YMaxGlobal{}, ZMinGlobal{}, ZMaxGlobal{};
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
    // void SetAtomInFullAtomSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
    Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) override;
public:
    void ClearVoxelSpaceAndParticles() override;
public:
    void ClearWholeVoxelSpace();
    void ClearSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_64 NewPointElement) override;
    bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
    bool MoveParticleByVectorIfSpaceIsEmpty(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ) override;
    bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, UnsignedInt VectorX, UnsignedInt VectorY, UnsignedInt VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) override;
protected:
    void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, UnsignedInt VectorX, UnsignedInt VectorY, UnsignedInt VectorZ) override;
public:
    explicit CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, bool GetMemoryForVoxelSpace, ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPos);
    ~CellEngineFullAtomSimulationSpace() override;
};

#endif
