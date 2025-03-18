
#ifndef CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_H
#define CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_H

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineIllinoisDataCreator.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"
#include "CellEngineVoxelSimulationSpaceStatistics.h"
#include "CellEngineChemicalReactionsInFullAtomSimulationSpace.h"
#include "CellEngineRealRandomParticlesInFullAtomSpaceGenerator.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineGenomeNucleicAcidsParticlesInFullAtomSpaceGenerator.h"

class CellEngineFullAtomSimulationSpace : public CellEngineSimulationSpace, public CellEngineGenomeNucleicAcidsParticlesInFullAtomSpaceGenerator, public CellEngineChemicalReactionsInFullAtomSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
{
protected:
    RealType XMinGlobal{}, XMaxGlobal{}, YMinGlobal{}, YMaxGlobal{}, ZMinGlobal{}, ZMaxGlobal{};
public:
    Particle& GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex);
private:
    ParticlesContainer<Particle>& Particles;
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) override;
public:
    void ClearFullAtomSpaceAndParticles() override;
public:
    static void ClearWholeFullAtomSpace();
    static void ClearSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, RealType VectorX, RealType VectorY, RealType VectorZ) override;
    void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_Real32 NewPointElement) override;
    bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, RealType VectorX, RealType VectorY, RealType VectorZ, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) override;
    bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, RealType VectorX, RealType VectorY, RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) override;
protected:
    bool CheckInsertOfParticle(const MPIParticleSenderStruct& MPIParticleSenderToInsert) override;
public:
    void WriteNumberOfParticlesInEachSectorToFile() const;
    static void WriteNumberOfParticlesKindsWithoutAtoms();
protected:
    SimulationSpaceSectorBounds GetBoundsForThreadSector() override;
public:
    void GenerateNStepsOfDiffusionForWholeCellSpace(bool InBounds, RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps) override;
    void GenerateNStepsOfOneRandomReactionForWholeCellSpace(RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps) override;
    void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) override;
    void GenerateOneRandomReactionForSelectedSpace(RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam, bool FindParticlesInProximityBool) override;
public:
    explicit CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, bool GetMemoryForFullAtomSpace, ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPos);
    ~CellEngineFullAtomSimulationSpace() override;
};

#endif
