
#ifndef CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H
#define CELL_ENGINE_NUCLEIC_ACIDS_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H

#include <unordered_set>

#include "CellEngineExecutionTimeStatistics.h"

#include "CellEngineChemicalReaction.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineNucleicAcidsComplexOperations.h"
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"
#include "CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.h"

class CellEngineChemicalReactionsInSimulationSpace : public CellEngineNucleicAcidsChemicalReactionsInSimulationSpace
{
protected:
    UnsignedInt GenerateCombinationsStateNumber{};
protected:
    virtual void ClearSpaceForParticle(Particle& ParticleObject, bool ClearVoxels) = 0;
protected:
    void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) override;
protected:
    void MakingZeroSizeForContainersForFoundParticlesInProximity(const CurrentThreadPosType& CurrentThreadPos);
    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleIndex);
    void PrintInformationAboutFoundParticlesInProximity(const CurrentThreadPosType& CurrentThreadPos);
protected:
    bool FindParticlesInProximityOfSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor);

    bool FindParticlesInProximityOfSimulationSpaceForSelectedSpace(bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);

    virtual void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;

    void SaveParticleFoundInProximity(UniqueIdInt ParticleIndex, unordered_set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides);
protected:
    explicit CellEngineChemicalReactionsInSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(ParticlesParam)
    {
    }
};

#endif
