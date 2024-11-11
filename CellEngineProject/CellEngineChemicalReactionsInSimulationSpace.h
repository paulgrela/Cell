
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
    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleIndex, const CurrentThreadPosType& CurrentThreadPos);
    void PrintInformationAboutFoundParticlesInProximity(const CurrentThreadPosType& CurrentThreadPos);
protected:
    bool FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, UnsignedInt AdditionalBoundFactor, const CurrentThreadPosType& CurrentThreadPos);

    bool FindParticlesInProximityOfSimulationSpaceForSelectedSpace(bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);

    virtual void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos) = 0;

    void SaveParticleFoundInProximity(const UniqueIdInt ParticleIndex, unordered_set<UnsignedInt>& FoundParticleIndexes, const bool UpdateNucleotides, const CurrentThreadPosType& CurrentThreadPos);
protected:
    explicit CellEngineChemicalReactionsInSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(ParticlesParam)
    {
    }
};

#endif
