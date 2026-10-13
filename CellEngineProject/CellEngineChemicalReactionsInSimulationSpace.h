
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_SIMULATION_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_SIMULATION_SPACE_H

#include "CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.h"

class CellEngineChemicalReactionsInSimulationSpace : public CellEngineNucleicAcidsChemicalReactionsInSimulationSpace
{
protected:
    UnsignedInt GenerateCombinationsStateNumber{};
protected:
    virtual void ClearSpaceForParticle(Particle& ParticleObject, bool ClearVoxels) = 0;
public:
    void RemoveParticle(UniqueIdUInt64Type ParticleUniqueId, bool ClearVoxels, bool CountAllParticlesBool) override;
protected:
    void MakingZeroSizeForContainersForFoundParticlesInProximity(const ThreadPosType& CurrentThreadPos);
    void UpdateFoundNucleotidesForFoundParticlesInProximity(UnsignedInt ParticleUniqueId);
    void PrintInformationAboutFoundParticlesInProximity(const ThreadPosType& CurrentThreadPos);
protected:
    bool FindParticlesInProximityOfSimulationSpaceForSelectedSpace(bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);

    virtual void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(MainSetType<UnsignedInt>& FoundParticleUniqueIdes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;

    void SaveParticleFoundInProximity(UniqueIdUInt64Type ParticleUniqueId, MainSetType<UnsignedInt>& FoundParticleUniqueIdes, bool UpdateNucleotides);
protected:
    explicit CellEngineChemicalReactionsInSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineNucleicAcidsChemicalReactionsInSimulationSpace(ParticlesParam)
    {
    }
};

#endif
