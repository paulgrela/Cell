
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_FULL_ATOM_SIMULATION_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_FULL_ATOM_SIMULATION_SPACE_H

#include "CellEngineParticle.h"

#include "CellEngineParticlesFullAtomOperations.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

class CellEngineChemicalReactionsInFullAtomSimulationSpace : virtual public CellEngineChemicalReactionsInSimulationSpace, virtual public CellEngineParticlesFullAtomOperations
{
protected:
    void ClearSpaceForParticle(Particle& ParticleObject, bool ClearVoxels) override;
    void MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, float AddX, float AddY, float AddZ) override;
    void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
protected:
    explicit CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
