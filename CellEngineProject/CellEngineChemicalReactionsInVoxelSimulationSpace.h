
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SIMULATION_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SIMULATION_SPACE_H

#include "CellEngineParticle.h"

#include "CellEngineChemicalReactionsInSimulationSpace.h"

class CellEngineChemicalReactionsInVoxelSimulationSpace : virtual public CellEngineChemicalReactionsInSimulationSpace, virtual public CellEngineParticlesVoxelsOperations
{
protected:
    void ClearSpaceForParticle(Particle& ParticleObject, bool ClearVoxels) override;
    void MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, SignedInt AddX, SignedInt AddY, SignedInt AddZ) override;
    void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::set<UnsignedInt> &FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
protected:
    explicit CellEngineChemicalReactionsInVoxelSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
