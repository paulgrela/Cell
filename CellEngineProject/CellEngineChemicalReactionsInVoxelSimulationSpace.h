
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SIMULATION_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SIMULATION_SPACE_H

#include "CellEngineParticle.h"

#include "CellEngineChemicalReactionsInSimulationSpace.h"

class CellEngineChemicalReactionsInVoxelSimulationSpace : virtual public CellEngineChemicalReactionsInSimulationSpace, virtual public CellEngineParticlesVoxelsOperations
{
protected:
    void ClearSpaceForParticle(Particle& ParticleObject, bool ClearVoxels) override;
    void MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, RealType AddX, RealType AddY, RealType AddZ) override;
    void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(MainSetType<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
protected:
    explicit CellEngineChemicalReactionsInVoxelSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
