
#include "CellEngineExecutionTimeStatistics.h"
#include "CellEngineChemicalReactionsInFullAtomSimulationSpace.h"

void CellEngineChemicalReactionsInFullAtomSimulationSpace::FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, const bool UpdateNucleotides, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        for (const auto& ParticleObject : Particles[StartXPosParam][StartYPosParam][StartZPosParam].Particles)
        {
            LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.emplace_back(ParticleObject.second.Index);
            LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ParticleObject.second.EntityId]++;
        }

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForFindingParticles += chrono::duration(stop_time - start_time);
    }
    CATCH("finding particles in proximity of simulation space for selected local space")
};

void CellEngineChemicalReactionsInFullAtomSimulationSpace::MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const RealType AddX, const RealType AddY, const RealType AddZ)
{
    try
    {
        MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(ParticleObject, Particles, CurrentSectorPos, NewPositionParticleObject, AddX, AddY, AddZ, CurrentThreadPos);
    }
    CATCH("moving particle near other particle if space is empty or to near space")
}

void CellEngineChemicalReactionsInFullAtomSimulationSpace::ClearSpaceForParticle(Particle& ParticleObject, const bool ClearFullAtoms)
{
}