
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"

using namespace std;

bool CellEngineChemicalReactionsInBasicSimulationSpace::CompareFitnessOfParticle(const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction)
{
    return
        (ParticleKindForReactionObject.LinkedParticleTypes.empty() == true ||
         (ParticleKindForReactionObject.LinkedParticleTypes.empty() == false && ParticleObjectForReaction.LinkedParticlesPointersList.size() == ParticleKindForReactionObject.LinkedParticleTypes.size()
         && all_of(ParticleObjectForReaction.LinkedParticlesPointersList.begin(), ParticleObjectForReaction.LinkedParticlesPointersList.end(), [](const Particle* PointerToParticle){ return PointerToParticle != nullptr; })
         && equal(ParticleObjectForReaction.LinkedParticlesPointersList.begin(), ParticleObjectForReaction.LinkedParticlesPointersList.end(), ParticleKindForReactionObject.LinkedParticleTypes.begin(), [](const Particle* PointerToParticle, UniqueIdInt ParticleType){ return PointerToParticle->EntityId == ParticleType; })));
}

void CellEngineChemicalReactionsInBasicSimulationSpace::EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(const UnsignedInt ParticleIndexChosenForReaction, vector<vector3_16>& Centers)
{
    try
    {
        auto& ParticleObjectToBeErased = GetParticleFromIndex(ParticleIndexChosenForReaction);
        Centers.emplace_back(ParticleObjectToBeErased.Center.X, ParticleObjectToBeErased.Center.Y, ParticleObjectToBeErased.Center.Z);
        LoggersManagerObject.Log(STREAM("Centers - X = " << to_string(ParticleObjectToBeErased.Center.X) << " Y = " << to_string(ParticleObjectToBeErased.Center.Y) << " Z = " << to_string(ParticleObjectToBeErased.Center.Z) << endl));

        RemoveParticle(ParticleIndexChosenForReaction, true);
    }
    CATCH("erasing particles chosen for reaction and get centers for new products of reaction")
}