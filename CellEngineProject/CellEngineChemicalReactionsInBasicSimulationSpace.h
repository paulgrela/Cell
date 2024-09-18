
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H

#include "CellEngineBasicParticlesOperations.h"

class CellEngineChemicalReactionsInBasicSimulationSpace : virtual public CellEngineBasicParticlesOperations
{
protected:
    std::map<EntityIdInt, UnsignedInt> ParticlesKindsFoundInProximity;
    std::vector<UniqueIdInt> ParticlesSortedByCapacityFoundInProximity;
protected:
    std::vector<UniqueIdInt> NucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> NucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreePrevEndingsFoundInProximity;
protected:
    std::vector<UniqueIdInt> NucleotidesFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFoundInProximity;
protected:
    std::vector<UniqueIdInt> DNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdInt> tRNAChargedFoundInProximity;
    std::vector<UniqueIdInt> tRNAUnchargedFoundInProximity;
protected:
    static bool CompareFitnessOfParticle(const ParticleKindForChemicalReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    void EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, std::vector <vector3_16> &Centers);
protected:
    explicit CellEngineChemicalReactionsInBasicSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
    virtual ~CellEngineChemicalReactionsInBasicSimulationSpace() = default;
};

#endif
