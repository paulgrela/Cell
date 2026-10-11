
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_BASIC_SIMULATION_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_BASIC_SIMULATION_SPACE_H

#include "CellEngineBasicParticlesOperations.h"

struct ThreadLocalParticlesInProximity
{
public:
    MainMapType<UniqueIdUInt32Type, UnsignedInt> ParticlesKindsFoundInProximity;
    std::vector<UniqueIdUInt64Type> ParticlesSortedByCapacityFoundInProximity;
public:
    std::vector<UniqueIdUInt64Type> NucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdUInt64Type> NucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdUInt64Type> DNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdUInt64Type> DNANucleotidesWithFreePrevEndingsFoundInProximity;
public:
    std::vector<UniqueIdUInt64Type> NucleotidesFreeFoundInProximity;
    std::vector<UniqueIdUInt64Type> RNANucleotidesFreeFoundInProximity;
    std::vector<UniqueIdUInt64Type> RNANucleotidesFoundInProximity;
public:
    std::vector<UniqueIdUInt64Type> DNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdUInt64Type> RNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdUInt64Type> tRNAChargedFoundInProximity;
    std::vector<UniqueIdUInt64Type> tRNAUnchargedFoundInProximity;
};

class CellEngineChemicalReactionsInBasicSimulationSpace : virtual public CellEngineBasicParticlesOperations
{
protected:
    ThreadLocalParticlesInProximity LocalThreadParticlesInProximityObject;
protected:
    static bool CompareFitnessOfParticle(const ParticleKindForChemicalReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    void EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, ListOfCentersType &Centers, std::vector<Particle>& ParticlesBackup);
protected:
    explicit CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
        #ifdef CONTAINERS_FOR_SPEED
        LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.reserve(1000);
        #endif
        LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.reserve(10000);
    }
    ~CellEngineChemicalReactionsInBasicSimulationSpace() override = default;
};

#endif
