
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
    void EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleUniqueIdChosenForReaction, ListOfCentersType &Centers, std::vector<Particle>& ParticlesBackup);
protected:
    explicit CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
        #ifdef CONTAINERS_FOR_SPEED
        LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.reserve(1024);
        #endif
        LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.reserve(10240);

        LocalThreadParticlesInProximityObject.NucleotidesWithFreeNextEndingsFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.NucleotidesWithFreePrevEndingsFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.DNANucleotidesWithFreeNextEndingsFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.DNANucleotidesWithFreePrevEndingsFoundInProximity.reserve(1024);

        LocalThreadParticlesInProximityObject.NucleotidesFreeFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.RNANucleotidesFreeFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.reserve(1024);

        LocalThreadParticlesInProximityObject.DNANucleotidesFullFreeFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.RNANucleotidesFullFreeFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.tRNAChargedFoundInProximity.reserve(1024);
        LocalThreadParticlesInProximityObject.tRNAUnchargedFoundInProximity.reserve(1024);
    }
    ~CellEngineChemicalReactionsInBasicSimulationSpace() override = default;
};

#endif
