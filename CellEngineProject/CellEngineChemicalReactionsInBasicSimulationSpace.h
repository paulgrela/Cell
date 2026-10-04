
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_BASIC_SIMULATION_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_BASIC_SIMULATION_SPACE_H

#include "CellEngineBasicParticlesOperations.h"
#include "CellEngineConfigData.h"
#include "CellEngineConfigurationFileReaderWriter.h"

struct ThreadLocalParticlesInProximity
{
public:
    MainMapType<EntityIdInt, UnsignedInt> ParticlesKindsFoundInProximity;
    std::vector<UniqueIdUnsignedInt> ParticlesSortedByCapacityFoundInProximity;
public:
    std::vector<UniqueIdUnsignedInt> NucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdUnsignedInt> NucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdUnsignedInt> DNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdUnsignedInt> DNANucleotidesWithFreePrevEndingsFoundInProximity;
public:
    std::vector<UniqueIdUnsignedInt> NucleotidesFreeFoundInProximity;
    std::vector<UniqueIdUnsignedInt> RNANucleotidesFreeFoundInProximity;
    std::vector<UniqueIdUnsignedInt> RNANucleotidesFoundInProximity;
public:
    std::vector<UniqueIdUnsignedInt> DNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdUnsignedInt> RNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdUnsignedInt> tRNAChargedFoundInProximity;
    std::vector<UniqueIdUnsignedInt> tRNAUnchargedFoundInProximity;
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
