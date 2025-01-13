
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_IN_VOXEL_SPACE_H

#include "CellEngineBasicParticlesOperations.h"
#include "CellEngineConfigData.h"
#include "CellEngineConfigurationFileReaderWriter.h"

struct ThreadLocalParticlesInProximity
{
public:
    std::map<EntityIdInt, UnsignedInt> ParticlesKindsFoundInProximity;
    std::vector<UniqueIdInt> ParticlesSortedByCapacityFoundInProximity;
public:
    std::vector<UniqueIdInt> NucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> NucleotidesWithFreePrevEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreeNextEndingsFoundInProximity;
    std::vector<UniqueIdInt> DNANucleotidesWithFreePrevEndingsFoundInProximity;
public:
    std::vector<UniqueIdInt> NucleotidesFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFoundInProximity;
public:
    std::vector<UniqueIdInt> DNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdInt> RNANucleotidesFullFreeFoundInProximity;
    std::vector<UniqueIdInt> tRNAChargedFoundInProximity;
    std::vector<UniqueIdInt> tRNAUnchargedFoundInProximity;
};

class CellEngineChemicalReactionsInBasicSimulationSpace : virtual public CellEngineBasicParticlesOperations
{
protected:
    ThreadLocalParticlesInProximity LocalThreadParticlesInProximityObject;
protected:
    static bool CompareFitnessOfParticle(const ParticleKindForChemicalReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    void EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, ListOfElements &Centers);
protected:
    explicit CellEngineChemicalReactionsInBasicSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
    ~CellEngineChemicalReactionsInBasicSimulationSpace() override = default;
};

#endif
