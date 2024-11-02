
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
    std::vector<std::vector<std::vector<ThreadLocalParticlesInProximity>>> ThreadsLocalParticlesInProximity;
protected:
    ThreadLocalParticlesInProximity& GetThreadsLocalParticlesInProximity(const CurrentThreadPosType& CurrentThreadPos)
    {
        return ThreadsLocalParticlesInProximity[CurrentThreadPos.ThreadPosX][CurrentThreadPos.ThreadPosY][CurrentThreadPos.ThreadPozZ];
    }
protected:
    void ConstructDataForMultiThreadedExecution()
    {
        ThreadsLocalParticlesInProximity.clear();
        ThreadsLocalParticlesInProximity.resize(CellEngineConfigDataObject.NumberOfXThreads);
        for (auto& ThreadLocalParticlesInProximityXPos : ThreadsLocalParticlesInProximity)
        {
            ThreadLocalParticlesInProximityXPos.resize(CellEngineConfigDataObject.NumberOfYThreads);
            for (auto& ThreadLocalParticlesInProximityYPos : ThreadLocalParticlesInProximityXPos)
                ThreadLocalParticlesInProximityYPos.resize(CellEngineConfigDataObject.NumberOfZThreads);
        }
    }
protected:
    static bool CompareFitnessOfParticle(const ParticleKindForChemicalReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction);
    void EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(UnsignedInt ParticleIndexChosenForReaction, std::vector <vector3_16> &Centers);
protected:
    explicit CellEngineChemicalReactionsInBasicSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
    ~CellEngineChemicalReactionsInBasicSimulationSpace() override = default;
};

#endif
