
#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineCompiledDataCreator.h"
#include "CellEngineIllinoisDataCreator.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"
#include "CellEngineSimulationParallelExecutionManager.h"
#include "CellEngineSimulationSpaceStatistics.h"

using IndexesChosenForReactionType = vector<pair<UniqueIdUInt64Type, UnsignedInt>>;

class CellEngineSimulationSpace : public CellEngineChemicalReactionsEngine, public CellEngineIllinoisDataCreator, public CellEngineCompiledDataCreator, virtual public CellEngineChemicalReactionsInSimulationSpace, public CellEngineSimulationSpaceStatistics, public CellEngineSimulationParallelExecutionManager
{
private:
    ParticlesContainer<Particle>& Particles;
protected:
    virtual void FillParticleElementsInSpace(UniqueIdUInt64Type ParticleUniqueId, const ParticleKind& ParticleKindObjectForProduct, RealType VectorX, RealType VectorY, RealType VectorZ) = 0;
    virtual void FillParticleElementInSpace(UniqueIdUInt64Type ParticleUniqueId, vector3_Real32 NewPointElement) = 0;
protected:
    virtual bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, RealType VectorX, RealType VectorY, RealType VectorZ, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) = 0;

    virtual bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsReactions(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, RealType VectorX, RealType VectorY, RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) = 0;
    virtual bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsDiffusion(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, RealType VectorX, RealType VectorY, RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) = 0;
protected:
    virtual SimulationSpaceSectorBounds GetBoundsForThreadSector() = 0;
public:
    void GenerateOneStepOfElectricDiffusionForOneParticle(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdUInt64Type ParticleUniqueIdParam, ElectricChargeType (*NeighbourPoints)[3][3][3], UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(const Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], double MultiplyElectricChargeFactor);
public:
    void SetMakeSimulationStepNumberZero();
    void SetIncSimulationStepNumber();
    void SaveParticlesStatisticsOnce();
    void SaveReactionsStatisticsToFile() const;
public:
    void SaveHistogramOfParticlesStatisticsToFile() const;
    void SaveNumberOfParticlesStatisticsToFile();
public:
    virtual void GenerateNStepsOfDiffusionForWholeCellSpace(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, bool InBounds, RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps) = 0;
    virtual void GenerateNStepsOfOneRandomReactionForWholeCellSpace(RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps) = 0;
protected:
    std::tuple<IndexesChosenForReactionType, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject);
protected:
    void PrepareRandomReaction();
    void FindAndExecuteRandomReactionVersion3(UnsignedInt MaxNumberOfReactantsParam);
private:
    set<UnsignedInt> GetAllPossibleReactionsFromParticlesInProximity();
public:
    void FindAndExecuteRandomReaction(UnsignedInt MaxNumberOfReactantsParam);
    bool FindAndExecuteChosenReaction(UnsignedInt ReactionId);
protected:
    bool CancelChemicalReaction(const vector<UniqueIdUInt64Type>& CreatedParticlesIndexes, const ListOfCentersType& Centers, const vector<Particle>& ParticlesBackup, chrono::high_resolution_clock::time_point start_time, const ParticleKind& ParticleKindObjectForProduct, char PlaceStr);
    bool PlaceNewProductParticleInSpaceDeterminedFromPositionOfFormerReactantParticleOrCancelReaction(UniqueIdUInt64Type ParticleUniqueId, const vector<Particle>& ParticlesBackup, const vector<UniqueIdUInt64Type>& CreatedParticlesIndexes, UnsignedInt CenterIndex, const ListOfCentersType& Centers, ParticleKind& ParticleKindObjectForProduct, chrono::high_resolution_clock::time_point start_time);
    bool PlaceNewProductParticleInSpaceInNewRandomPositionOrCancelReaction(UniqueIdUInt64Type ParticleUniqueId, const vector<Particle>& ParticlesBackup, const vector<UniqueIdUInt64Type>& CreatedParticlesIndexes, UnsignedInt CenterIndex, const ListOfCentersType& Centers, ParticleKind& ParticleKindObjectForProduct, chrono::high_resolution_clock::time_point start_time);
protected:
    bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject) override;
    bool MakeChemicalReaction(ChemicalReaction& ReactionObject) override;
protected:
    #ifdef WELL_STIRRED
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants) override;
    std::vector<UnsignedInt> GetRandomParticlesVersion3() const;
    #endif
public:
    explicit CellEngineSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineSimulationSpaceStatistics(), Particles(ParticlesParam), CellEngineSimulationParallelExecutionManager()
    {
        SetMakeSimulationStepNumberZero();
    }
    ~CellEngineSimulationSpace() override = default;
};

#endif
