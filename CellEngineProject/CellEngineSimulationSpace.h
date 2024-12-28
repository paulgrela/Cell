#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include <stack>
#include <barrier>
#include <unordered_set>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineConfigData.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineCompiledDataCreator.h"
#include "CellEngineIllinoisDataCreator.h"
#include "CellEngineBasicVoxelsOperations.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineBasicParticlesOperations.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"
#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"
#include "CellEngineTestParticlesInVoxelSpaceGenerator.h"
#include "CellEngineChemicalReactionsInBasicSimulationSpace.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"
#include "CellEngineVoxelSimulationSpaceStatistics.h"
#include "CellEngineSimulationSpaceStatistics.h"

#define SIMULATION_DETAILED_LOG

class CellEngineSimulationSpace : public CellEngineChemicalReactionsEngine, public CellEngineIllinoisDataCreator, public CellEngineCompiledDataCreator, virtual public CellEngineChemicalReactionsInSimulationSpace, public CellEngineSimulationSpaceStatistics
{
private:
    UnsignedInt ErrorCounter = 0;
    UnsignedInt AddedParticlesInReactions = 0;
private:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
protected:
    virtual bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;
    virtual bool MoveParticleByVectorIfSpaceIsEmpty(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ) = 0;
    virtual bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, UnsignedInt VectorX, UnsignedInt VectorY, UnsignedInt VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) = 0;
protected:
    virtual void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_64 NewPointElement) = 0;
public:
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateOneStepOfElectricDiffusionForOneParticle(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt ParticleIndexParam, ElectricChargeType (*NeighbourPoints)[3][3][3], UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedSpace(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], double MultiplyElectricChargeFactor);
public:
    void SetMakeSimulationStepNumberZero();
    void SetIncSimulationStepNumber();
    void SaveParticlesStatisticsOnce();
    void SaveReactionsStatisticsToFile() const;
public:
    void SaveNumberOfParticlesStatisticsToFile();
public:
    void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateNStepsOfDiffusionForBigPartOfCellSpace(bool InBounds, UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfDiffusionForWholeCellSpace(bool InBounds, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
public:
    void GenerateNStepsOfOneChosenReactionForWholeCellSpace(UnsignedInt ReactionId, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfOneRandomReactionForWholeCellSpace(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateOneStepOfRandomReactionsForAllParticles(const CurrentThreadPosType& CurrentThreadPos);
public:
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
    void GenerateOneStepOfRandomReactionsForOneParticleFromRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt ShiftIndexOfChosenParticle);
    void GenerateOneRandomReactionForChosenParticle(const Particle& ParticleObject);
    void GenerateOneRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, bool FindParticlesInProximityBool);
    void GenerateOneChosenReactionForSelectedSpace(UnsignedInt ReactionId, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex, bool StateOfVoxelSpaceDivisionForThreads);
    void GenerateNStepsOfSimulationForWholeCellSpaceInThreads(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
public:
    void GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside, bool PrintTime);
public:
    void ExchangeParticlesBetweenThreads(UnsignedInt StepOutside, bool StateOfVoxelSpaceDivisionForThreads, bool PrintInfo) const;
    void ExchangeParticlesBetweenThreadsParallelInsert(UnsignedInt StepOutside, bool StateOfVoxelSpaceDivisionForThreads, bool PrintInfo) const;
    void ExchangeParticlesBetweenThreadsParallelExtract(UnsignedInt StepOutside, bool StateOfVoxelSpaceDivisionForThreads, bool PrintInfo) const;
    void GatherParticlesFromThreadsToParticlesInMainThread();
    void FirstSendParticlesForThreads(bool PrintCenterOfParticleWithThreadIndex, bool PrintTime);
    void JoinStatisticsFromThreads() const;
public:
    void CheckParticlesCenters(bool PrintAllParticles) const;
protected:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject);
    SimulationSpaceSectorBounds GetBoundsForThreadSector() const;
protected:
    void PrepareRandomReaction();
    void FindAndExecuteRandomReactionVersion1(UnsignedInt MaxNumberOfReactantsParam);
    void FindAndExecuteRandomReactionVersion2(UnsignedInt MaxNumberOfReactantsParam);
    void FindAndExecuteRandomReactionVersion3(UnsignedInt MaxNumberOfReactantsParam);
    void FindAndExecuteRandomReactionVersion4(UnsignedInt MaxNumberOfReactantsParam);
protected:
    std::vector<UnsignedInt> GetRandomParticlesVersion1(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants);
    std::vector<UnsignedInt> GetRandomParticlesVersion2(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants);
    std::vector<UnsignedInt> GetRandomParticlesVersion3(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants);
public:
    void FindAndExecuteRandomReaction(UnsignedInt MaxNumberOfReactantsParam);
    bool FindAndExecuteChosenReaction(UnsignedInt ReactionId);
protected:
    bool CancelChemicalReaction(const vector<UniqueIdInt>& CreatedParticlesIndexes, chrono::high_resolution_clock::time_point start_time, const ParticleKind& ParticleKindObjectForProduct, char PlaceStr);
    void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, const ParticleKind& ParticleKindObjectForProduct, UnsignedInt VectorX, UnsignedInt VectorY, UnsignedInt VectorZ);
    bool PlaceProductParticleInSpaceInDeterminedPositionOrCancelReaction(UniqueIdInt ParticleIndex, const vector<UniqueIdInt>& CreatedParticlesIndexes, UnsignedInt CenterIndex, const vector<vector3_16>& Centers, const ParticleKind& ParticleKindObjectForProduct, chrono::high_resolution_clock::time_point start_time);
    bool PlaceProductParticleInSpaceInRandomPositionOrCancelReaction(UniqueIdInt ParticleIndex, const vector<UniqueIdInt>& CreatedParticlesIndexes, UnsignedInt CenterIndex, const vector<vector3_16>& Centers, const ParticleKind& ParticleKindObjectForProduct, chrono::high_resolution_clock::time_point start_time);
protected:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants) override;
    bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject) override;
    bool MakeChemicalReaction(ChemicalReaction& ReactionObject) override;
public:
    explicit CellEngineSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineSimulationSpaceStatistics(), Particles(ParticlesParam)
    {
        SetMakeSimulationStepNumberZero();
    }
    ~CellEngineSimulationSpace() override = default;
};

#endif
