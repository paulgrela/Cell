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
    std::unordered_map<UniqueIdInt, Particle>& Particles;
protected:
    virtual bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) = 0;
    virtual bool MoveParticleByVectorIfSpaceIsEmpty(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ) = 0;
protected:
    virtual void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_64 NewVoxel) = 0;
public:
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateOneStepOfElectricDiffusionForOneParticle(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt ParticleIndexParam, ElectricChargeType (*NeighbourPoints)[3][3][3], UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);
    void GenerateOneStepOfElectricDiffusionForSelectedSpace(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);
    void GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);
protected:
    void UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], double MultiplyElectricChargeFactor, const CurrentThreadPosType& CurrentThreadPos);
public:
    void SetMakeSimulationStepNumberZero();
    void SetIncSimulationStepNumber();
    void SaveParticlesStatisticsOnce();
    void SaveReactionsStatisticsToFile();
public:
    void SaveNumberOfParticlesStatisticsToFile();
public:
    void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);
public:
    void GenerateNStepsOfDiffusionForBigPartOfCellSpace(bool InBounds, UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, const CurrentThreadPosType& CurrentThreadPos, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfDiffusionForWholeCellSpace(bool InBounds, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, const CurrentThreadPosType& CurrentThreadPos, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, const CurrentThreadPosType& CurrentThreadPos, UnsignedInt NumberOfSimulationSteps);
public:
    void GenerateNStepsOfOneChosenReactionForWholeCellSpace(UnsignedInt ReactionId, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, const CurrentThreadPosType& CurrentThreadPos, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfOneRandomReactionForWholeCellSpace(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, const CurrentThreadPosType& CurrentThreadPos, UnsignedInt NumberOfSimulationSteps);
    void GenerateOneStepOfRandomReactionsForAllParticles(const CurrentThreadPosType& CurrentThreadPos);
public:
    void GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(UnsignedInt NumberOfStepsInside, UnsignedInt StepOutside, UnsignedInt ThreadXIndex, UnsignedInt ThreadYIndex, UnsignedInt ThreadZIndex);
    void GenerateNStepsOfSimulationForWholeCellSpaceInThreads(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside);
public:
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const CurrentThreadPosType& CurrentThreadPos);
    void GenerateOneStepOfRandomReactionsForOneParticleFromRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt ShiftIndexOfChosenParticle, const CurrentThreadPosType& CurrentThreadPos);
    void GenerateOneRandomReactionForChosenParticle(const Particle& ParticleObject, const CurrentThreadPosType& CurrentThreadPos);
    void GenerateOneRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);
    void GenerateOneChosenReactionForSelectedSpace(UnsignedInt ReactionId, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos);
protected:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos);
protected:
    void PrepareRandomReaction();
    void FindAndExecuteRandomReactionVersion1(UnsignedInt MaxNumberOfReactantsParam, const CurrentThreadPosType& CurrentThreadPos);
    void FindAndExecuteRandomReactionVersion2(UnsignedInt MaxNumberOfReactantsParam, const CurrentThreadPosType& CurrentThreadPos);
    void FindAndExecuteRandomReactionVersion3(UnsignedInt MaxNumberOfReactantsParam, const CurrentThreadPosType& CurrentThreadPos);
    void FindAndExecuteRandomReactionVersion4(UnsignedInt MaxNumberOfReactantsParam, const CurrentThreadPosType& CurrentThreadPos);
protected:
    std::vector<UnsignedInt> GetRandomParticlesVersion1(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants, const CurrentThreadPosType& CurrentThreadPos);
    std::vector<UnsignedInt> GetRandomParticlesVersion2(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants, const CurrentThreadPosType& CurrentThreadPos);
    std::vector<UnsignedInt> GetRandomParticlesVersion3(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants, const CurrentThreadPosType& CurrentThreadPos);
public:
    void FindAndExecuteRandomReaction(UnsignedInt MaxNumberOfReactantsParam, const CurrentThreadPosType& CurrentThreadPos);
    bool FindAndExecuteChosenReaction(UnsignedInt ReactionId, const CurrentThreadPosType& CurrentThreadPos);
protected:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants, const CurrentThreadPosType& CurrentThreadPos) override;
    bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos) override;
    bool MakeChemicalReaction(ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos) override;
public:
    explicit CellEngineSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineSimulationSpaceStatistics(), Particles(ParticlesParam)
    {
        SetMakeSimulationStepNumberZero();
        ConstructDataForMultiThreadedExecution();
    }
    ~CellEngineSimulationSpace() override = default;
};

#endif
