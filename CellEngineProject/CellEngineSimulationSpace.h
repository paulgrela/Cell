#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include <stack>
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
    void GenerateOneStepOfDiffusionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
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
    void SaveReactionsStatisticsToFile();
public:
    void GenerateDiffusionForBigPartOfCellSpace(UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
public:
    void GenerateChosenReactionsForWholeCellSpace(UnsignedInt ReactionId, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
    void GenerateRandomReactionsForWholeCellSpace(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
    void GenerateRandomReactionsForAllParticles();
public:
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
    void GenerateOneStepOfRandomReactionsForOneParticle(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
    void GenerateRandomReactionForParticle(Particle& ParticleObject);
    void GenerateRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateChosenReactionForSelectedSpace(UnsignedInt ReactionId, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
protected:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject);
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
    //void FindAndExecuteChosenReaction(UnsignedInt ReactionId);
    bool FindAndExecuteChosenReaction(UnsignedInt ReactionId);
protected:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants) override;
    bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject) override;
    bool MakeChemicalReaction(ChemicalReaction& ReactionObject) override;
public:
    explicit CellEngineSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : Particles(ParticlesParam)
    {
        SetMakeSimulationStepNumberZero();
    }
};

#endif
