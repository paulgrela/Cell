#pragma once

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

#define SIMULATION_DETAILED_LOG

class CellEngineSimulationSpace : public CellEngineChemicalReactionsEngine, public CellEngineIllinoisDataCreator, public CellEngineCompiledDataCreator, virtual public CellEngineChemicalReactionsInSimulationSpace, public CellEngineSimulationSpaceStatistics, public CellEngineSimulationParallelExecutionManager
{
private:
    ParticlesContainer<Particle>& Particles;
protected:
    virtual void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, RealType VectorX, RealType VectorY, RealType VectorZ) = 0;
    virtual void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_Real32 NewPointElement) = 0;
protected:
    virtual bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, RealType VectorX, RealType VectorY, RealType VectorZ, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) = 0;
    virtual bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, RealType VectorX, RealType VectorY, RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) = 0;
protected:
    virtual SimulationSpaceSectorBounds GetBoundsForThreadSector() = 0;
public:
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateOneStepOfElectricDiffusionForOneParticle(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt ParticleIndexParam, ElectricChargeType (*NeighbourPoints)[3][3][3], UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedSpace(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(const Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], double MultiplyElectricChargeFactor);
public:
    void SetMakeSimulationStepNumberZero();
    void SetIncSimulationStepNumber();
    void SaveParticlesStatisticsOnce();
    void SaveReactionsStatisticsToFile() const;
public:
    void SaveHistogramOfParticlesStatisticsToFile();
    void SaveNumberOfParticlesStatisticsToFile();
public:
    void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
public:
    void GenerateNStepsOfDiffusionForBigPartOfCellSpace(bool InBounds, UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfDiffusionForWholeCellSpace(bool InBounds, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);

    void GenerateOneStepOfDiffusionForSelectedSpaceFullAtom(bool InBounds, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam);
    void GenerateNStepsOfDiffusionForWholeCellSpaceFullAtom(bool InBounds, RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps);

    void GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
public:
    void GenerateNStepsOfOneChosenReactionForWholeCellSpace(UnsignedInt ReactionId, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfOneRandomReactionForWholeCellSpace(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateOneStepOfRandomReactionsForAllParticles();
public:
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
    void GenerateOneStepOfRandomReactionsForOneParticleFromRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt ShiftIndexOfChosenParticle);
    void GenerateOneRandomReactionForChosenParticle(const Particle& ParticleObject);
    void GenerateOneRandomReactionForSelectedSpace(RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam, bool FindParticlesInProximityBool) override;
    void GenerateOneChosenReactionForSelectedSpace(RealType ReactionId, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam);
protected:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject);
protected:
    void PrepareRandomReaction();
    void FindAndExecuteRandomReactionVersion1(UnsignedInt MaxNumberOfReactantsParam);
    void FindAndExecuteRandomReactionVersion2(UnsignedInt MaxNumberOfReactantsParam);
    void FindAndExecuteRandomReactionVersion3(UnsignedInt MaxNumberOfReactantsParam);
    void FindAndExecuteRandomReactionVersion4(UnsignedInt MaxNumberOfReactantsParam);
private:
    set<UnsignedInt> GetAllPossibleReactionsFromParticlesInProximity();
protected:
    std::vector<UnsignedInt> GetRandomParticlesVersion1(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants);
    std::vector<UnsignedInt> GetRandomParticlesVersion2(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants);
    std::vector<UnsignedInt> GetRandomParticlesVersion3(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants) const;
public:
    void FindAndExecuteRandomReaction(UnsignedInt MaxNumberOfReactantsParam);
    bool FindAndExecuteChosenReaction(UnsignedInt ReactionId);
protected:
    bool CancelChemicalReaction(const vector<UniqueIdInt>& CreatedParticlesIndexes, chrono::high_resolution_clock::time_point start_time, const ParticleKind& ParticleKindObjectForProduct, char PlaceStr);
    bool PlaceProductParticleInSpaceInDeterminedPositionOrCancelReaction(UniqueIdInt ParticleIndex, const vector<UniqueIdInt>& CreatedParticlesIndexes, UnsignedInt CenterIndex, const ListOfAtomsType& Centers, ParticleKind& ParticleKindObjectForProduct, chrono::high_resolution_clock::time_point start_time);
    bool PlaceProductParticleInSpaceInRandomPositionOrCancelReaction(UniqueIdInt ParticleIndex, const vector<UniqueIdInt>& CreatedParticlesIndexes, UnsignedInt CenterIndex, const ListOfAtomsType& Centers, ParticleKind& ParticleKindObjectForProduct, chrono::high_resolution_clock::time_point start_time);
protected:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants) override;
    bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject) override;
    bool MakeChemicalReaction(ChemicalReaction& ReactionObject) override;
public:
    explicit CellEngineSimulationSpace(ParticlesContainer<Particle>& ParticlesParam) : CellEngineSimulationSpaceStatistics(), Particles(ParticlesParam), CellEngineSimulationParallelExecutionManager()
    {
        SetMakeSimulationStepNumberZero();
    }
    ~CellEngineSimulationSpace() override = default;
};

#endif
