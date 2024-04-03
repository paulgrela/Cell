#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include <stack>
#include <unordered_set>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
#include "CellEngineConfigData.h"
#include "CellEngineChemicalReactions.h"
#include "CellEngineChemicalReactionsCreator.h"
#include "CellEngineBasicVoxelsOperations.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineBasicParticlesOperations.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"
#include "CellEngineRealRandomParticlesGenerator.h"
#include "CellEngineNucleicAcidsBasicOperations.h"
#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"
#include "CellEngineTestParticlesGenerator.h"
#include "CellEngineChemicalReactionsInVoxelSpace.h"
#include "CellEngineNucleicAcidsChemicalReactionsInVoxelSpace.h"
#include "CellEngineVoxelSimulationSpaceStatistics.h"

#define SIMULATION_DETAILED_LOG

class CellEngineAllTypesOfParticlesGenerator : public CellEngineTestParticlesGenerator
{
protected:
    explicit CellEngineAllTypesOfParticlesGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineTestParticlesGenerator(ParticlesParam)
    {
    }
};

class CellEngineAllTypesOfChemicalReactionsInVoxelSpace : public CellEngineNucleicAcidsChemicalReactionsInVoxelSpace
{
protected:
    explicit CellEngineAllTypesOfChemicalReactionsInVoxelSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineNucleicAcidsChemicalReactionsInVoxelSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

class CellEngineVoxelSimulationSpaceForOuterClass : virtual public CellEngineBasicParticlesOperations
{
public:
    SimulationSpaceVoxel GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
    Particle& GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex);
};

class CellEngineVoxelSimulationSpace : public CellEngineChemicalReactionsCreator, public CellEngineAllTypesOfChemicalReactionsInVoxelSpace, public CellEngineAllTypesOfParticlesGenerator, public CellEngineVoxelSimulationSpaceForOuterClass, public CellEngineVoxelSimulationSpaceStatistics
{
private:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    void SetAtomInVoxelSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
public:
    void ClearVoxelSpaceAndParticles() override;
public:
    void ClearWholeVoxelSpace();
    void ClearSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], double MultiplyElectricChargeFactor);
public:
    void GenerateRandomReactionsForAllParticles();
    void GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
    void GenerateOneStepOfRandomReactionsForOneParticle(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam);
public:
    void GenerateRandomReactionForParticle(Particle& ParticleObject);
    void GenerateRandomReactionForSelectedVoxelSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateChosenReactionForSelectedVoxelSpace(UnsignedInt ReactionId, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    std::tuple<std::vector<std::pair<UniqueIdInt, UnsignedInt>>, bool> ChooseParticlesForReactionFromAllParticlesInProximity(const Reaction& ReactionObject);
protected:
    void PrepareRandomReaction();
    void FindAndExecuteRandomReaction();
    void FindAndExecuteChosenReaction(UnsignedInt ReactionId);
protected:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) override;
    bool IsChemicalReactionPossible(const Reaction& ReactionObject) override;
    bool MakeChemicalReaction(Reaction& ReactionObject) override;
public:
    explicit CellEngineVoxelSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam);
    ~CellEngineVoxelSimulationSpace();
};

#endif
