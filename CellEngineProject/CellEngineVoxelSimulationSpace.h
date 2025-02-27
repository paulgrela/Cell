#pragma once

#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_H

#include <stack>
#include <unordered_set>

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineConfigData.h"
#include "CellEngineChemicalReactionsEngine.h"
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
#include "CellEngineChemicalReactionsInVoxelSimulationSpace.h"
#include "CellEngineSimulationSpace.h"

class CellEngineVoxelSimulationSpace : public CellEngineSimulationSpace, public CellEngineTestParticlesInVoxelSpaceGenerator, public CellEngineChemicalReactionsInVoxelSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
{
protected:
    UnsignedInt XMinGlobal{}, XMaxGlobal{}, YMinGlobal{}, YMaxGlobal{}, ZMinGlobal{}, ZMaxGlobal{};
public:
    SimulationSpaceVoxel GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
    Particle& GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex);
private:
    ParticlesContainer<Particle>& Particles;
public:
    [[nodiscard]] static RealType ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    void SetAtomInVoxelSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
    Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) override;
public:
    void ClearVoxelSpaceAndParticles() override;
public:
    void ClearWholeVoxelSpace();
    void ClearSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void FillParticleElementsInSpace(UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, RealType VectorX, RealType VectorY, RealType VectorZ) override;
    void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_Real32 NewPointElement) override;
    bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, RealType VectorX, RealType VectorY, RealType VectorZ, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) override;
    bool CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, RealType VectorX, RealType VectorY, RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam) override;
protected:
    SimulationSpaceSectorBounds GetBoundsForThreadSector() override;
public:
    void GenerateNStepsOfDiffusionForWholeCellSpace(bool InBounds, RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps) override;
    void GenerateNStepsOfOneRandomReactionForWholeCellSpace(RealType XStartParam, RealType YStartParam, RealType ZStartParam, RealType XStepParam, RealType YStepParam, RealType ZStepParam, RealType XSizeParam, RealType YSizeParam, RealType ZSizeParam, RealType NumberOfSimulationSteps) override;
    void GenerateOneStepOfDiffusionForSelectedSpace(bool InBounds, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam) override;
    void GenerateOneRandomReactionForSelectedSpace(RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam, bool FindParticlesInProximityBool) override;
public:
    void GenerateOneStepOfElectricDiffusionForSelectedSpace(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, UnsignedInt AdditionalSpaceBoundFactor, double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    void GenerateNStepsOfDiffusionForBigPartOfCellSpace(bool InBounds, UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(UnsignedInt SizeNMultiplyFactor, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateOneChosenReactionForSelectedSpace(RealType ReactionId, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam);
    void GenerateNStepsOfOneChosenReactionForWholeCellSpace(UnsignedInt ReactionId, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt NumberOfSimulationSteps);
    void GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
public:
    template <class T>
    static void CheckParticlesIndexes(ParticlesDetailedContainer<T>& FormerParticlesIndexes, const string& FormerState);
    void CheckCancelledParticlesIndexes();
    void CheckFormerExistedParticlesIndexes();
public:
    explicit CellEngineVoxelSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, bool GetMemoryForVoxelSpace, ThreadIdType ThreadIndexParam, const CurrentThreadPosType& CurrentThreadPos);
    ~CellEngineVoxelSimulationSpace() override;
};

#endif
