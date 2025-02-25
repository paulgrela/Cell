
#include <set>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"
#include "CellEngineChemicalReactionsManager.h"

#include "CellEngineDataBuilderForVoxelSimulationSpace.h"

using namespace std;

SimulationSpaceVoxel CellEngineVoxelSimulationSpace::GetSpaceVoxelForOuterClass(const UnsignedInt X, const UnsignedInt Y, const UnsignedInt Z)
{
    return GetSpaceVoxel(X, Y, Z);
}

Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndexForOuterClass(const UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

[[nodiscard]] RealType CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(const UnsignedInt CoordinateParam)
{
    return static_cast<RealType>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForSimulationSpace;
};

[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(const double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForSimulationSpace) + (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2);
};

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, const bool GetMemoryForVoxelSpace, const ThreadIdType ThreadIndexParam, const CurrentThreadPosType& CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInVoxelSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineTestParticlesInVoxelSpaceGenerator(ParticlesParam), CellEngineVoxelSimulationSpaceStatistics()
{
    try
    {
        CurrentThreadIndex = ThreadIndexParam;
        CurrentThreadPos = CurrentThreadPosParam;

        LoggersManagerObject.Log(STREAM("CurrentThreadIndexes = " << CurrentThreadIndex << " (" << CurrentThreadPos.ThreadPosX << "," << CurrentThreadPos.ThreadPosY << "," << CurrentThreadPos.ThreadPosZ << ")" << std::endl));

        if (GetMemoryForVoxelSpace == true)
        {
            SpacePointer = malloc(sizeof(Space_2048_2048_2048));

            XMinGlobal = YMinGlobal = ZMinGlobal = 10000;
            XMaxGlobal = YMaxGlobal = ZMaxGlobal = 0;

            SetValueToVoxelsForCuboidSelectedSpace(nullptr, GetZeroSimulationSpaceVoxel(), 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension);

            Genomes.resize(2);
            GenomesLines.resize(2);
        }
        else
        {
            SpacePointer = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SpacePointer;
        }
    }
    CATCH("execution of constructor of voxel simulation space")
}

CellEngineVoxelSimulationSpace::~CellEngineVoxelSimulationSpace()
{
    try
    {
        free(SpacePointer);
    }
    CATCH("execution of destructor of voxel simulation space")
}

[[nodiscard]] stringstream CellEngineVoxelSimulationSpace::PrintSpaceMinMaxValues() const
{
    stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMinGlobal) << " ][ Xmax = " << to_string(XMaxGlobal) << " ][ Ymin = " << to_string(YMinGlobal) << " ][ Ymax = " << to_string(YMaxGlobal) << " ][ Zmin = " << to_string(ZMinGlobal) << " ][ Zmax = " << to_string(XMaxGlobal) << " ] " << endl;
    return ss;
}

void CellEngineVoxelSimulationSpace::FillParticleElementsInSpace(const UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
{
    try
    {
        GetParticleFromIndex(ParticleIndex).ListOfVoxels.clear();

        for (const auto& NewPointElement : ParticleKindObjectForProduct.ListOfVoxels)
            FillParticleElementInSpace(ParticleIndex, { NewPointElement.X + VectorX, NewPointElement.Y + VectorY, NewPointElement.Z + VectorZ });

        GetMinMaxCoordinatesForParticle<UnsignedInt, vector3_16>(GetParticleFromIndex(ParticleIndex), &Particle::ListOfVoxels, &ParticleKind::ListOfVoxels, false);
    }
    CATCH("filling particle elements in space")
}

void CellEngineVoxelSimulationSpace::FillParticleElementInSpace(const UniqueIdInt ParticleIndex, const vector3_Real32 NewPointElement)
{
    try
    {
        SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, NewPointElement.X, NewPointElement.Y, NewPointElement.Z);
    }
    CATCH("filling particle element in space")
}

Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndexForGenerator(const UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

void CellEngineVoxelSimulationSpace::SetAtomInVoxelSimulationSpace(const UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom)
{
    try
    {
        const UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
        const UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        const UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates<UnsignedInt>(PosX, PosY, PosZ, XMinGlobal, XMaxGlobal, YMinGlobal, YMaxGlobal, ZMinGlobal, ZMaxGlobal);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == GetZeroSimulationSpaceVoxel())
            SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, PosX, PosY, PosZ);
    }
    CATCH("setting atom in voxel simulation space")
}

void CellEngineVoxelSimulationSpace::ClearSelectedSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        SetValueToVoxelsForCuboidSelectedSpace(nullptr, GetZeroSimulationSpaceVoxel(), StartXPosParam, StartYPosParam, StartZPosParam, StepXParam, StepYParam, StepZParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("clearing selected space")
}

void CellEngineVoxelSimulationSpace::ClearWholeVoxelSpace()
{
    try
    {
        SetValueToVoxelsForCuboidSelectedSpace(nullptr, GetZeroSimulationSpaceVoxel(), 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension);
    }
    CATCH("clearing whole voxel space")
};

void CellEngineVoxelSimulationSpace::ClearVoxelSpaceAndParticles()
{
    try
    {
        GetParticles().clear();
        ClearWholeVoxelSpace();
    }
    CATCH("clearing voxel space and particles")
}

SimulationSpaceSectorBounds CellEngineVoxelSimulationSpace::GetBoundsForThreadSector()
{
    SimulationSpaceSectorBounds SimulationSpaceSectorBoundsObject{};

    try
    {
        if (CurrentThreadIndex == 0)
            SimulationSpaceSectorBoundsObject.SetParameters(0, 0, 0, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension);
        else
            SimulationSpaceSectorBoundsObject.SetParametersForParallelExecutionSectors(CurrentThreadPos, CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);
    }
    CATCH("getting bounds for thread sector")

    return SimulationSpaceSectorBoundsObject;
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfElectricDiffusionForSelectedSpace(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        ElectricChargeType NeighbourPoints[3][3][3];

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        auto ParticlesSortedByCapacityFoundInProximityCopy(LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity);

        for (auto& ParticleInProximityIndex : ParticlesSortedByCapacityFoundInProximityCopy)
            if (CellEngineUseful::IsDNA(GetParticleFromIndex(ParticleInProximityIndex).EntityId) == false)
                GenerateOneStepOfElectricDiffusionForOneParticle(TypeOfLookingForParticles, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, ParticleInProximityIndex, &NeighbourPoints, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating one step of electric diffusion")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpace(const bool InBounds, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-1, 1);

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
        for (auto& ParticleInProximityIndex : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
            if (CellEngineUseful::IsDNA(GetParticleFromIndex(ParticleInProximityIndex).EntityId) == false)
                MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(GetParticleFromIndex(ParticleInProximityIndex), Particles, CurrentSectorPos, UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("generating one step of diffusion for selected space")
}

void CellEngineVoxelSimulationSpace::GenerateNStepsOfDiffusionForWholeCellSpace(const bool InBounds, const RealType XStartParam, const RealType YStartParam, const RealType ZStartParam, const RealType XStepParam, const RealType YStepParam, const RealType ZStepParam, const RealType XSizeParam, RealType YSizeParam, const RealType ZSizeParam, const RealType NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                        GenerateOneStepOfDiffusionForSelectedSpace(InBounds, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating diffusion for whole cell space")
}

void CellEngineVoxelSimulationSpace::GenerateOneRandomReactionForSelectedSpace(const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const bool FindParticlesInProximityBool)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityBool == false || (FindParticlesInProximityBool == true && FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true))
        {
            ActualSimulationSpaceSectorBoundsObject = { StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam - 1, StartYPosParam + SizeYParam - 1, StartZPosParam + SizeZParam - 1 };
            FindAndExecuteRandomReaction(min(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
        }
    }
    CATCH("generating random reaction for selected space")
}

void CellEngineVoxelSimulationSpace::GenerateNStepsOfOneRandomReactionForWholeCellSpace(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineExecutionTimeStatisticsObject.PrintMeasureTime();

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                        GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam, true);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating random reactions in whole cell space has taken time: ","Execution in threads")));

        CellEngineExecutionTimeStatisticsObject.PrintMeasureTime();
    }
    CATCH("generating n steps of random reactions for whole cell space")
}

void CellEngineVoxelSimulationSpace::GenerateOneChosenReactionForSelectedSpace(const RealType ReactionId, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    try
    {
        if (FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            if (ReactionId != 0)
            {
                ActualSimulationSpaceSectorBoundsObject = { StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam - 1, StartYPosParam + SizeYParam - 1, StartZPosParam + SizeZParam - 1 };
                FindAndExecuteChosenReaction(ReactionId);
            }
    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::GenerateNStepsOfOneChosenReactionForWholeCellSpace(const UnsignedInt ReactionId, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                        GenerateOneChosenReactionForSelectedSpace(ReactionId, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for whole cell space")
}








inline void GetRangeOfParticlesForRandomParticles(UniqueIdInt& StartParticleIndexParam, UniqueIdInt& EndParticleIndexParam, UniqueIdInt MaxParticleIndex)
{
    if (EndParticleIndexParam == 0)
    {
        StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
        EndParticleIndexParam = MaxParticleIndex;
    }
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-1, 1);

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(GetParticleFromIndex(ParticleIndex), Particles, CurrentSectorPos, UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("generating one step of diffusion for selected range of particles")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        ElectricChargeType NeighbourPoints[3][3][3];

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            GenerateOneStepOfElectricDiffusionForOneParticle(TypeOfLookingForParticles, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, ParticleIndex, &NeighbourPoints, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating one step of electric diffusion")
}

void CellEngineVoxelSimulationSpace::GenerateNStepsOfDiffusionForBigPartOfCellSpace(const bool InBounds, const UnsignedInt SizeNMultiplyFactor, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam - SizeNMultiplyFactor * XStepParam; PosX <= XStartParam + SizeNMultiplyFactor * XStepParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam - SizeNMultiplyFactor * YStepParam; PosY <= YStartParam + SizeNMultiplyFactor * YStepParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam - SizeNMultiplyFactor * ZStepParam; PosZ <= ZStartParam + SizeNMultiplyFactor * ZStepParam; PosZ += ZStepParam)
                        GenerateOneStepOfDiffusionForSelectedSpace(InBounds, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating diffusion for big part of cell")
}

void CellEngineVoxelSimulationSpace::GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(const UnsignedInt SizeNMultiplyFactor, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam - SizeNMultiplyFactor * XStepParam; PosX <= XStartParam + SizeNMultiplyFactor * XStepParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam - SizeNMultiplyFactor * YStepParam; PosY <= YStartParam + SizeNMultiplyFactor * YStepParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam - SizeNMultiplyFactor * ZStepParam; PosZ <= ZStartParam + SizeNMultiplyFactor * ZStepParam; PosZ += ZStepParam)
                        GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam, true);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for big part of cell")
}

























bool CellEngineVoxelSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    return MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(ParticleObject, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineVoxelSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfVoxels(ParticleKindObjectForProduct.ListOfVoxels, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam);
}

template <class T>
void CellEngineVoxelSimulationSpace::CheckParticlesIndexes(ParticlesDetailedContainer<T>& FormerParticlesIndexes, const string& FormerState)
{
    try
    {
        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension; PosZ++)
                {
                    SimulationSpaceVoxel SimulationSpaceVoxelObject = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSpaceVoxelForOuterClass(PosX, PosY, PosZ);
                    if (SimulationSpaceVoxelObject != CellEngineParticlesVoxelsOperations::GetZeroSimulationSpaceVoxel())
                        if (auto FoundParticleIter = CellEngineDataFileObjectPointer->GetParticleIteratorFromIndex(SimulationSpaceVoxelObject); FoundParticleIter == CellEngineDataFileObjectPointer->GetParticleEnd())
                            if (auto FormerParticleIter = FormerParticlesIndexes.find(SimulationSpaceVoxelObject); FormerParticleIter != FormerParticlesIndexes.end())
                                LoggersManagerObject.LogInformation(STREAM("Try to draw the particle from not existing index but formerly " << FormerState << " = " << SimulationSpaceVoxelObject));
                            else
                                LoggersManagerObject.LogInformation(STREAM("Try to draw the particle from not existing index but formerly not " << FormerState <<  " = " << SimulationSpaceVoxelObject));
                }
    }
    CATCH("checking particles indexes")
}

template void CellEngineVoxelSimulationSpace::CheckParticlesIndexes(ParticlesDetailedContainer<UniqueIdInt>& FormerParticlesIndexes, const string& FormerState);
template void CellEngineVoxelSimulationSpace::CheckParticlesIndexes(ParticlesDetailedContainer<Particle>& FormerParticlesIndexes, const string& FormerState);

void CellEngineVoxelSimulationSpace::CheckCancelledParticlesIndexes()
{
    CheckParticlesIndexes<UniqueIdInt>(CancelledParticlesIndexes, "cancelled");
}

void CellEngineVoxelSimulationSpace::CheckFormerExistedParticlesIndexes()
{
    CheckParticlesIndexes<Particle>(FormerParticlesIndexes, "existed");
}
