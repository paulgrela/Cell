
#include <set>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"

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

bool CellEngineVoxelSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    return MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(ParticleObject, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineVoxelSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfVoxels(ParticleKindObjectForProduct.ListOfVoxels, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam);
}

template <class T>
void CellEngineVoxelSimulationSpace::CheckParticlesIndexes(ParticlesContainerInternal<T>& FormerParticlesIndexes, const string& FormerState)
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

template void CellEngineVoxelSimulationSpace::CheckParticlesIndexes(ParticlesContainerInternal<UniqueIdInt>& FormerParticlesIndexes, const string& FormerState);
template void CellEngineVoxelSimulationSpace::CheckParticlesIndexes(ParticlesContainerInternal<Particle>& FormerParticlesIndexes, const string& FormerState);

void CellEngineVoxelSimulationSpace::CheckCancelledParticlesIndexes()
{
    CheckParticlesIndexes<UniqueIdInt>(CancelledParticlesIndexes, "cancelled");
}

void CellEngineVoxelSimulationSpace::CheckFormerExistedParticlesIndexes()
{
    CheckParticlesIndexes<Particle>(FormerParticlesIndexes, "existed");
}
