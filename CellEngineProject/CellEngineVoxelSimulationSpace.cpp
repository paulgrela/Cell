
#include <set>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"

#include "CellEngineDataBuilderForVoxelSimulationSpace.h"

using namespace std;

SimulationSpaceVoxel CellEngineVoxelSimulationSpace::GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z)
{
    return GetSpaceVoxel(X, Y, Z);
}

Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

[[nodiscard]] float CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
{
    return static_cast<float>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForSimulationSpace;
};

[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForSimulationSpace) + (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2);
};

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, bool GetMemoryForVoxelSpace, ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInVoxelSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineTestParticlesInVoxelSpaceGenerator(ParticlesParam)
{
    try
    {
        CurrentThreadIndex = ThreadIndexParam;
        CurrentThreadPos = CurrentThreadPosParam;

        LoggersManagerObject.Log(STREAM("CurrentThreadIndexes = " << CurrentThreadIndex << " (" << CurrentThreadPos.ThreadPosX << "," << CurrentThreadPos.ThreadPosY << "," << CurrentThreadPos.ThreadPosZ << ")" << std::endl));

        if (GetMemoryForVoxelSpace == true)
        {
            SpacePointer = malloc(sizeof(Space_2048_2048_2048));

            SetStartValuesForSpaceMinMax();

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
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMin) << " ][ Xmax = " << to_string(XMax) << " ][ Ymin = " << to_string(YMin) << " ][ Ymax = " << to_string(YMax) << " ][ Zmin = " << to_string(ZMin) << " ][ Zmax = " << to_string(XMax) << " ] " << endl;
    return ss;
}

Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

void CellEngineVoxelSimulationSpace::SetAtomInVoxelSimulationSpace(const UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom)
{
    try
    {
        UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates(PosX, PosY, PosZ, XMin, XMax, YMin, YMax, ZMin, ZMax);

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
        Particles.clear();
        ClearWholeVoxelSpace();
    }
    CATCH("clearing voxel space and particles")
}

void CellEngineVoxelSimulationSpace::FillParticleElementInSpace(const UniqueIdInt ParticleIndex, const vector3_64 NewPointElement)
{
    try
    {
        SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, NewPointElement.X, NewPointElement.Y, NewPointElement.Z);
    }
    CATCH("filling particle element in space")
}

bool CellEngineVoxelSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    return MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(ParticleObject, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineVoxelSimulationSpace::MoveParticleByVectorIfSpaceIsEmpty(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ)
{
    return MoveParticleByVectorIfVoxelSpaceIsEmpty(ParticleObject, VectorX, VectorY, VectorZ);
}

bool CellEngineVoxelSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfVoxels(ParticleKindObjectForProduct.ListOfVoxels, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam);
}
