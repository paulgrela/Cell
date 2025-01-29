
#include <set>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineFullAtomSimulationSpace.h"

#include "CellEngineDataBuilderForFullAtomSimulationSpace.h"

using namespace std;

SimulationSpaceVoxel CellEngineFullAtomSimulationSpace::GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z)
{
    //return GetSpaceVoxel(X, Y, Z);
    return 1;
}

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

// [[nodiscard]] float CellEngineFullAtomSimulationSpace::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
// {
//     return static_cast<float>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForSimulationSpace;
// };
//
// [[nodiscard]] UnsignedInt CellEngineFullAtomSimulationSpace::ConvertToSpaceCoordinate(double CoordinateParam)
// {
//     return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForSimulationSpace) + (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2);
// };

// ZAMIENIC PONIZEJ NA CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesParam)
CellEngineFullAtomSimulationSpace::CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, const bool GetMemoryForVoxelSpace, const ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineTestParticlesInVoxelSpaceGenerator(ParticlesParam)
{
    try
    {
        CurrentThreadIndex = ThreadIndexParam;
        CurrentThreadPos = CurrentThreadPosParam;

        LoggersManagerObject.Log(STREAM("CurrentThreadIndexes = " << CurrentThreadIndex << " (" << CurrentThreadPos.ThreadPosX << "," << CurrentThreadPos.ThreadPosY << "," << CurrentThreadPos.ThreadPosZ << ")" << std::endl));

        Genomes.resize(2);
        GenomesLines.resize(2);
    }
    CATCH("execution of constructor of voxel simulation space")
}

CellEngineFullAtomSimulationSpace::~CellEngineFullAtomSimulationSpace()
{
    try
    {
        //free(SpacePointer);
    }
    CATCH("execution of destructor of voxel simulation space")
}

[[nodiscard]] stringstream CellEngineFullAtomSimulationSpace::PrintSpaceMinMaxValues() const
{
    stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMinGlobal) << " ][ Xmax = " << to_string(XMaxGlobal) << " ][ Ymin = " << to_string(YMinGlobal) << " ][ Ymax = " << to_string(YMaxGlobal) << " ][ Zmin = " << to_string(ZMinGlobal) << " ][ Zmax = " << to_string(XMaxGlobal) << " ] " << endl;
    return ss;
}

void CellEngineFullAtomSimulationSpace::FillParticleElementsInSpace(const UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ)
{
    try
    {
        GetParticleFromIndex(ParticleIndex).ListOfVoxels.clear();

        for (const auto& NewPointElement : ParticleKindObjectForProduct.ListOfVoxels)
            //FillParticleElementInSpace(ParticleIndex, { static_cast<UnsignedInt>(NewPointElement.X) + VectorX, static_cast<UnsignedInt>(NewPointElement.Y) + VectorY, static_cast<UnsignedInt>(NewPointElement.Z) + VectorZ });
            FillParticleElementInSpace(ParticleIndex, { NewPointElement.X + VectorX, NewPointElement.Y + VectorY, NewPointElement.Z + VectorZ });

        GetMinMaxCoordinatesForParticle<float>(GetParticleFromIndex(ParticleIndex), false);
    }
    CATCH("filling particle elements in space")
}

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForGenerator(const UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

// void CellEngineFullAtomSimulationSpace::SetAtomInFullAtomSimulationSpace(const UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom)
// {
//     try
//     {
//         UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
//         UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
//         UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);
//
//         GetMinMaxOfCoordinates(PosX, PosY, PosZ, XMin, XMax, YMin, YMax, ZMin, ZMax);
//
//         if (GetSpaceVoxel(PosX, PosY, PosZ) == GetZeroSimulationSpaceVoxel())
//             SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, PosX, PosY, PosZ);
//     }
//     CATCH("setting atom in voxel simulation space")
// }

void CellEngineFullAtomSimulationSpace::ClearSelectedSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        SetValueToVoxelsForCuboidSelectedSpace(nullptr, 0, StartXPosParam, StartYPosParam, StartZPosParam, StepXParam, StepYParam, StepZParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("clearing selected space")
}

void CellEngineFullAtomSimulationSpace::ClearWholeVoxelSpace()
{
    try
    {
        SetValueToVoxelsForCuboidSelectedSpace(nullptr, 0, 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension);
    }
    CATCH("clearing whole voxel space")
};

void CellEngineFullAtomSimulationSpace::ClearVoxelSpaceAndParticles()
{
    try
    {
        GetParticles().clear();
        ClearWholeVoxelSpace();
    }
    CATCH("clearing voxel space and particles")
}

void CellEngineFullAtomSimulationSpace::FillParticleElementInSpace(const UniqueIdInt ParticleIndex, const vector3_64 NewPointElement)
{
    try
    {
        SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, NewPointElement.X, NewPointElement.Y, NewPointElement.Z);
    }
    CATCH("filling particle element in space")
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    return MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(ParticleObject, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmpty(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ)
{
    return MoveParticleByVectorIfFullAtomSpaceIsEmpty(ParticleObject, VectorX, VectorY, VectorZ);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtoms(ParticleKindObjectForProduct.ListOfVoxels, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam);
}
