
#include <set>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineFullAtomSimulationSpace.h"

#include "CellEngineDataBuilderForFullAtomSimulationSpace.h"

using namespace std;

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

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
        GetParticleFromIndex(ParticleIndex).ListOfAtoms.clear();

        // for (const auto& NewPointElement : ParticleKindObjectForProduct.ListOfVoxels)
        //     //FillParticleElementInSpace(ParticleIndex, { static_cast<UnsignedInt>(NewPointElement.X) + VectorX, static_cast<UnsignedInt>(NewPointElement.Y) + VectorY, static_cast<UnsignedInt>(NewPointElement.Z) + VectorZ });
        //     FillParticleElementInSpace(ParticleIndex, { NewPointElement.X + VectorX, NewPointElement.Y + VectorY, NewPointElement.Z + VectorZ });

        GetMinMaxCoordinatesForParticle<float, CellEngineAtom>(GetParticleFromIndex(ParticleIndex), &Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, false);
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
//         GetMinMaxOfCoordinates(PosX, PosY, PosZ, XMin, XMax, YMin, YMax, ZMin, ZMax);
//
//         if (GetSpaceVoxel(PosX, PosY, PosZ) == GetZeroSimulationSpaceVoxel())
//             SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, PosX, PosY, PosZ);
//     }
//     CATCH("setting atom in voxel simulation space")
// }

void CellEngineFullAtomSimulationSpace::ClearSelectedSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    // try
    // {
    //     SetValueToVoxelsForCuboidSelectedSpace(nullptr, 0, StartXPosParam, StartYPosParam, StartZPosParam, StepXParam, StepYParam, StepZParam, SizeXParam, SizeYParam, SizeZParam);
    // }
    // CATCH("clearing selected space")
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
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ, float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam)
{
    return MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(ParticleObject, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtoms(ParticleKindObjectForProduct.ListOfAtoms, ParticlesInSector, CurrentSectorPos, ParticleKindObjectForProduct.Radius, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters);
}
