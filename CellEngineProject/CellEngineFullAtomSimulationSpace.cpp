
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

CellEngineFullAtomSimulationSpace::CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, const bool GetMemoryForFullAtomSpace, const ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineRealRandomParticlesInFullAtomSpaceGenerator(ParticlesParam), CellEngineVoxelSimulationSpaceStatistics()
{
    try
    {
        CurrentThreadIndex = ThreadIndexParam;
        CurrentThreadPos = CurrentThreadPosParam;

        LoggersManagerObject.Log(STREAM("CurrentThreadIndexes = " << CurrentThreadIndex << " (" << CurrentThreadPos.ThreadPosX << "," << CurrentThreadPos.ThreadPosY << "," << CurrentThreadPos.ThreadPosZ << ")" << std::endl));

        // Genomes.resize(2);
        // GenomesLines.resize(2);
    }
    CATCH("execution of constructor of full atom simulation space")
}

CellEngineFullAtomSimulationSpace::~CellEngineFullAtomSimulationSpace()
{
    try
    {
    }
    CATCH("execution of destructor of full atom simulation space")
}

[[nodiscard]] stringstream CellEngineFullAtomSimulationSpace::PrintSpaceMinMaxValues() const
{
    stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMinGlobal) << " ][ Xmax = " << to_string(XMaxGlobal) << " ][ Ymin = " << to_string(YMinGlobal) << " ][ Ymax = " << to_string(YMaxGlobal) << " ][ Zmin = " << to_string(ZMinGlobal) << " ][ Zmax = " << to_string(XMaxGlobal) << " ] " << endl;
    return ss;
}

void CellEngineFullAtomSimulationSpace::FillParticleElementsInSpace(const UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
{
    try
    {
        GetParticleFromIndex(ParticleIndex).ListOfAtoms.clear();

        for (const auto& NewPointElement : ParticleKindObjectForProduct.ListOfAtoms)
            GetParticleFromIndex(ParticleIndex).ListOfAtoms.emplace_back(NewPointElement);

        GetMinMaxCoordinatesForParticle<RealType, CellEngineAtom>(GetParticleFromIndex(ParticleIndex), &Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, false);
    }
    CATCH("filling particle elements in space")
}

void CellEngineFullAtomSimulationSpace::FillParticleElementInSpace(const UniqueIdInt ParticleIndex, const vector3_Real32 NewPointElement)
{
}

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForGenerator(const UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

void CellEngineFullAtomSimulationSpace::ClearSelectedSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
}

void CellEngineFullAtomSimulationSpace::ClearWholeFullAtomSpace()
{
    try
    {
    }
    CATCH("clearing whole full atom space")
};

void CellEngineFullAtomSimulationSpace::ClearFullAtomSpaceAndParticles()
{
    try
    {
        GetParticles().clear();
    }
    CATCH("clearing full atom space and particles")
}

SimulationSpaceSectorBounds CellEngineFullAtomSimulationSpace::GetBoundsForThreadSector()
{
    SimulationSpaceSectorBounds SimulationSpaceSectorBoundsObject{};

    try
    {
        if (CurrentThreadIndex == 0)
            SimulationSpaceSectorBoundsObject.SetParameters(-1 * CellEngineConfigDataObject.ShiftCenterX, -1 * CellEngineConfigDataObject.ShiftCenterY, -1 * CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ);
        else
            SimulationSpaceSectorBoundsObject.SetParametersForParallelExecutionSectors(CurrentThreadPos, CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);
    }
    CATCH("getting bounds for thread sector")

    return SimulationSpaceSectorBoundsObject;
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    return MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(ParticleObject, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtoms(ParticleKindObjectForProduct.ListOfAtoms, ParticlesInSector, CurrentSectorPos, ParticleKindObjectForProduct.Radius, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters);
}
