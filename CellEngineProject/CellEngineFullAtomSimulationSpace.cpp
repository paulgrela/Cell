
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

CellEngineFullAtomSimulationSpace::CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, const bool GetMemoryForFullAtomSpace, const ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineRealRandomParticlesInFullAtomSpaceGenerator(ParticlesParam)
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

void CellEngineFullAtomSimulationSpace::FillParticleElementsInSpace(const UniqueIdInt ParticleIndex, ParticleKind& ParticleKindObjectForProduct, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ)
{
    try
    {
        GetParticleFromIndex(ParticleIndex).ListOfAtoms.clear();

        GetMinMaxCoordinatesForParticle<float, CellEngineAtom>(GetParticleFromIndex(ParticleIndex), &Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, false);
    }
    CATCH("filling particle elements in space")
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
        //SetValueToVoxelsForCuboidSelectedSpace(nullptr, 0, 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension);
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

void CellEngineFullAtomSimulationSpace::FillParticleElementInSpace(const UniqueIdInt ParticleIndex, const vector3_64 NewPointElement)
{
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, float VectorX, float VectorY, float VectorZ, float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam)
{
    //cout << "CCC" << endl;
    //std::cout << "Particle Drawn Inside " << ParticleObject.Index << " " << ParticleObject.EntityId << std::endl;

    return MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(ParticleObject, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtoms(ParticleKindObjectForProduct.ListOfAtoms, ParticlesInSector, CurrentSectorPos, ParticleKindObjectForProduct.Radius, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters);
}
