
#include <set>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineFullAtomSimulationSpace.h"

#include "CellEngineDataBuilderForFullAtomSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsManager.h"

using namespace std;

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

CellEngineFullAtomSimulationSpace::CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, const bool GetMemoryForFullAtomSpace, const ThreadIdType ThreadIndexParam, CurrentThreadPosType CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineGenomeNucleicAcidsParticlesInFullAtomSpaceGenerator(ParticlesParam), CellEngineVoxelSimulationSpaceStatistics()
{
    try
    {
        CurrentMPIProcessIndex = ThreadIndexParam;

        CurrentThreadIndex = ThreadIndexParam;
        CurrentThreadPos = CurrentThreadPosParam;

        LoggersManagerObject.Log(STREAM("CurrentThreadIndexes = " << CurrentThreadIndex << " (" << CurrentThreadPos.ThreadPosX << "," << CurrentThreadPos.ThreadPosY << "," << CurrentThreadPos.ThreadPosZ << ")" << std::endl));

        Genomes.resize(2);
        GenomesLines.resize(2);
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

        if (ParticleKindObjectForProduct.ListOfAtoms.empty() == false)
            for (const auto& NewPointElement : ParticleKindObjectForProduct.ListOfAtoms)
                GetParticleFromIndex(ParticleIndex).ListOfAtoms.emplace_back(NewPointElement.X + VectorX, NewPointElement.Y + VectorY, NewPointElement.Z + VectorZ);
        else
            GetParticleFromIndex(ParticleIndex).ListOfAtoms.emplace_back(VectorX, VectorY, VectorZ);

        GetParticleFromIndex(ParticleIndex).Radius = ParticleKindObjectForProduct.Radius;
        GetParticleFromIndex(ParticleIndex).Center = { VectorX, VectorY, VectorZ };
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
    return ActualSimulationSpaceSectorBoundsObject;
}

void CellEngineFullAtomSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpace(const bool InBounds, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-10, 10);

        for (auto& ParticleInProximityObject : Particles[StartXPosParam][StartYPosParam][StartZPosParam].Particles)
        {
            CurrentSectorPos = SectorPosType{ static_cast<UnsignedInt>(StartXPosParam), static_cast<UnsignedInt>(StartYPosParam), static_cast<UnsignedInt>(StartZPosParam) };
            if (CellEngineUseful::IsDNA(ParticleInProximityObject.second.EntityId) == false)
                MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleInProximityObject.second, Particles, CurrentSectorPos, UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), 0, 0, 0, SizeXParam, SizeYParam, SizeZParam);
        }
    }
    CATCH("generating one step of diffusion for selected space")
}

void CellEngineFullAtomSimulationSpace::GenerateNStepsOfDiffusionForWholeCellSpace(const bool InBounds, const RealType XStartParam, const RealType YStartParam, const RealType ZStartParam, const RealType XStepParam, const RealType YStepParam, const RealType ZStepParam, const RealType XSizeParam, RealType YSizeParam, const RealType ZSizeParam, const RealType NumberOfSimulationSteps)
{
    try
    {
        CellEngineExecutionTimeStatisticsObject.ZeroMeasureTime();

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            FOR_EACH_PARTICLE_IN_XYZ_ONLY
                GenerateOneStepOfDiffusionForSelectedSpace(InBounds, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, XSizeParam, YSizeParam, ZSizeParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating random reactions in whole cell space has taken time: ","Execution in threads")));

        CellEngineExecutionTimeStatisticsObject.PrintMeasureTime();
    }
    CATCH("generating diffusion for whole cell space full atom")
}

void CellEngineFullAtomSimulationSpace::GenerateOneRandomReactionForSelectedSpace(const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const bool FindParticlesInProximityBool)
{
    try
    {
        PrepareRandomReaction();

        SetCurrentSectorPos(SectorPosType{ static_cast<UnsignedInt>(StartXPosParam), static_cast<UnsignedInt>(StartYPosParam), static_cast<UnsignedInt>(StartZPosParam) });

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        ActualSimulationSpaceSectorBoundsObject.SetParametersForChosenSector(StartXPosParam, StartYPosParam, StartZPosParam, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.SizeOfParticlesSectorX, CellEngineConfigDataObject.SizeOfParticlesSectorY, CellEngineConfigDataObject.SizeOfParticlesSectorZ);

        FindAndExecuteRandomReaction(min(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
    }
    CATCH("generating random reaction for selected space")
}

void CellEngineFullAtomSimulationSpace::GenerateNStepsOfOneRandomReactionForWholeCellSpace(const RealType XStartParam, const RealType YStartParam, const RealType ZStartParam, const RealType XStepParam, const RealType YStepParam, const RealType ZStepParam, const RealType XSizeParam, RealType YSizeParam, const RealType ZSizeParam, const RealType NumberOfSimulationSteps)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = false;

        CellEngineExecutionTimeStatisticsObject.ZeroMeasureTime();

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            FOR_EACH_PARTICLE_IN_XYZ_ONLY
                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, XSizeParam, YSizeParam, ZSizeParam, false);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating random reactions in whole cell space has taken time: ","Execution in threads")));

        CellEngineExecutionTimeStatisticsObject.PrintMeasureTime();
    }
    CATCH("generating random reactions for whole cell space full atom")
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    return MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(ParticleObject, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtoms(ParticleKindObjectForProduct.ListOfAtoms, ParticlesInSector, CurrentSectorPos, ParticleKindObjectForProduct.Radius, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters);
}

bool CellEngineFullAtomSimulationSpace::CheckInsertOfParticle(const MPIParticleSenderStruct& MPIParticleSenderToInsert)
{
    auto SimulationSpaceSectorBoundsObject = SimulationSpaceSectorBounds().SetParametersForChosenSector(MPIParticleSenderToInsert.SectorPos.X, MPIParticleSenderToInsert.SectorPos.Y, MPIParticleSenderToInsert.SectorPos.Z, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.SizeOfParticlesSectorX, CellEngineConfigDataObject.SizeOfParticlesSectorY, CellEngineConfigDataObject.SizeOfParticlesSectorZ);
    auto& ParticleKindToCheck = ParticlesKindsManagerObject.GetParticleKind(MPIParticleSenderToInsert.ParticleKindId);

    if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(ParticleKindToCheck, Particles, { MPIParticleSenderToInsert.SectorPos.X, MPIParticleSenderToInsert.SectorPos.Y, MPIParticleSenderToInsert.SectorPos.Z }, MPIParticleSenderToInsert.NewPosition.X, MPIParticleSenderToInsert.NewPosition.Y, MPIParticleSenderToInsert.NewPosition.Z, SimulationSpaceSectorBoundsObject) == true)
    {
        UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), MPIParticleSenderToInsert.ParticleKindId, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillParticleElementsInSpace(ParticleIndex, ParticleKindToCheck, MPIParticleSenderToInsert.NewPosition.X, MPIParticleSenderToInsert.NewPosition.Y, MPIParticleSenderToInsert.NewPosition.Z);

        return true;
    }

    return false;
}

void CellEngineFullAtomSimulationSpace::WriteNumberOfParticlesInEachSectorToFile() const
{
    try
    {
        UnsignedInt SectorsHistogram[1000];
        for (auto& SectorHistogram : SectorsHistogram)
            SectorHistogram = 0;

        FOR_EACH_PARTICLE_IN_XYZ_ONLY
        {
            if (Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size() > 1000)
                LoggersManagerObject.Log(STREAM("SECTOR = (" << ParticleSectorXIndex << "," << ParticleSectorYIndex  << "," << ParticleSectorZIndex << ") SIZE = " << Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size()));

            SectorsHistogram[Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size()]++;
        }

        UnsignedInt NumberOfSector = 0;
        UnsignedInt NumberOfAllSectorsWithParticles = 0;
        for (const auto& SectorHistogram : SectorsHistogram)
        {
            LoggersManagerObject.Log(STREAM("SECTOR SIZE = " << NumberOfSector << " NUMBER OF SECTORS OF THIS SIZE = " << SectorHistogram));
            NumberOfAllSectorsWithParticles += SectorHistogram;
            NumberOfSector++;
        }
        LoggersManagerObject.Log(STREAM("NUMBER OF ALL SECTORS = " << NumberOfAllSectorsWithParticles));
    }
    CATCH("writing number of particles in each sector to file")
}

void CellEngineFullAtomSimulationSpace::WriteNumberOfParticlesKindsWithoutAtoms()
{
    try
    {
        UnsignedInt CountParticleKindsWithoutAtoms = 0;

        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.ListOfAtoms.empty() == true)
            {
                LoggersManagerObject.Log(STREAM("Particle Kind Without Atoms = " << ParticleKindObject.second.EntityId << " " << ParticleKindObject.second.IdStr));
                CountParticleKindsWithoutAtoms++;
            }

        LoggersManagerObject.Log(STREAM("NUMBER OF PARTICLE KINDS WITHOUT VOXELS = " << CountParticleKindsWithoutAtoms << " " << ParticlesKindsManagerObject.ParticlesKinds.size()));
    }
    CATCH("getting number of particles kinds without atoms")
}