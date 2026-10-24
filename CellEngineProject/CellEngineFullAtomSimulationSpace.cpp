
#include <unordered_map>

#include "CellEngineAtom.h"
#include "CellEngineUseful.h"
#include "CellEngineMPIProcess.h"
#include "CellEngineSimulationExecutionMeasureTimeStatisticsData.h"
#include "CellEngineFullAtomSimulationSpace.h"

#include "CellEngineDataBuilderForFullAtomSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsManager.h"

using namespace std;

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForOuterClass(UniqueIdUInt64Type ParticleUniqueId)
{
    return GetParticleFromUniqueId(ParticleUniqueId);
}

void CellEngineFullAtomSimulationSpace::CountMinMaxValuesOfBordersAtomsPositionsInTheCell()
{
    try
    {
        XMinGlobal = YMinGlobal = ZMinGlobal = 10000;
        XMaxGlobal = YMaxGlobal = ZMaxGlobal = -10000;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles | views::values)
                for (const auto& AtomObject : ParticleObject.ListOfAtoms)
                    GetMinMaxOfCoordinates<RealType>(AtomObject.X, AtomObject.Y, AtomObject.Z, XMinGlobal, XMaxGlobal, YMinGlobal, YMaxGlobal, ZMinGlobal, ZMaxGlobal);
    }
    CATCH("counting min max values of borders atoms positions in the cell")
}

CellEngineFullAtomSimulationSpace::CellEngineFullAtomSimulationSpace(ParticlesContainer<Particle>& ParticlesParam, const bool GetMemoryForFullAtomSpace, const ThreadIdType ThreadIndexParam, const ThreadPosType& CurrentThreadPosParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineChemicalReactionsInFullAtomSimulationSpace(ParticlesParam), CellEngineSimulationSpace(ParticlesParam), CellEngineGenomeNucleicAcidsParticlesInFullAtomSpaceGenerator(ParticlesParam), CellEngineVoxelSimulationSpaceStatistics()
{
    try
    {
        XMinGlobal = YMinGlobal = ZMinGlobal = 10000;
        XMaxGlobal = YMaxGlobal = ZMaxGlobal = -10000;

        MPIProcessIndex = ThreadIndexParam;

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
    stringstream SpaceMinMaxValuesStringStream;
    SpaceMinMaxValuesStringStream << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMinGlobal) << " ][ Xmax = " << to_string(XMaxGlobal) << " ][ Ymin = " << to_string(YMinGlobal) << " ][ Ymax = " << to_string(YMaxGlobal) << " ][ Zmin = " << to_string(ZMinGlobal) << " ][ Zmax = " << to_string(XMaxGlobal) << " ] " << endl;
    return SpaceMinMaxValuesStringStream;
}

void CellEngineFullAtomSimulationSpace::FillParticleElementsInSpace(const UniqueIdUInt64Type ParticleUniqueId, const ParticleKind& ParticleKindObjectForProduct, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
{
    try
    {
        GetParticleFromUniqueId(ParticleUniqueId).ListOfAtoms.clear();

        if (ParticleKindObjectForProduct.ListOfAtoms.empty() == false)
            for (const auto& NewPointElement : ParticleKindObjectForProduct.ListOfAtoms)
                GetParticleFromUniqueId(ParticleUniqueId).ListOfAtoms.emplace_back(NewPointElement.X + VectorX, NewPointElement.Y + VectorY, NewPointElement.Z + VectorZ);
        else
            GetParticleFromUniqueId(ParticleUniqueId).ListOfAtoms.emplace_back(VectorX, VectorY, VectorZ);

        GetParticleFromUniqueId(ParticleUniqueId).Radius = ParticleKindObjectForProduct.Radius;
        GetParticleFromUniqueId(ParticleUniqueId).Center = { .X = VectorX, .Y = VectorY, .Z = VectorZ };

        #if defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION) || defined(DEBUG_LOG_CHEMICAL_REACTION)
        if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false)
            CellEngineUseful::CheckCenterForSector<Particle>(GetParticleFromUniqueId(ParticleUniqueId), Particles, "CF1");
        #endif
    }
    CATCH("filling particle elements in space")
}

void CellEngineFullAtomSimulationSpace::FillParticleElementInSpace(const UniqueIdUInt64Type ParticleUniqueId, const vector3_Real32 NewPointElement)
{
}

Particle& CellEngineFullAtomSimulationSpace::GetParticleFromIndexForGenerator(const UniqueIdUInt64Type ParticleUniqueId)
{
    return GetParticleFromUniqueId(ParticleUniqueId);
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

void CellEngineFullAtomSimulationSpace::TransferParticlesFromOneSectorToAnotherSectorInOneStepOfDiffusionForSelectedSpaceForParallelExecutionOrForContainersForSpeedAndForNonParallelExecution(vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, UnsignedInt& NumberOfParticlesMovedInDiffusionParam) const
{
    try
    {
        NumberOfParticlesMovedInDiffusionParam += ListOfParticlesToChangeSectors.size();

        for (const auto& [ParticleUniqueId, SectorPosSource, SectorPosTarget] : ListOfParticlesToChangeSectors)
            Particles[SectorPosTarget.SectorPosX][SectorPosTarget.SectorPosY][SectorPosTarget.SectorPosZ].Particles.insert(Particles[SectorPosSource.SectorPosX][SectorPosSource.SectorPosY][SectorPosSource.SectorPosZ].Particles.extract(ParticleUniqueId));

        ListOfParticlesToChangeSectors.clear();
    }
    CATCH("")
}

void CellEngineFullAtomSimulationSpace::TransferParticlesFromOneSectortoAnotherSectorInOneStepOfDiffusionForSelectedSpaceForContainersNotForSpeedAndForNonParallelExecution(std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors) const
{
    try
    {
        for (const auto& [ParticleUniqueId, SectorPosSource, SectorPosTarget] : ListOfParticlesToChangeSectors)
        {
            auto& LocalSourceParticles = Particles[SectorPosSource.SectorPosX][SectorPosSource.SectorPosY][SectorPosSource.SectorPosZ].Particles;

            if (const auto ParticleFromSourceToMoveToTargetIterator = LocalSourceParticles.find(ParticleUniqueId); ParticleFromSourceToMoveToTargetIterator != LocalSourceParticles.end())
            {
                Particles[SectorPosTarget.SectorPosX][SectorPosTarget.SectorPosY][SectorPosTarget.SectorPosZ].Particles.insert_or_assign(ParticleFromSourceToMoveToTargetIterator->first, std::move(ParticleFromSourceToMoveToTargetIterator->second));

                DEBUG_LOG_DIFFUSION(CellEngineUseful::CheckCenterForSector<Particle>(ParticleFromSourceToMoveToTargetIterator->second, Particles, "C2");)

                LocalSourceParticles.erase(ParticleFromSourceToMoveToTargetIterator);
            }
            else
            {
                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("LACK OF PARTICLE INDEX IN SOURCE MAP DATA = " << ParticleUniqueId << " TARGET = " << SectorPosTarget.SectorPosX << "," << SectorPosTarget.SectorPosY << "," << SectorPosTarget.SectorPosZ << " FROM " << SectorPosSource.SectorPosX << "," << SectorPosSource.SectorPosY << "," << SectorPosSource.SectorPosZ));)
            }
        }

        ListOfParticlesToChangeSectors.clear();
    }
    CATCH("")
}

void CellEngineFullAtomSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpace(vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, UnsignedInt& NumberOfParticlesMovedInDiffusionParam, UnsignedInt& NumberOfParticlesToBeMovedInDiffusionRejectedLocallyParam, UnsignedInt& NumberOfParticlesSentToCheckInAnotherThreadToBeMovedInDiffusion, UnsignedInt& NumberOfDNAParticlesParam, const bool InBounds, const UnsignedInt StartSectorXPosParam, const UnsignedInt StartSectorYPosParam, const UnsignedInt StartSectorZPosParam, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const bool ContainersNotForSpeedUnparalled)
{
    try
    {
        uniform_real_distribution<RealType> UniformDistributionObjectMoveParticleDirection_RealType(-ParticleMoveSizeInDiffusionRealType, ParticleMoveSizeInDiffusionRealType);

        for (auto& ParticleInProximityObject : Particles[StartSectorXPosParam][StartSectorYPosParam][StartSectorZPosParam].Particles | views::values)
                                                                                                                        //if (ParticleInProximityObject.ParticleKindId != 100042) //gdy reakcje wlaczone
        {
            DEBUG_LOG_DIFFUSION(CellEngineUseful::CheckCenterForSector<Particle>(ParticleInProximityObject, Particles, "C0");)

            SetCurrentSectorPos(SectorPosType{ .SectorPosX = static_cast<SignedInt>(StartSectorXPosParam), .SectorPosY = static_cast<SignedInt>(StartSectorYPosParam), .SectorPosZ = static_cast<SignedInt>(StartSectorZPosParam) });

            if (CellEngineUseful::IsDNA(ParticleInProximityObject.ParticleKindId) == false)
            {
                if (MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleInProximityObject, Particles, ListOfParticlesToChangeSectors, GetRandomValue<uniform_real_distribution, RealType>(UniformDistributionObjectMoveParticleDirection_RealType), GetRandomValue<uniform_real_distribution, RealType>(UniformDistributionObjectMoveParticleDirection_RealType), GetRandomValue<uniform_real_distribution, RealType>(UniformDistributionObjectMoveParticleDirection_RealType), StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == false)
                    NumberOfParticlesToBeMovedInDiffusionRejectedLocallyParam++;
            }
            else
                NumberOfDNAParticlesParam++;
        }

        if (ContainersNotForSpeedUnparalled == false)
            TransferParticlesFromOneSectorToAnotherSectorInOneStepOfDiffusionForSelectedSpaceForParallelExecutionOrForContainersForSpeedAndForNonParallelExecution(ListOfParticlesToChangeSectors, NumberOfParticlesMovedInDiffusionParam);
        else
        {
            #ifdef CONTAINERS_FOR_SPEED
            TransferParticlesFromOneSectorToAnotherSectorInOneStepOfDiffusionForSelectedSpaceForParallelExecutionOrForContainersForSpeedAndForNonParallelExecution(ListOfParticlesToChangeSectors, NumberOfParticlesMovedInDiffusionParam);
            #else
            TransferParticlesFromOneSectorToAnotherSectorInOneStepOfDiffusionForSelectedSpaceForContainersNotForSpeedAndForNonParallelExecution(ListOfParticlesToChangeSectors, InBounds, StartSectorXPosParam, StartSectorYPosParam, StartSectorZPosParam, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, true);
            #endif
        }
    }
    CATCH("generating one step of diffusion for selected space")
}

void CellEngineFullAtomSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInMPIProcesses(const bool InBounds, const UnsignedInt StartSectorXPosParam, const UnsignedInt StartSectorYPosParam, const UnsignedInt StartSectorZPosParam, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    GenerateOneStepOfDiffusionForSelectedSpace(ListOfParticlesToChangeSectors, NumberOfParticlesMovedInDiffusion, NumberOfParticlesToBeMovedInDiffusionRejectedLocally, NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread, NumberOfDNAParticles, InBounds, StartSectorXPosParam, StartSectorYPosParam, StartSectorZPosParam, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, false);
}

void CellEngineFullAtomSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInThreads(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const bool InBounds, const UnsignedInt StartSectorXPosParam, const UnsignedInt StartSectorYPosParam, const UnsignedInt StartSectorZPosParam, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("INSIDE XStart = " << StartXPosParam << " YStart = " << StartYPosParam << " ZStart = " << StartZPosParam << " XEnd = " << SizeXParam << " YEnd = " << SizeYParam << " ZEnd = " << SizeZParam << " PosX = " << StartSectorXPosParam << " PosY = " << StartSectorXPosParam << " PosZ = " << StartSectorXPosParam));)

    GenerateOneStepOfDiffusionForSelectedSpace(CurrentThreadLocalSimulationSpaceData->ListOfParticlesToChangeSectors, CurrentThreadLocalSimulationSpaceData->NumberOfParticlesMovedInDiffusion, CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionRejectedLocally, CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread, CurrentThreadLocalSimulationSpaceData->NumberOfDNAParticles, InBounds, StartSectorXPosParam, StartSectorYPosParam, StartSectorZPosParam, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, false);
}

void CellEngineFullAtomSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpaceForExecutionNonParallel(const bool InBounds, const UnsignedInt StartSectorXPosParam, const UnsignedInt StartSectorYPosParam, const UnsignedInt StartSectorZPosParam, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    GenerateOneStepOfDiffusionForSelectedSpace(ListOfParticlesToChangeSectors, NumberOfParticlesMovedInDiffusion, NumberOfParticlesToBeMovedInDiffusionRejectedLocally, NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread, NumberOfDNAParticles, InBounds, StartSectorXPosParam, StartSectorYPosParam, StartSectorZPosParam, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, true);
}

void CellEngineFullAtomSimulationSpace::GenerateNStepsOfDiffusionForWholeCellSpace(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const bool InBounds, const RealType XStartParam, const RealType YStartParam, const RealType ZStartParam, const RealType XStepParam, const RealType YStepParam, const RealType ZStepParam, const RealType XSizeParam, RealType YSizeParam, const RealType ZSizeParam, const RealType NumberOfSimulationSteps)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = false;

        CellEngineExecutionTimeStatisticsObject.ZeroMeasureTime();

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            FOR_EACH_SECTOR_IN_XYZ_ONLY
                GenerateOneStepOfDiffusionForSelectedSpaceForExecutionNonParallel(InBounds, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, -CellEngineConfigDataObject.ShiftCenterX, -CellEngineConfigDataObject.ShiftCenterY, -CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating random diffusion in whole cell space has taken time: ","Execution in threads")));

        CellEngineExecutionTimeStatisticsObject.PrintMeasureTime();
    }
    CATCH("generating diffusion for whole cell space full atom")
}

void CellEngineFullAtomSimulationSpace::GenerateOneRandomReactionForSelectedSpace(const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const bool FindParticlesInProximityBool)
{
    try
    {
        PrepareRandomReaction();

        SetCurrentSectorPos(SectorPosType{ .SectorPosX = static_cast<SignedInt>(StartXPosParam), .SectorPosY = static_cast<SignedInt>(StartYPosParam), .SectorPosZ = static_cast<SignedInt>(StartZPosParam) });

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType);

        ActualSimulationSpaceSectorBoundsObject.SetParametersForChosenSector(StartXPosParam, StartYPosParam, StartZPosParam, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.SizeOfParticlesSectorX, CellEngineConfigDataObject.SizeOfParticlesSectorY, CellEngineConfigDataObject.SizeOfParticlesSectorZ);

        FindAndExecuteRandomReaction(min(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
    }
    CATCH("generating random reaction for selected space")
}

void CellEngineFullAtomSimulationSpace::GenerateNStepsOfOneRandomReactionForWholeCellSpace(const RealType XStartParam, const RealType YStartParam, const RealType ZStartParam, const RealType XStepParam, const RealType YStepParam, const RealType ZStepParam, const RealType XSizeParam, const RealType YSizeParam, const RealType ZSizeParam, const RealType NumberOfSimulationSteps)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = false;

        CellEngineExecutionTimeStatisticsObject.ZeroMeasureTime();

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            FOR_EACH_SECTOR_IN_XYZ_ONLY
                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, XSizeParam, YSizeParam, ZSizeParam, false);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating random reactions in whole cell space has taken time: ","Execution in threads")));

        CellEngineExecutionTimeStatisticsObject.PrintMeasureTime();
    }
    CATCH("generating random reactions for whole cell space full atom")
}

bool CellEngineFullAtomSimulationSpace::MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false)
        return MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBoundsForThreadsDiffusion(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, CurrentThreadPos);
    else
        return MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBoundsForMPIProcessesDiffusion(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, NeighborProcessesIndexes, VectorOfParticlesToSendToNeighborProcesses, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsReactions(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtomsReactions(ParticleKindObjectForProduct.ListOfAtoms, ParticlesInSector, ParticleKindObjectForProduct.Radius, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam);
}

bool CellEngineFullAtomSimulationSpace::CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsDiffusion(const ParticleKind& ParticleKindObjectForProduct, ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
{
    return CheckFreeSpaceAndBoundsForListOfAtomsDiffusion(ParticleKindObjectForProduct.ListOfAtoms, ParticlesInSector, ParticleKindObjectForProduct.Radius, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam);
}

bool CellEngineFullAtomSimulationSpace::CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMPIMultiProcessDiffusion(const ParticleSenderStructMPIMultiProcess& ParticleSenderToInsert, UniqueIdUInt64Type& UniqueIdOfNewCreatedParticle)
{
    try
    {
        const auto SimulationSpaceSectorBoundsObject = SimulationSpaceSectorBounds().SetParametersForChosenSector(ParticleSenderToInsert.ReceiverSectorPos.X, ParticleSenderToInsert.ReceiverSectorPos.Y, ParticleSenderToInsert.ReceiverSectorPos.Z, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.SizeOfParticlesSectorX, CellEngineConfigDataObject.SizeOfParticlesSectorY, CellEngineConfigDataObject.SizeOfParticlesSectorZ);
        const auto& ParticleKindToCheck = ParticlesKindsManagerObject.GetParticleKind(ParticleSenderToInsert.ParticleKindId);

        if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsDiffusion(ParticleKindToCheck, Particles, ParticleSenderToInsert.NewParticleCenter.X, ParticleSenderToInsert.NewParticleCenter.Y, ParticleSenderToInsert.NewParticleCenter.Z, SimulationSpaceSectorBoundsObject) == true)
        {
            UniqueIdOfNewCreatedParticle = AddNewParticle(Particle(GetNewFreeParticleUniqueId(), ParticleSenderToInsert.ParticleKindId, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
            FillParticleElementsInSpace(UniqueIdOfNewCreatedParticle, ParticleKindToCheck, ParticleSenderToInsert.NewParticleCenter.X, ParticleSenderToInsert.NewParticleCenter.Y, ParticleSenderToInsert.NewParticleCenter.Z);

            return true;
        }
    }
    CATCH("")

    return false;
}

bool CellEngineFullAtomSimulationSpace::CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(const ParticleSenderStructMultiThreaded& ParticleSenderToInsert, const Particle& ParticleObject)
{
    try
    {
        const auto SimulationSpaceSectorBoundsObject = SimulationSpaceSectorBounds().SetParametersForChosenSector(ParticleSenderToInsert.ReceiverSectorPos.X, ParticleSenderToInsert.ReceiverSectorPos.Y, ParticleSenderToInsert.ReceiverSectorPos.Z, CellEngineConfigDataObject.ShiftCenterX, CellEngineConfigDataObject.ShiftCenterY, CellEngineConfigDataObject.ShiftCenterZ, CellEngineConfigDataObject.SizeOfParticlesSectorX, CellEngineConfigDataObject.SizeOfParticlesSectorY, CellEngineConfigDataObject.SizeOfParticlesSectorZ);
        const auto& ParticleKindToCheck = ParticlesKindsManagerObject.GetParticleKind(ParticleSenderToInsert.ParticleKindId);

        if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsDiffusion(ParticleKindToCheck, Particles, ParticleSenderToInsert.NewParticleCenter.X, ParticleSenderToInsert.NewParticleCenter.Y, ParticleSenderToInsert.NewParticleCenter.Z, SimulationSpaceSectorBoundsObject) == true)
        {
            #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
            auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
            LoggersManagerObject.LogUnconditional(STREAM("Particle in thread = (" << CurrentThreadPos.ThreadPosX << "," << CurrentThreadPos.ThreadPosY << "," << CurrentThreadPos.ThreadPosZ << ") CurrentThreadIndex = " << CurrentThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << ParticleSenderToInsert.NewParticleCenter.X << ", " << ParticleSenderToInsert.NewParticleCenter.Y << ", " << ParticleSenderToInsert.NewParticleCenter.Z << ") PARTCLE_INDEX = " << ParticleObject.ParticleUniqueId));
            #endif

            GetParticles().insert(std::move(make_pair(ParticleObject.ParticleUniqueId, std::move(ParticleObject))));

            return true;
        }
    }
    CATCH("")

    return false;
}

void CellEngineFullAtomSimulationSpace::WriteNumberOfParticlesToFile(const bool PrintEverySectorSize) const
{
    try
    {
        UnsignedInt TotalNumberOfParticles = 0;
        UnsignedInt NumberOfNotEmptySectors = 0;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
        {
            TotalNumberOfParticles += Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size();

            if (Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size() > 0)
                NumberOfNotEmptySectors++;

            if (PrintEverySectorSize == true)
                LoggersManagerObject.Log(STREAM("SECTOR = (" << ParticleSectorXIndex << "," << ParticleSectorYIndex  << "," << ParticleSectorZIndex << ") SIZE = " << Particles[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.size()));
        }

        LoggersManagerObject.Log(STREAM("NUMBER OF PARTICLES = " << TotalNumberOfParticles << " NUMBER OF NOT EMPTY SECTORS = " << NumberOfNotEmptySectors));
    }
    CATCH("writing number of particles to file")
}

void CellEngineFullAtomSimulationSpace::WriteNumberOfParticlesInEachSectorToFile() const
{
    try
    {
        UnsignedInt SectorsHistogram[1000];
        for (auto& SectorHistogram : SectorsHistogram)
            SectorHistogram = 0;

        FOR_EACH_SECTOR_IN_XYZ_ONLY
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

        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds | views::values)
            if (ParticleKindObject.ListOfAtoms.empty() == true)
            {
                LoggersManagerObject.Log(STREAM("Particle Kind Without Atoms = " << ParticleKindObject.ParticleKindId << " " << ParticleKindObject.IdStr));
                CountParticleKindsWithoutAtoms++;
            }

        LoggersManagerObject.Log(STREAM("NUMBER OF PARTICLE KINDS WITHOUT VOXELS = " << CountParticleKindsWithoutAtoms << " " << ParticlesKindsManagerObject.ParticlesKinds.size()));
    }
    CATCH("getting number of particles kinds without atoms")
}