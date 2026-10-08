
#include <set>
#include <map>
#include <algorithm>

#include "Combinatorics.h"
#include "CellEngineTypes.h"
#include "CellEngineUseful.h"

#include "CellEngineDataFile.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineExecutionTimeStatistics.h"
#include "CellEngineMPIProcess.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

constexpr bool PrintDetailsOfMakingReaction = true;

using namespace std;

template <class T>
static T sqr(T A)
{
    return A * A;
}

template <class T>
static inline void UpdateNeighbourPointsForChosenElement(T UpdateFunction)
{
    try
    {
        for (SignedInt XPos = 0; XPos <= 2; XPos++)
            for (SignedInt YPos = 0; YPos <= 2; YPos++)
                for (SignedInt ZPos = 0; ZPos <= 2; ZPos++)
                    UpdateFunction(XPos, YPos, ZPos);
    }
    CATCH("updating probability of move from electric interaction for selected particle")
}

void CellEngineSimulationSpace::UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(const Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], const double MultiplyElectricChargeFactor)
{
    try
    {
        UpdateNeighbourPointsForChosenElement([&NeighbourPoints](SignedInt X, SignedInt Y, SignedInt Z){ (*NeighbourPoints)[X][Y][Z] = 0; });

        for (const auto& NeighbourParticleIndexObjectToWrite : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
        {
            if (const Particle& NeighbourParticleObject = GetParticleFromIndex(NeighbourParticleIndexObjectToWrite); NeighbourParticleObject.ElectricCharge != 0)
            {
                for (SignedInt X = 0; X <= 2; X++)
                    for (SignedInt Y = 0; Y <= 2; Y++)
                        for (SignedInt Z = 0; Z <= 2; Z++)
                            if (X != 1 && Y != 1 && Z != 1)
                                if ((NeighbourParticleObject.Center.X < ParticleObject.Center.X && ParticleObject.Center.X + (X - 1) < ParticleObject.Center.X) ||
                                    (NeighbourParticleObject.Center.Y < ParticleObject.Center.Y && ParticleObject.Center.Y + (Y - 1) < ParticleObject.Center.Y) ||
                                    (NeighbourParticleObject.Center.Z < ParticleObject.Center.Z && ParticleObject.Center.Z + (Z - 1) < ParticleObject.Center.Z) ||
                                    (NeighbourParticleObject.Center.X > ParticleObject.Center.X && ParticleObject.Center.X + (X - 1) > ParticleObject.Center.X) ||
                                    (NeighbourParticleObject.Center.Y > ParticleObject.Center.Y && ParticleObject.Center.Y + (Y - 1) > ParticleObject.Center.Y) ||
                                    (NeighbourParticleObject.Center.Z > ParticleObject.Center.Z && ParticleObject.Center.Z + (Z - 1) > ParticleObject.Center.Z)
                                )
                                {
                                    (*NeighbourPoints)[X][Y][Z] += static_cast<ElectricChargeType>((-1.0 * NeighbourParticleObject.ElectricCharge * ParticleObject.ElectricCharge) * MultiplyElectricChargeFactor / sqr(DistanceOfParticles(ParticleObject, NeighbourParticleObject)));
                                    (*NeighbourPoints)[X][Y][Z] = (*NeighbourPoints)[X][Y][Z] < 0 ? 0 : (*NeighbourPoints)[X][Y][Z];

                                    DEBUGLOG(LoggersManagerObject.Log(STREAM("new value after change from neighbour = " << to_string((*NeighbourPoints)[X][Y][Z]) << " " << to_string(static_cast<ElectricChargeType>(X - 1)) << " "<< to_string(static_cast<ElectricChargeType>(Y - 1)) << " " << to_string(static_cast<ElectricChargeType>(Z - 1))));)
                                }

                DEBUGLOG(LoggersManagerObject.Log(STREAM("ParticleIndex of neighbour particle = " << to_string(NeighbourParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(NeighbourParticleIndexObjectToWrite).EntityId) << " Electric Charge = " << to_string(NeighbourParticleObject.ElectricCharge) << " Electric Charge = " << to_string(ParticleObject.ElectricCharge) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(NeighbourParticleIndexObjectToWrite).EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(NeighbourParticleObject.ChainId) : '0') << " GENOME INDEX = " << NeighbourParticleObject.GenomeIndex));)
            }
        }
    }
    CATCH("updating probability of move from electric interaction for selected particle")
}

void CellEngineSimulationSpace::GenerateOneStepOfElectricDiffusionForOneParticle(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdUnsignedInt ParticleIndexParam, ElectricChargeType (*NeighbourPoints)[3][3][3], const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        if (GetParticleFromIndex(ParticleIndexParam).ElectricCharge != 0)
        {
            Particle& ParticleObject = GetParticleFromIndex(ParticleIndexParam);

            DEBUGLOG(LoggersManagerObject.Log(STREAM("EntityId = " << to_string(ParticleObject.EntityId) << " ElectricCharge = " << to_string(ParticleObject.ElectricCharge)));)

            const auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);

            switch (TypeOfLookingForParticles)
            {
                case TypesOfLookingForParticlesInProximity::FromChosenParticleAsCenter : FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, ParticleObject.Center.X - ParticleKindObject.XSizeDiv2 - AdditionalSpaceBoundFactor, ParticleObject.Center.Y - ParticleKindObject.YSizeDiv2 - AdditionalSpaceBoundFactor, ParticleObject.Center.Z - ParticleKindObject.ZSizeDiv2 - AdditionalSpaceBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalSpaceBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalSpaceBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalSpaceBoundFactor); break;
                case TypesOfLookingForParticlesInProximity::InChosenSectorOfSimulationSpace : FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam); break;
                default: break;
            }

            UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(ParticleObject, NeighbourPoints, MultiplyElectricChargeFactor);

            vector<vector3<SignedInt>> MoveVectors;
            UpdateNeighbourPointsForChosenElement([&MoveVectors](const SignedInt X, const SignedInt Y, const SignedInt Z){ MoveVectors.emplace_back(X - 1, Y - 1, Z - 1); });

            vector<int> DiscreteDistribution;
            DiscreteDistribution.reserve(9);

            UpdateNeighbourPointsForChosenElement([&NeighbourPoints, &DiscreteDistribution](const SignedInt X, const SignedInt Y, const SignedInt Z){ DiscreteDistribution.emplace_back((*NeighbourPoints)[X][Y][Z]); });

            #ifdef SIMULATION_DETAILED_DEBUG_LOG
            UnsignedInt NumberOfElement = 0;
            UpdateNeighbourPointsForChosenElement([&NeighbourPoints, &MoveVectors, &NumberOfElement](const SignedInt X, const SignedInt Y, const SignedInt Z){ LoggersManagerObject.Log(STREAM("Element[" << NumberOfElement << "] = " << to_string((*NeighbourPoints)[X][Y][Z]) + " for (X,Y,Z) = (" << to_string(MoveVectors[NumberOfElement].X) << "," << to_string(MoveVectors[NumberOfElement].Y) << "," << to_string(MoveVectors[NumberOfElement].Z) << ")")); NumberOfElement++; });
            #endif

            discrete_distribution<int> UniformDiscreteDistributionMoveParticleDirectionObject(DiscreteDistribution.begin(), DiscreteDistribution.end());

            const UnsignedInt RandomMoveVectorIndex = UniformDiscreteDistributionMoveParticleDirectionObject(mt64R);
            auto EmptyParticlesIter = GetParticles().end();

            if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false && CellEngineConfigDataObject.NonParallelProcessesExecution == false)
                MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleObject, Particles, EmptyParticlesIter, CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ListOfParticlesToChangeSectors, CurrentSectorPos, MoveVectors[RandomMoveVectorIndex].X, MoveVectors[RandomMoveVectorIndex].Y, MoveVectors[RandomMoveVectorIndex].Z, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
            else
                MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleObject, Particles, EmptyParticlesIter, ListOfParticlesToChangeSectors, CurrentSectorPos, MoveVectors[RandomMoveVectorIndex].X, MoveVectors[RandomMoveVectorIndex].Y, MoveVectors[RandomMoveVectorIndex].Z, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

            DEBUGLOG(LoggersManagerObject.Log(STREAM("Random Index = " << to_string(RandomMoveVectorIndex) << " " << to_string(MoveVectors[RandomMoveVectorIndex].X) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Y) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Z) << endl));)
        }
    }
    CATCH("generating one step of electric diffusion for one particle")
}

tuple<IndexesChosenForReactionType, bool> CellEngineSimulationSpace::ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject)
{
    const auto start_time1 = chrono::high_resolution_clock::now();

    bool AllAreZero = false;

    IndexesChosenForReactionType ParticlesIndexesChosenForReaction;

    vector<UnsignedInt> ReactantsCounters(ReactionObject.Reactants.size());

    try
    {
        IndexesChosenForReactionType AllParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction;

        for (UnsignedInt ReactantIndex = 0; ReactantIndex < ReactionObject.Reactants.size(); ReactantIndex++)
            ReactantsCounters[ReactantIndex] = ReactionObject.Reactants[ReactantIndex].Counter;

        for (const auto& ParticleObjectIndex : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
        {
            auto& ParticleObjectTestedForReaction = GetParticleFromIndex(ParticleObjectIndex);

            DEBUGLOG(LoggersManagerObject.Log(STREAM("ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z)));)

            vector<ParticleKindForChemicalReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = ranges::find_if(ReactionObject.Reactants, [&ParticleObjectTestedForReaction](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return ParticleKindForReactionObjectParam.EntityId == ParticleObjectTestedForReaction.EntityId && CompareFitnessOfParticle(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
            else
                ReactantIterator = ranges::find_if(ReactionObject.Reactants, [&ParticleObjectTestedForReaction, this](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsDNA(ParticleKindForReactionObjectParam.EntityId) == true && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByVectorLoop, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

            auto PositionInReactants = distance(ReactionObject.Reactants.cbegin(), ReactantIterator);

            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == true)
                if (ReactantIterator != ReactionObject.Reactants.cend() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == false)
                    NucleotidesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.cend() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == true)
                ParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.cend() && ReactantsCounters[PositionInReactants] > 0)
            {
                AllParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);
                DEBUGLOG(LoggersManagerObject.Log(STREAM("CHOSEN ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z) << endl));)
                ReactantsCounters[PositionInReactants]--;
            }

            AllAreZero = ranges::all_of(std::as_const(ReactantsCounters), [](const UnsignedInt& Counter){ return Counter == 0; });
            if (AllAreZero == true)
            {
                DEBUGLOG(LoggersManagerObject.Log(STREAM("ALL ARE ZERO"));)
                break;
            }

            DEBUGLOG(LoggersManagerObject.Log(STREAM(""));)
        }

        const auto start_time2 = chrono::high_resolution_clock::now();

        if (ReactionObject.SpecialReactionFunction != nullptr)
            ReactionObject.SpecialReactionFunction(this, AllParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);

        const auto stop_time2 = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions += chrono::duration(stop_time2 - start_time2);
    }
    CATCH("choosing particles for reaction from all particles in proximity")

    const auto stop_time1 = chrono::high_resolution_clock::now();

    CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions += chrono::duration(stop_time1 - start_time1);

    if (AllAreZero == true)
    {
        DEBUGLOG(LoggersManagerObject.Log(STREAM("ALL ARE ZERO AT END = " << to_string(ParticlesIndexesChosenForReaction.size()))));
        return { ParticlesIndexesChosenForReaction, true };
    }
    else
        return { IndexesChosenForReactionType(), false };
}

#ifdef SIMULATION_DETAILED_DEBUG_LOG
static void LogParticleData(const UniqueIdUnsignedInt ParticleIndex, const UnsignedInt CenterIndex, const ListOfAtomsType& Centers, const ParticleKind& ParticleKindObjectForProduct, const vector3_16& ParticleKindElement)
{
    LoggersManagerObject.Log(STREAM(endl));
    LoggersManagerObject.Log(STREAM("I " << ParticleIndex << " " << Centers.size() << " " << CenterIndex << endl));
    LoggersManagerObject.Log(STREAM("C " << Centers.size() << " " << CenterIndex << " " << Centers[CenterIndex].X << " " << Centers[CenterIndex].Y << " " << Centers[CenterIndex].Z << endl));
    LoggersManagerObject.Log(STREAM("P " << ParticleKindObjectForProduct.XSizeDiv2 << " " << ParticleKindObjectForProduct.YSizeDiv2 << " " << ParticleKindObjectForProduct.ZSizeDiv2 << endl));
    LoggersManagerObject.Log(STREAM("K " << ParticleKindElement.X << " " << ParticleKindElement.Y << " " << ParticleKindElement.Z << endl));
}
#endif

bool CellEngineSimulationSpace::CancelChemicalReaction(const vector<UniqueIdUnsignedInt>& CreatedParticlesIndexes, const ListOfCentersType& Centers, const vector<Particle>& ParticlesBackup, const chrono::high_resolution_clock::time_point start_time, const ParticleKind& ParticleKindObjectForProduct, const char PlaceStr)
{
    try
    {
        DEBUGLOG(LoggersManagerObject.Log(STREAM("CANCELLED REACTION IN BOUNDS " << PlaceStr << " = " << ActualSimulationSpaceSectorBoundsObject.StartXPos << " " << ActualSimulationSpaceSectorBoundsObject.StartYPos << " "  << ActualSimulationSpaceSectorBoundsObject.StartZPos << " " << ActualSimulationSpaceSectorBoundsObject.EndXPos << " " << ActualSimulationSpaceSectorBoundsObject.EndYPos << " " << ActualSimulationSpaceSectorBoundsObject.EndZPos << " " << ParticleKindObjectForProduct.ListOfVoxels.size() << " " << ParticleKindObjectForProduct.ListOfAtoms.size() << " " << ParticleKindObjectForProduct.EntityId));)

        for (const auto& CreatedParticleIndex : CreatedParticlesIndexes)
        {
            RemoveParticle(CreatedParticleIndex, true, false);

            DEBUGLOG(LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("Current Thread = " << CurrentThreadPos.ThreadPosX << " " << CurrentThreadPos.ThreadPosY << " "  << CurrentThreadPos.ThreadPosZ));)

            if (CellEngineConfigDataObject.GatherCancelledParticlesIndexes == true)
            {
                if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false && CellEngineConfigDataObject.NonParallelProcessesExecution == false)
                    CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->CancelledParticlesIndexes.insert(pair(CreatedParticleIndex, CreatedParticleIndex));
                else
                    CancelledParticlesIndexes.insert(pair(CreatedParticleIndex, CreatedParticleIndex));
            }
        }

        NumberOfCancelledReactions++;
        RemovedParticlesInReactions -= ParticlesBackup.size();
        RestoredParticlesInCancelledReactions += ParticlesBackup.size();

        PlaceStr == 'A' ? NumberOfCancelledAReactions++ : NumberOfCancelledBReactions++;

        for (auto& Particle : ParticlesBackup)
            AddNewParticle(std::move(Particle));

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingCancelledChemicalReactions += chrono::duration(stop_time - start_time);
    }
    CATCH("cancelling chemical reaction")

    return false;
}

bool CellEngineSimulationSpace::PlaceNewProductParticleInSpaceDeterminedFromPositionOfFormerReactantParticleOrCancelReaction(const UniqueIdUnsignedInt ParticleIndex, const vector<Particle>& ParticlesBackup, const vector<UniqueIdUnsignedInt>& CreatedParticlesIndexes, const UnsignedInt CenterIndex, const ListOfCentersType& Centers, ParticleKind& ParticleKindObjectForProduct, const chrono::high_resolution_clock::time_point start_time)
{
    try
    {
        const vector3_Real32 NewCenter(Centers[CenterIndex].X - ParticleKindObjectForProduct.XSizeDiv2, Centers[CenterIndex].Y - ParticleKindObjectForProduct.YSizeDiv2, Centers[CenterIndex].Z - ParticleKindObjectForProduct.ZSizeDiv2);

        if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsReactions(ParticleKindObjectForProduct, Particles, CurrentSectorPos, NewCenter.X, NewCenter.Y, NewCenter.Z, GetBoundsForThreadSector()) == true)
        {
            FillParticleElementsInSpace(ParticleIndex, ParticleKindObjectForProduct, NewCenter.X, NewCenter.Y, NewCenter.Z);
            AddedParticlesInReactions++;
        }
        else
            return CancelChemicalReaction(CreatedParticlesIndexes, Centers, ParticlesBackup, start_time, ParticleKindObjectForProduct, 'A');
    }
    CATCH("placing product particle in space in determined position or cancel reaction")

    return true;
}

bool CellEngineSimulationSpace::PlaceNewProductParticleInSpaceInNewRandomPositionOrCancelReaction(const UniqueIdUnsignedInt ParticleIndex, const vector<Particle>& ParticlesBackup, const vector<UniqueIdUnsignedInt>& CreatedParticlesIndexes, const UnsignedInt CenterIndex, const ListOfCentersType& Centers, ParticleKind& ParticleKindObjectForProduct, const chrono::high_resolution_clock::time_point start_time)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirectionX_int64t(ActualSimulationSpaceSectorBoundsObject.StartXPos, ActualSimulationSpaceSectorBoundsObject.EndXPos);
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirectionY_int64t(ActualSimulationSpaceSectorBoundsObject.StartYPos, ActualSimulationSpaceSectorBoundsObject.EndYPos);
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirectionZ_int64t(ActualSimulationSpaceSectorBoundsObject.StartZPos, ActualSimulationSpaceSectorBoundsObject.EndZPos);

        const auto SimulationSpaceSectorBoundsObject = GetBoundsForThreadSector();

        bool FoundFreePlace = false;

        UnsignedInt NumberOfTries = 0;
        while (NumberOfTries < 1000)
        {
            NumberOfTries++;

            const auto RandomVectorX = GetRandomValue<uniform_int_distribution, SignedInt>(UniformDistributionObjectMoveParticleDirectionX_int64t);
            const auto RandomVectorY = GetRandomValue<uniform_int_distribution, SignedInt>(UniformDistributionObjectMoveParticleDirectionY_int64t);
            const auto RandomVectorZ = GetRandomValue<uniform_int_distribution, SignedInt>(UniformDistributionObjectMoveParticleDirectionZ_int64t);

            DEBUGLOG(LoggersManagerObject.Log(STREAM("R1 = (" << RandomVectorX << " " << RandomVectorY << " " << RandomVectorZ << ") (" << SimulationSpaceSectorBoundsObject.StartXPos << "," << SimulationSpaceSectorBoundsObject.EndXPos << ") (" << SimulationSpaceSectorBoundsObject.StartYPos << "," << SimulationSpaceSectorBoundsObject.EndYPos << ") (" << SimulationSpaceSectorBoundsObject.StartZPos << "," << SimulationSpaceSectorBoundsObject.EndZPos << ")"));)

            if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElementsReactions(ParticleKindObjectForProduct, Particles, CurrentSectorPos, RandomVectorX, RandomVectorY, RandomVectorZ, SimulationSpaceSectorBoundsObject) == true)
            {
                #ifdef SIMULATION_DETAILED_DEBUG_LOG
                LoggersManagerObject.Log(STREAM("R2 = (" << RandomVectorX << " " << RandomVectorY << " " << RandomVectorZ << ") (" << SimulationSpaceSectorBoundsObject.StartXPos << "," << SimulationSpaceSectorBoundsObject.EndXPos << ") (" << SimulationSpaceSectorBoundsObject.StartYPos << "," << SimulationSpaceSectorBoundsObject.EndYPos << ") (" << SimulationSpaceSectorBoundsObject.StartZPos << "," << SimulationSpaceSectorBoundsObject.EndZPos << ")"));

                const auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(RandomVectorX, RandomVectorY, RandomVectorZ);
                if (SectorPosX != CurrentSectorPos.SectorPosX || SectorPosY != CurrentSectorPos.SectorPosY || SectorPosZ != CurrentSectorPos.SectorPosZ)
                    LoggersManagerObject.Log(STREAM("R4" << " Error particle not in proper sector " << " SectorPosX = " << SectorPosX << " SectorPosY = " << SectorPosY << " SectorPosZ = " << SectorPosZ << " SectorPosXS = " << CurrentSectorPos.SectorPosX << " SectorPosYS = " << CurrentSectorPos.SectorPosY << " SectorPosZS = " << CurrentSectorPos.SectorPosZ));
                #endif

                FoundFreePlace = true;

                FillParticleElementsInSpace(ParticleIndex, ParticleKindObjectForProduct, RandomVectorX, RandomVectorY, RandomVectorZ);

                #ifdef SIMULATION_DETAILED_DEBUG_LOG
                auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::CheckCenterForSector<Particle>(GetParticleFromIndex(ParticleIndex), Particles, "CK1");
                if (SectorPosX1 != CurrentSectorPos.SectorPosX || SectorPosY1 != CurrentSectorPos.SectorPosY || SectorPosZ1 != CurrentSectorPos.SectorPosZ)
                    LoggersManagerObject.Log(STREAM("R5" << " Error particle not in proper sector " << " SectorPosX = " << SectorPosX1 << " SectorPosY = " << SectorPosY1 << " SectorPosZ = " << SectorPosZ1 << " SectorPosXS = " << CurrentSectorPos.SectorPosX << " SectorPosYS = " << CurrentSectorPos.SectorPosY << " SectorPosZS = " << CurrentSectorPos.SectorPosZ));

                const auto [SectorPosX11, SectorPosY11, SectorPosZ11] = CellEngineUseful::GetSectorPos(GetParticleFromIndex(ParticleIndex).Center.X, GetParticleFromIndex(ParticleIndex).Center.Y, GetParticleFromIndex(ParticleIndex).Center.Z);

                LoggersManagerObject.Log(STREAM("R2 = (" << RandomVectorX << " " << RandomVectorY << " " << RandomVectorZ << ") (" << SimulationSpaceSectorBoundsObject.StartXPos << "," << SimulationSpaceSectorBoundsObject.EndXPos << ") (" << SimulationSpaceSectorBoundsObject.StartYPos << "," << SimulationSpaceSectorBoundsObject.EndYPos << ") (" << SimulationSpaceSectorBoundsObject.StartZPos << "," << SimulationSpaceSectorBoundsObject.EndZPos << ")" << " ListOfAtoms.size() = " << ParticleKindObjectForProduct.ListOfAtoms.size()));
                LoggersManagerObject.Log(STREAM("R2A = " << ParticleIndex << " " << GetParticleFromIndex(ParticleIndex).Center.X << " " << GetParticleFromIndex(ParticleIndex).Center.Y << " " << GetParticleFromIndex(ParticleIndex).Center.Z << " " << RandomVectorX << " " << RandomVectorY << " " << RandomVectorZ << " " << SimulationSpaceSectorBoundsObject.StartXPos << " " << SimulationSpaceSectorBoundsObject.EndXPos << " " << SimulationSpaceSectorBoundsObject.StartYPos << " " << SimulationSpaceSectorBoundsObject.EndYPos << " " << SimulationSpaceSectorBoundsObject.StartZPos << " " << SimulationSpaceSectorBoundsObject.EndZPos << " ListOfAtoms.size() = " << ParticleKindObjectForProduct.ListOfAtoms.size()));
                LoggersManagerObject.Log(STREAM("R2B = " << ParticleIndex << " " << SectorPosX << " " << SectorPosY << " " << SectorPosZ << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX11 << " " << SectorPosY11 << " " << SectorPosZ11));
                #endif

                break;
            }
        }
        if (FoundFreePlace == false)
            return CancelChemicalReaction(CreatedParticlesIndexes, Centers, ParticlesBackup, start_time, ParticleKindObjectForProduct, 'B');
        else
            AddedParticlesInReactions++;
    }
    CATCH("placing product particle in space in random position or cancel reaction")

    return true;
}

bool CellEngineSimulationSpace::MakeChemicalReaction(ChemicalReaction& ReactionObject)
{
    const auto start_time = chrono::high_resolution_clock::now();

    try
    {
        auto [ParticlesIndexesChosenForReaction, FoundInProximity] = ChooseParticlesForReactionFromAllParticlesInProximity(ReactionObject);

        if (FoundInProximity == false)
            return false;

        DEBUGLOG(LoggersManagerObject.Log(STREAM("Reaction Step 1 - chosen particles for reaction from all particles in proximity" << endl));)

        ListOfCentersType Centers;
        vector<Particle> ParticlesBackup;
        for (const auto& ParticleIndexChosenForReaction : ParticlesIndexesChosenForReaction | views::keys)
            EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(ParticleIndexChosenForReaction, Centers, ParticlesBackup);

        DEBUGLOG(LoggersManagerObject.Log(STREAM("Reaction Step 2 - erasing particles chosen for reaction" << endl));)

        DEBUGLOG(LoggersManagerObject.Log(STREAM("Centers size = " << to_string(Centers.size()) << endl));)

        vector<UniqueIdUnsignedInt> CreatedParticlesIndexes;

        UnsignedInt CenterIndex = 0;
        for (const auto& ReactionProduct : ReactionObject.Products)
        {
            UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ReactionProduct.EntityId, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));

            CreatedParticlesIndexes.emplace_back(ParticleIndex);

            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

            if (CellEngineConfigDataObject.TypeOfReactionsDeterminedByFindingNewPositionForPlacingNewParticles == CellEngineConfigData::TypesOfReactionsDeterminedByFindingNewPositionForPlacingNewParticles::NewParticlePlacedInNewRandomPosition)
            {
                if (PlaceNewProductParticleInSpaceInNewRandomPositionOrCancelReaction(ParticleIndex, ParticlesBackup, CreatedParticlesIndexes, CenterIndex, Centers, ParticleKindObjectForProduct, start_time) == false)
                    return false;
            }
            else
            {
                if (PlaceNewProductParticleInSpaceDeterminedFromPositionOfFormerReactantParticleOrCancelReaction(ParticleIndex, ParticlesBackup, CreatedParticlesIndexes, CenterIndex, Centers, ParticleKindObjectForProduct, start_time) == false)
                    return false;
            }

            CenterIndex++;
        }

        NumberOfExecutedReactions++;

        DEBUGLOG(LoggersManagerObject.Log(STREAM("Reaction Step 3 - Reaction finished" << endl));)

        if (SaveReactionsStatisticsBool == true)
            SaveReactionForStatistics(ReactionObject);
    }
    CATCH("making chemical reaction")

    const auto stop_time = chrono::high_resolution_clock::now();

    CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactions += chrono::duration(stop_time - start_time);

    return true;
};

bool CellEngineSimulationSpace::IsChemicalReactionPossible(const ChemicalReaction& ReactionObject)
{
    return ranges::all_of(ReactionObject.Reactants, [this](const ParticleKindForChemicalReaction& ReactionReactant){ return ReactionReactant.Counter <= LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ReactionReactant.EntityId]; });
}

void CellEngineSimulationSpace::PrepareRandomReaction()
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMainRandomCondition_Uint64t(0, 1);
        if (GetRandomValue<uniform_int_distribution, UnsignedInt>(UniformDistributionObjectMainRandomCondition_Uint64t) == 0)
            return;
    }
    CATCH("preparing random reaction")
}

set<UnsignedInt> CellEngineSimulationSpace::GetAllPossibleReactionsFromParticlesInProximity()
{
    set<UnsignedInt> PossibleReactionsIdNums;

    try
    {
        for (const auto& ParticleKindFoundInProximityObject : LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity | views::keys)
            for (const auto& ReactionIdNum : ParticlesKindsManagerObject.GetParticleKind(ParticleKindFoundInProximityObject).AssociatedChemicalReactions)
                if (auto ReactionIter = ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.find(ReactionIdNum); ReactionIter != ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.end())
                    if (IsChemicalReactionPossible(ChemicalReactionsManagerObject.ChemicalReactions[ReactionIter->second]) == true)
                        PossibleReactionsIdNums.insert(ReactionIdNum);
    }
    CATCH("finding and executing random reaction v3")

    return PossibleReactionsIdNums;
}

void CellEngineSimulationSpace::FindAndExecuteRandomReactionVersion3(const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        if (set<UnsignedInt> PossibleReactionsIdNums = GetAllPossibleReactionsFromParticlesInProximity(); PossibleReactionsIdNums.empty() == false)
        {
            string ListOfPossibleReactions;
            for (const auto& PossibleReactionsIdNum : PossibleReactionsIdNums)
                ListOfPossibleReactions += to_string(PossibleReactionsIdNum) + ",";

            DEBUGLOG(LoggersManagerObject.Log(STREAM("ListOfPossibleReactions = " << ListOfPossibleReactions));)

            std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, PossibleReactionsIdNums.size() - 1);
            const auto ReactionIdNum = *std::next(std::begin(PossibleReactionsIdNums), static_cast<int>(GetRandomValue<uniform_int_distribution, UnsignedInt>(UniformDistributionObjectUint64t)));

            DEBUGLOG(LoggersManagerObject.Log(STREAM("Random ReactionIdNum = " << ReactionIdNum));)

            FindAndExecuteChosenReaction(ReactionIdNum);
        }
        else
        {
            DEBUGLOG(LoggersManagerObject.Log(STREAM("NONE REACTION FOUND for particles kinds in proximity "));)
        }
    }
    CATCH("finding and executing random reaction v3")
}

void CellEngineSimulationSpace::FindAndExecuteRandomReaction(const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        FindAndExecuteRandomReactionVersion3(MaxNumberOfReactants);
    }
    CATCH("finding and executing random reaction")
}

bool CellEngineSimulationSpace::FindAndExecuteChosenReaction(const UnsignedInt ReactionId)
{
    try
    {
        if (const auto ReactionIter = ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.find(ReactionId); ReactionIter != ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.end())
        {
            auto& ReactionObject = ChemicalReactionsManagerObject.ChemicalReactions[ReactionIter->second];

            if (const bool IsPossible = IsChemicalReactionPossible(ReactionObject); IsPossible == true)
            {
                DEBUGLOG(LoggersManagerObject.Log(STREAM("CHOSEN REACTION POSSIBLE" << endl));)

                if (MakeChemicalReaction(ReactionObject) == false)
                {
                    DEBUGLOG(LoggersManagerObject.Log(STREAM("Chosen reaction not executed!"));)
                    return false;
                }
                else
                    return true;
            }
            else
            {
                DEBUGLOG(LoggersManagerObject.Log(STREAM("Chosen reaction impossible!"));)
            }
        }
        else
        {
            DEBUGLOG(LoggersManagerObject.Log(STREAM("Chosen reaction Id not found!"));)
        }
    }
    CATCH("finding and executing chosen reaction")

    return false;
}

void CellEngineSimulationSpace::SaveHistogramOfParticlesStatisticsToFile() const
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("SORTED HISTOGRAM OF PARTICLES"));

        for (const auto& [EntityId, Counter1, Counter2, Difference] : ParticlesKindsHistogramComparisons.back())
            if (EntityId != 0)
                LoggersManagerObject.LogStatistics(STREAM("PARTICLE KIND = " << EntityId << " NAME = " << ParticlesKindsManagerObject.GetParticleKind(EntityId).IdStr << " DIFFERENCE = " << Difference << " " << Counter1 << " " << Counter2 << endl));
    }
    CATCH("saving histograms of particles statistics to file")
}

void CellEngineSimulationSpace::SaveNumberOfParticlesStatisticsToFile()
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("THREAD = " << CurrentThreadIndex << " X = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInX << " Y = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInY << " Z = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInZ));

        GetNumberOfParticlesFromParticleKind(ParticlesKindsManagerObject.GetParticleKindFromStrId("M_glc__D_e")->EntityId);

        int CounterOfParticles = 0;
        FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
            if (ParticleObject.first != 0)
                if (auto IPK = ParticlesKindsManagerObject.ParticlesKinds.find(GetParticleFromIndex(ParticleObject.first).EntityId); IPK != ParticlesKindsManagerObject.ParticlesKinds.end() && IPK->second.IdStr == "M_glc__D_e")
                    CounterOfParticles++;

        LoggersManagerObject.LogStatistics(STREAM("Particle Name = " << "D-Glucose" << " Number of Particles = " << CounterOfParticles));
    }
    CATCH("saving number of particles statistics to file")
}

#ifdef SHORTER_CODE
void CellEngineSimulationSpace::SaveReactionsStatisticsToFile() const
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("NUMBER OF REACTIONS = " << SavedReactionsMap[SimulationStepNumber - 1].size() << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
        for (const auto& [ReactionDataFirst, ReactionDataSecond] : SavedReactionsMap[SimulationStepNumber - 1])
        {
            LoggersManagerObject.LogStatistics(STREAM("REACTION ID = " << ReactionDataFirst << " REACTION NAME = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionDataFirst).ReactionName << " REACTION ID_STR = #" << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionDataFirst).ReactionIdStr << "# REACTION COUNTER = " << ReactionDataSecond.Counter));
            LoggersManagerObject.LogStatistics(STREAM("REACTANTS_STR = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionDataFirst).ReactantsStr));
            LoggersManagerObject.LogStatistics(STREAM("PRODUCTS = " << ChemicalReactionsManager::GetStringOfSortedParticlesDataNames(ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionDataFirst).Products) << endl));
        }
    }
    CATCH("saving reactions statistics to file")
}
#else
void CellEngineSimulationSpace::SaveReactionsStatisticsToFileExtended() const
{
    try
    {
        for (const auto& ReactionData : SavedReactionsMap[SimulationStepNumber - 1])
        {
            LoggersManagerObject.LogStatistics(STREAM("REACTION ID = " << ReactionData.second.ReactionId << " REACTION NAME = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactionName << " REACTION ID_STR = #" << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactionIdStr << "# REACTION COUNTER = " << ReactionData.second.Counter));
            LoggersManagerObject.LogStatistics(STREAM("REACTANTS_STR = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactantsStr << endl));
        }
    }
    CATCH("saving reactions statistics to file extended")
}
#endif

void CellEngineSimulationSpace::SetMakeSimulationStepNumberZero()
{
    MakeSimulationStepNumberZeroForStatistics();
    IncSimulationStepNumberForStatistics();
    GenerateNewEmptyElementsForContainersForStatistics();
}

void CellEngineSimulationSpace::SetIncSimulationStepNumber()
{
    IncSimulationStepNumberForStatistics();
    GenerateNewEmptyElementsForContainersForStatistics();
}

void CellEngineSimulationSpace::SaveParticlesStatisticsOnce()
{
    SaveParticlesStatistics();
}

#ifdef WELL_STIRRED
std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticlesVersion3() const
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size() - 1);
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticles(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    return GetRandomParticlesVersion3();
}
#endif
