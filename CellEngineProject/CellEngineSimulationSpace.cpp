
#include <set>
#include <map>
#include <thread>
#include <barrier>
#include <algorithm>

#include "FileUtils.h"
#include "DateTimeUtils.h"
#include "Combinatorics.h"
#include "DoublyLinkedList.h"

#include "CellEngineTypes.h"
#include "CellEngineUseful.h"
#include "CellEngineConstants.h"

#include "CellEngineDataFile.h"
#include "CellEngineSimulationSpace.h"

#include <condition_variable>

#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineExecutionTimeStatistics.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

constexpr bool PrintDetailsOfMakingReaction = true;

using namespace std;

template <class T>
T sqr(T A)
{
    return A * A;
}

template <class T>
inline void UpdateNeighbourPointsForChosenElement(T UpdateFunction)
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
            Particle& NeighbourParticleObject = GetParticleFromIndex(NeighbourParticleIndexObjectToWrite);
            if (NeighbourParticleObject.ElectricCharge != 0)
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
                                    #ifdef SIMULATION_DETAILED_LOG
                                    LoggersManagerObject.Log(STREAM("new value after change from neighbour = " << to_string((*NeighbourPoints)[X][Y][Z]) << " " << to_string(static_cast<ElectricChargeType>(X - 1)) << " "<< to_string(static_cast<ElectricChargeType>(Y - 1)) << " " << to_string(static_cast<ElectricChargeType>(Z - 1))));
                                    #endif
                                }

                #ifdef SIMULATION_DETAILED_LOG
                LoggersManagerObject.Log(STREAM("ParticleIndex of neighbour particle = " << to_string(NeighbourParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(NeighbourParticleIndexObjectToWrite).EntityId) << " Electric Charge = " << to_string(NeighbourParticleObject.ElectricCharge) << " Electric Charge = " << to_string(ParticleObject.ElectricCharge) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(NeighbourParticleIndexObjectToWrite).EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(NeighbourParticleObject.ChainId) : '0') << " GENOME INDEX = " << NeighbourParticleObject.GenomeIndex));
                #endif
            }
        }
    }
    CATCH("updating probability of move from electric interaction for selected particle")
}

void CellEngineSimulationSpace::GenerateOneStepOfElectricDiffusionForOneParticle(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt ParticleIndexParam, ElectricChargeType (*NeighbourPoints)[3][3][3], const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        if (GetParticleFromIndex(ParticleIndexParam).ElectricCharge != 0)
        {
            Particle& ParticleObject = GetParticleFromIndex(ParticleIndexParam);

            #ifdef SIMULATION_DETAILED_LOG
            LoggersManagerObject.Log(STREAM("EntityId = " << to_string(ParticleObject.EntityId) << " ElectricCharge = " << to_string(ParticleObject.ElectricCharge)));
            #endif

            auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);

            switch (TypeOfLookingForParticles)
            {
                case TypesOfLookingForParticlesInProximity::FromChosenParticleAsCenter : FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, ParticleObject.Center.X - ParticleKindObject.XSizeDiv2 - AdditionalSpaceBoundFactor, ParticleObject.Center.Y - ParticleKindObject.YSizeDiv2 - AdditionalSpaceBoundFactor, ParticleObject.Center.Z - ParticleKindObject.ZSizeDiv2 - AdditionalSpaceBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalSpaceBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalSpaceBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalSpaceBoundFactor); break;
                case TypesOfLookingForParticlesInProximity::InChosenSectorOfSimulationSpace : FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam); break;
                default: break;
            }

            UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(ParticleObject, NeighbourPoints, MultiplyElectricChargeFactor);

            vector<vector3<SignedInt>> MoveVectors;
            UpdateNeighbourPointsForChosenElement([&MoveVectors](SignedInt X, SignedInt Y, SignedInt Z){ MoveVectors.emplace_back(X - 1, Y - 1, Z - 1); });

            vector<int> DiscreteDistribution;
            DiscreteDistribution.reserve(9);

            UpdateNeighbourPointsForChosenElement([&NeighbourPoints, &DiscreteDistribution](SignedInt X, SignedInt Y, SignedInt Z){ DiscreteDistribution.emplace_back((*NeighbourPoints)[X][Y][Z]); });
            #ifdef SIMULATION_DETAILED_LOG
            UnsignedInt NumberOfElement = 0;
            UpdateNeighbourPointsForChosenElement([&NeighbourPoints, &MoveVectors, &NumberOfElement](SignedInt X, SignedInt Y, SignedInt Z){ LoggersManagerObject.Log(STREAM("Element[" << NumberOfElement << "] = " << to_string((*NeighbourPoints)[X][Y][Z]) + " for (X,Y,Z) = (" << to_string(MoveVectors[NumberOfElement].X) << "," << to_string(MoveVectors[NumberOfElement].Y) << "," << to_string(MoveVectors[NumberOfElement].Z) << ")")); NumberOfElement++; });
            #endif

            discrete_distribution<int> UniformDiscreteDistributionMoveParticleDirectionObject(DiscreteDistribution.begin(), DiscreteDistribution.end());

            UnsignedInt RandomMoveVectorIndex = UniformDiscreteDistributionMoveParticleDirectionObject(mt64R);
            MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleObject, Particles, CurrentSectorPos, MoveVectors[RandomMoveVectorIndex].X, MoveVectors[RandomMoveVectorIndex].Y, MoveVectors[RandomMoveVectorIndex].Z, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

            #ifdef SIMULATION_DETAILED_LOG
            LoggersManagerObject.Log(STREAM("Random Index = " << to_string(RandomMoveVectorIndex) << " " << to_string(MoveVectors[RandomMoveVectorIndex].X) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Y) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Z) << endl));
            #endif
        }
    }
    CATCH("generating one step of electric diffusion for one particle")
}

tuple<vector<pair<UniqueIdInt, UnsignedInt>>, bool> CellEngineSimulationSpace::ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject)
{
    const auto start_time1 = chrono::high_resolution_clock::now();

    bool AllAreZero = false;

    vector<pair<UniqueIdInt, UnsignedInt>> NucleotidesIndexesChosenForReaction, ParticlesIndexesChosenForReaction, AllParticlesIndexesChosenForReaction;

    vector<UnsignedInt> ReactantsCounters(ReactionObject.Reactants.size());

    try
    {
        for (UnsignedInt ReactantIndex = 0; ReactantIndex < ReactionObject.Reactants.size(); ReactantIndex++)
            ReactantsCounters[ReactantIndex] = ReactionObject.Reactants[ReactantIndex].Counter;

        for (const auto& ParticleObjectIndex : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
        {
            auto& ParticleObjectTestedForReaction = GetParticleFromIndex(ParticleObjectIndex);

            LoggersManagerObject.Log(STREAM("ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z)));

            vector<ParticleKindForChemicalReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return ParticleKindForReactionObjectParam.EntityId == ParticleObjectTestedForReaction.EntityId && CompareFitnessOfParticle(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
            else
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction, this](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsDNA(ParticleKindForReactionObjectParam.EntityId) == true && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByVectorLoop, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

            auto PositionInReactants = distance(ReactionObject.Reactants.cbegin(), ReactantIterator);

            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == true)
                if (ReactantIterator != ReactionObject.Reactants.cend() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == false)
                    NucleotidesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.cend() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == true)
                ParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.cend() && ReactantsCounters[PositionInReactants] > 0)
            {
                AllParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);
                LoggersManagerObject.Log(STREAM("CHOSEN ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z) << endl));
                ReactantsCounters[PositionInReactants]--;
            }

            AllAreZero = all_of(ReactantsCounters.cbegin(), ReactantsCounters.cend(), [this](const UnsignedInt& Counter){ return Counter == 0; });
            if (AllAreZero == true)
            {
                LoggersManagerObject.Log(STREAM("ALL ARE ZERO"));
                break;
            }
            LoggersManagerObject.Log(STREAM(""));
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
        LoggersManagerObject.Log(STREAM("ALL ARE ZERO AT END = " << to_string(ParticlesIndexesChosenForReaction.size())));
        return { ParticlesIndexesChosenForReaction, true };
    }
    else
        return { vector<pair<UniqueIdInt, UnsignedInt>>(), false };
}

void LogParticleData(const UniqueIdInt ParticleIndex, const UnsignedInt CenterIndex, const ListOfAtomsType& Centers, const ParticleKind& ParticleKindObjectForProduct, const vector3_16& ParticleKindElement)
{
    LoggersManagerObject.Log(STREAM(endl));
    LoggersManagerObject.Log(STREAM("I " << ParticleIndex << " " << Centers.size() << " " << CenterIndex << endl));
    LoggersManagerObject.Log(STREAM("C " << Centers.size() << " " << CenterIndex << " " << Centers[CenterIndex].X << " " << Centers[CenterIndex].Y << " " << Centers[CenterIndex].Z << endl));
    LoggersManagerObject.Log(STREAM("P " << ParticleKindObjectForProduct.XSizeDiv2 << " " << ParticleKindObjectForProduct.YSizeDiv2 << " " << ParticleKindObjectForProduct.ZSizeDiv2 << endl));
    LoggersManagerObject.Log(STREAM("K " << ParticleKindElement.X << " " << ParticleKindElement.Y << " " << ParticleKindElement.Z << endl));
}

bool CellEngineSimulationSpace::CancelChemicalReaction(const vector<UniqueIdInt>& CreatedParticlesIndexes, const ListOfCentersType& Centers, const vector<Particle>& ParticlesBackup, const chrono::high_resolution_clock::time_point start_time, const ParticleKind& ParticleKindObjectForProduct, const char PlaceStr)
{
    try
    {
        LoggersManagerObject.Log(STREAM("CANCELLED PARTICLE IN BOUNDS " << PlaceStr << " = " << ActualSimulationSpaceSectorBoundsObject.StartXPos << " " << ActualSimulationSpaceSectorBoundsObject.StartYPos << " "  << ActualSimulationSpaceSectorBoundsObject.StartZPos << " " << ActualSimulationSpaceSectorBoundsObject.EndXPos << " " << ActualSimulationSpaceSectorBoundsObject.EndYPos << " " << ActualSimulationSpaceSectorBoundsObject.EndZPos << " " << ParticleKindObjectForProduct.ListOfVoxels.size() << " " << ParticleKindObjectForProduct.ListOfAtoms.size() << " " << ParticleKindObjectForProduct.EntityId));

        for (const auto& CreatedParticleIndex : CreatedParticlesIndexes)
        {
            RemoveParticle(CreatedParticleIndex, true);

            CancelledParticlesIndexes.insert(pair(CreatedParticleIndex, CreatedParticleIndex));
        }

        RemovedParticlesInReactions -= ParticlesBackup.size();

        for (auto& Particle : ParticlesBackup)
            AddNewParticle(move(Particle));

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingCancelledChemicalReactions += chrono::duration(stop_time - start_time);
    }
    CATCH("cancelling chemical reaction")

    return false;
}

bool CellEngineSimulationSpace::PlaceProductParticleInSpaceInDeterminedPositionOrCancelReaction(const UniqueIdInt ParticleIndex, const vector<Particle>& ParticlesBackup, const vector<UniqueIdInt>& CreatedParticlesIndexes, const UnsignedInt CenterIndex, const ListOfCentersType& Centers, ParticleKind& ParticleKindObjectForProduct, const chrono::high_resolution_clock::time_point start_time)
{
    try
    {
        const vector3_Real32 NewCenter(Centers[CenterIndex].X - ParticleKindObjectForProduct.XSizeDiv2, Centers[CenterIndex].Y - ParticleKindObjectForProduct.YSizeDiv2, Centers[CenterIndex].Z - ParticleKindObjectForProduct.ZSizeDiv2);

        if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(ParticleKindObjectForProduct, Particles, CurrentSectorPos, NewCenter.X, NewCenter.Y, NewCenter.Z, GetBoundsForThreadSector()) == true)
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

bool CellEngineSimulationSpace::PlaceProductParticleInSpaceInRandomPositionOrCancelReaction(const UniqueIdInt ParticleIndex, const vector<Particle>& ParticlesBackup, const vector<UniqueIdInt>& CreatedParticlesIndexes, const UnsignedInt CenterIndex, const ListOfCentersType& Centers, ParticleKind& ParticleKindObjectForProduct, const chrono::high_resolution_clock::time_point start_time)
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

            auto RandomVectorX = UniformDistributionObjectMoveParticleDirectionX_int64t(mt64R);
            auto RandomVectorY = UniformDistributionObjectMoveParticleDirectionY_int64t(mt64R);
            auto RandomVectorZ = UniformDistributionObjectMoveParticleDirectionZ_int64t(mt64R);

            LoggersManagerObject.Log(STREAM("R1 = " << RandomVectorX << " " << RandomVectorY << " " << RandomVectorZ << " " << SimulationSpaceSectorBoundsObject.StartXPos << " " << SimulationSpaceSectorBoundsObject.EndXPos << " " << SimulationSpaceSectorBoundsObject.StartYPos << " " << SimulationSpaceSectorBoundsObject.EndYPos << " " << SimulationSpaceSectorBoundsObject.StartZPos << " " << SimulationSpaceSectorBoundsObject.EndZPos));

            if (CheckIfSpaceIsEmptyAndIsInBoundsForParticleElements(ParticleKindObjectForProduct, Particles, CurrentSectorPos, RandomVectorX, RandomVectorY, RandomVectorZ, SimulationSpaceSectorBoundsObject) == true)
            {
                FoundFreePlace = true;

                FillParticleElementsInSpace(ParticleIndex, ParticleKindObjectForProduct, RandomVectorX, RandomVectorY, RandomVectorZ);

                LoggersManagerObject.Log(STREAM("R2 = " << RandomVectorX << " " << RandomVectorY << " " << RandomVectorZ << " " << SimulationSpaceSectorBoundsObject.StartXPos << " " << SimulationSpaceSectorBoundsObject.EndXPos << " " << SimulationSpaceSectorBoundsObject.StartYPos << " " << SimulationSpaceSectorBoundsObject.EndYPos << " " << SimulationSpaceSectorBoundsObject.StartZPos << " " << SimulationSpaceSectorBoundsObject.EndZPos << " ListOfAtoms.size() = " << ParticleKindObjectForProduct.ListOfAtoms.size()));

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

        LoggersManagerObject.Log(STREAM("Reaction Step 1 - chosen particles for reaction from all particles in proximity" << endl));

        ListOfCentersType Centers;
        vector<Particle> ParticlesBackup;
        for (const auto& ParticleIndexChosenForReaction : ParticlesIndexesChosenForReaction)
            EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(ParticleIndexChosenForReaction.first, Centers, ParticlesBackup);

        LoggersManagerObject.Log(STREAM("Reaction Step 2 - erasing particles chosen for reaction" << endl));

        LoggersManagerObject.Log(STREAM("Centers size = " << to_string(Centers.size()) << endl));

        vector<UniqueIdInt> CreatedParticlesIndexes;

        UnsignedInt CenterIndex = 0;
        for (const auto& ReactionProduct : ReactionObject.Products)
        {
            UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ReactionProduct.EntityId, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));

            CreatedParticlesIndexes.emplace_back(ParticleIndex);

            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

            bool ResultOfPlacingParticle;
            if (CenterIndex < Centers.size())
                ResultOfPlacingParticle = PlaceProductParticleInSpaceInDeterminedPositionOrCancelReaction(ParticleIndex, ParticlesBackup, CreatedParticlesIndexes, CenterIndex, Centers, ParticleKindObjectForProduct, start_time);
            else
                ResultOfPlacingParticle = PlaceProductParticleInSpaceInRandomPositionOrCancelReaction(ParticleIndex, ParticlesBackup, CreatedParticlesIndexes, CenterIndex, Centers, ParticleKindObjectForProduct, start_time);

            if (ResultOfPlacingParticle == false)
                return false;

            CenterIndex++;
        }

        LoggersManagerObject.Log(STREAM("Reaction Step 3 - Reaction finished" << endl));

        if (SaveReactionsStatisticsBool == true)
            SaveReactionForStatistics(ReactionObject);
    }
    CATCH("making chemical reaction")

    const auto stop_time = chrono::high_resolution_clock::now();

    CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactions += chrono::duration(stop_time - start_time);

    return true;
};

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticlesVersion3(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants) const
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size() - 1);
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticlesVersion2(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        auto BitsValuesString = Combinations::CreateBoolStringFromInt64BitState(GenerateCombinationsStateNumber);
        LoggersManagerObject.Log(STREAM("GenerateCombinationsStateNumber NEXT = " << BitsValuesString));

        for (UnsignedInt ReactantNumberBitValuePos = 0; ReactantNumberBitValuePos < MaxNumberOfReactants; ReactantNumberBitValuePos++)
            if (BitsValuesString[ReactantNumberBitValuePos] == '1')
            {
                RandomParticlesTypes.emplace_back(std::next(std::begin(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity), static_cast<int>(ReactantNumberBitValuePos))->first);

                LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumberBitValuePos) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
            }

        GenerateCombinationsStateNumber = Combinations::NextNumberWithTheSameNumberOf1Bits(GenerateCombinationsStateNumber);
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticlesVersion1(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size() - 1);

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(std::next(std::begin(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity), static_cast<int>(UniformDistributionObjectUint64t(mt64R)))->first);

            LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
        }
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticles(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    return GetRandomParticlesVersion3(NumberOfReactants, MaxNumberOfReactants);
}

bool CellEngineSimulationSpace::IsChemicalReactionPossible(const ChemicalReaction& ReactionObject)
{
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForChemicalReaction& ReactionReactant){ return ReactionReactant.Counter <= LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ReactionReactant.EntityId]; });
}

void CellEngineSimulationSpace::PrepareRandomReaction()
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMainRandomCondition_Uint64t(0, 1);
        if (UniformDistributionObjectMainRandomCondition_Uint64t(mt64R) == 0)
            return;
    }
    CATCH("preparing random reaction")
}

void CellEngineSimulationSpace::FindAndExecuteRandomReactionVersion4(const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        for (const auto& ParticleKindFoundInProximityObject : LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity)
            for (const auto& ReactionIdNum : ParticlesKindsManagerObject.GetParticleKind(ParticleKindFoundInProximityObject.first).AssociatedChemicalReactions)
                if (FindAndExecuteChosenReaction(ReactionIdNum) == true)
                    goto EndLoops;

        EndLoops:;
    }
    CATCH("finding and executing random reaction v4")
}

set<UnsignedInt> CellEngineSimulationSpace::GetAllPossibleReactionsFromParticlesInProximity()
{
    set<UnsignedInt> PossibleReactionsIdNums;

    try
    {
        for (const auto& ParticleKindFoundInProximityObject : LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity)
            for (const auto& ReactionIdNum : ParticlesKindsManagerObject.GetParticleKind(ParticleKindFoundInProximityObject.first).AssociatedChemicalReactions)
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
        set<UnsignedInt> PossibleReactionsIdNums = GetAllPossibleReactionsFromParticlesInProximity();

        if (PossibleReactionsIdNums.empty() == false)
        {
            string ListOfPossibleReactions;
            for (const auto& PossibleReactionsIdNum : PossibleReactionsIdNums)
                ListOfPossibleReactions += to_string(PossibleReactionsIdNum) + ",";
            LoggersManagerObject.Log(STREAM("ListOfPossibleReactions = " << ListOfPossibleReactions));

            std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, PossibleReactionsIdNums.size() - 1);
            auto ReactionIdNum = *std::next(std::begin(PossibleReactionsIdNums), static_cast<int>(UniformDistributionObjectUint64t(mt64R)));
            LoggersManagerObject.Log(STREAM("Random ReactionIdNum = " << ReactionIdNum));
            FindAndExecuteChosenReaction(ReactionIdNum);
        }
        else
            LoggersManagerObject.Log(STREAM("NONE REACTION FOUND for particles kinds in proximity "));
    }
    CATCH("finding and executing random reaction v3")
}

void CellEngineSimulationSpace::FindAndExecuteRandomReactionVersion2(const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        LoggersManagerObject.Log(STREAM("MAX NUMBER OF REACTANTS = " << MaxNumberOfReactants));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectNumberOfReactants_Uint64t(1, MaxNumberOfReactants);

        UnsignedInt NumberOfRandom = 0;
        while (NumberOfRandom < 100)
        {
            NumberOfRandom++;

            UnsignedInt NumberOfReactants = UniformDistributionObjectNumberOfReactants_Uint64t(mt64R);

            LoggersManagerObject.Log(STREAM("NumberOfReactants = " << NumberOfReactants));
            LoggersManagerObject.Log(STREAM("MaxNumberOfReactants = " << MaxNumberOfReactants));
            UnsignedInt NumberOfCombinations = Combinations::NumberOfCombinations(MaxNumberOfReactants, NumberOfReactants);
            LoggersManagerObject.Log(STREAM("NumberOfCombinations = " << NumberOfCombinations));
            GenerateCombinationsStateNumber = Combinations::SetKBitsInNumber(MaxNumberOfReactants, NumberOfReactants);
            LoggersManagerObject.Log(STREAM("GenerateCombinationsStateNumber START = " << Combinations::CreateBoolStringFromInt64BitState(GenerateCombinationsStateNumber)));

            UnsignedInt NumberOfTries = 0;
            UnsignedInt NumberOfAllPossibleTries = min(UnsignedInt(100), NumberOfCombinations);
            while (NumberOfTries < NumberOfAllPossibleTries)
            {
                NumberOfTries++;
                LoggersManagerObject.Log(STREAM("Number Of Tries = " << NumberOfTries));

                if (TryToMakeRandomChemicalReaction(NumberOfReactants, MaxNumberOfReactants) == true)
                    goto BreakLoop;
            }
        }
        BreakLoop:;
    }
    CATCH("finding and executing random reaction v2")
}

void CellEngineSimulationSpace::FindAndExecuteRandomReactionVersion1(const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        LoggersManagerObject.Log(STREAM("MAX NUMBER OF REACTANTS = " << MaxNumberOfReactants));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectNumberOfReactants_Uint64t(1, MaxNumberOfReactants);

        UnsignedInt NumberOfTries = 0;
        while (NumberOfTries <= 100)
        {
            NumberOfTries++;

            UnsignedInt NumberOfReactants = UniformDistributionObjectNumberOfReactants_Uint64t(mt64R);

            if (TryToMakeRandomChemicalReaction(NumberOfReactants, MaxNumberOfReactants) == true)
                break;
        }
    }
    CATCH("finding and executing random reaction v1")
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
        if (auto ReactionIter = ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.find(ReactionId); ReactionIter != ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.end())
        {
            auto& ReactionObject = ChemicalReactionsManagerObject.ChemicalReactions[ReactionIter->second];

            bool IsPossible = IsChemicalReactionPossible(ReactionObject);
            if (IsPossible == true)
            {
                LoggersManagerObject.Log(STREAM("CHOSEN REACTION POSSIBLE" << endl));

                if (MakeChemicalReaction(ReactionObject) == false)
                {
                    LoggersManagerObject.Log(STREAM("Chosen reaction not executed!"));
                    return false;
                }
                else
                    return true;
            }
            else
                LoggersManagerObject.Log(STREAM("Chosen reaction impossible!"));
        }
        else
            LoggersManagerObject.Log(STREAM("Chosen reaction Id not found!"));
    }
    CATCH("finding and executing chosen reaction")

    return false;
}

void CellEngineSimulationSpace::SaveHistogramOfParticlesStatisticsToFile()
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("SORTED HISTOGRAM OF PARTICLES"));

        for (const auto& ParticleKindHistogramComparisonObject : ParticlesKindsHistogramComparisons.back())
            if (ParticleKindHistogramComparisonObject.EntityId != 0)
                LoggersManagerObject.LogStatistics(STREAM("PARTICLE KIND = " << ParticleKindHistogramComparisonObject.EntityId << " NAME = " << ParticlesKindsManagerObject.GetParticleKind(ParticleKindHistogramComparisonObject.EntityId).IdStr << " DIFFERENCE = " << ParticleKindHistogramComparisonObject.Difference << " " << ParticleKindHistogramComparisonObject.Counter1 << " " << ParticleKindHistogramComparisonObject.Counter2 << endl));
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
        FOR_EACH_PARTICLE_IN_XYZ_CONST
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
        for (const auto& ReactionData : SavedReactionsMap[SimulationStepNumber - 1])
        {
            LoggersManagerObject.LogStatistics(STREAM("REACTION ID = " << ReactionData.first << " REACTION NAME = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.first).ReactionName << " REACTION ID_STR = #" << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.first).ReactionIdStr << "# REACTION COUNTER = " << ReactionData.second.Counter));
            LoggersManagerObject.LogStatistics(STREAM("REACTANTS_STR = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.first).ReactantsStr));
            LoggersManagerObject.LogStatistics(STREAM("PRODUCTS = " << ChemicalReactionsManager::GetStringOfSortedParticlesDataNames(ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.first).Products) << endl));
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
