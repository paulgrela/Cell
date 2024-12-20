
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

#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

constexpr bool PrintDetailsOfMakingReaction = false;

using namespace std;

inline void GetRangeOfParticlesForRandomParticles(UniqueIdInt& StartParticleIndexParam, UniqueIdInt& EndParticleIndexParam, UniqueIdInt MaxParticleIndex)
{
    if (EndParticleIndexParam == 0)
    {
        StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
        EndParticleIndexParam = MaxParticleIndex;
    }
}

void CellEngineSimulationSpace::GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-1, 1);

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(GetParticleFromIndex(ParticleIndex), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("generating one step of diffusion for selected range of particles")
}

void CellEngineSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpace(const bool InBounds, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-1, 1);

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
        for (auto& ParticleInProximityIndex : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
            if (CellEngineUseful::IsDNA(GetParticleFromIndex(ParticleInProximityIndex).EntityId) == false)
                if (InBounds == false)
                    MoveParticleByVectorIfSpaceIsEmpty(GetParticleFromIndex(ParticleInProximityIndex), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R));
                else
                    MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(GetParticleFromIndex(ParticleInProximityIndex), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("generating one step of diffusion for selected space")
}

void CellEngineSimulationSpace::GenerateNStepsOfDiffusionForBigPartOfCellSpace(const bool InBounds, const UnsignedInt SizeNMultiplyFactor, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam - SizeNMultiplyFactor * XStepParam; PosX <= XStartParam + SizeNMultiplyFactor * XStepParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam - SizeNMultiplyFactor * YStepParam; PosY <= YStartParam + SizeNMultiplyFactor * YStepParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam - SizeNMultiplyFactor * ZStepParam; PosZ <= ZStartParam + SizeNMultiplyFactor * ZStepParam; PosZ += ZStepParam)
                        GenerateOneStepOfDiffusionForSelectedSpace(InBounds, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating diffusion for big part of cell")
}

void CellEngineSimulationSpace::GenerateNStepsOfDiffusionForWholeCellSpace(const bool InBounds, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                        GenerateOneStepOfDiffusionForSelectedSpace(InBounds, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating diffusion for whole cell space")
}

template <class T>
T sqr(T A)
{
    return A * A;
}

template <class T>
inline void UpdateNeighbourPointsForChosenVoxel(T UpdateFunction)
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

void CellEngineSimulationSpace::UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], const double MultiplyElectricChargeFactor)
{
    try
    {
        UpdateNeighbourPointsForChosenVoxel([&NeighbourPoints](SignedInt X, SignedInt Y, SignedInt Z){ (*NeighbourPoints)[X][Y][Z] = 0; });

        for (const auto& NeighbourParticleIndexObjectToWrite : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
        {
            Particle& NeighbourParticleObject = GetParticleFromIndex(NeighbourParticleIndexObjectToWrite);
            if (NeighbourParticleObject.ElectricCharge != 0)
            {
                for (UnsignedInt X = 0; X <= 2; X++)
                    for (UnsignedInt Y = 0; Y <= 2; Y++)
                        for (UnsignedInt Z = 0; Z <= 2; Z++)
                            if (X != 1 && Y != 1 && Z != 1)
                                if ((NeighbourParticleObject.Center.X < ParticleObject.Center.X && ParticleObject.Center.X +(X - 1) < ParticleObject.Center.X) ||
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
                case TypesOfLookingForParticlesInProximity::InChosenVoxelSpace : FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam); break;
                default: break;
            }

            UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(ParticleObject, NeighbourPoints, MultiplyElectricChargeFactor);

            vector<vector3<SignedInt>> MoveVectors;
            UpdateNeighbourPointsForChosenVoxel([&MoveVectors](SignedInt X, SignedInt Y, SignedInt Z){ MoveVectors.emplace_back(X - 1, Y - 1, Z - 1); });

            vector<int> DiscreteDistribution;
            DiscreteDistribution.reserve(9);

            UpdateNeighbourPointsForChosenVoxel([&NeighbourPoints, &DiscreteDistribution](SignedInt X, SignedInt Y, SignedInt Z){ DiscreteDistribution.emplace_back((*NeighbourPoints)[X][Y][Z]); });
            #ifdef SIMULATION_DETAILED_LOG
            UnsignedInt NumberOfElement = 0;
            UpdateNeighbourPointsForChosenVoxel([&NeighbourPoints, &MoveVectors, &NumberOfElement](SignedInt X, SignedInt Y, SignedInt Z){ LoggersManagerObject.Log(STREAM("Element[" << NumberOfElement << "] = " << to_string((*NeighbourPoints)[X][Y][Z]) + " for (X,Y,Z) = (" << to_string(MoveVectors[NumberOfElement].X) << "," << to_string(MoveVectors[NumberOfElement].Y) << "," << to_string(MoveVectors[NumberOfElement].Z) << ")")); NumberOfElement++; });
            #endif

            discrete_distribution<int> UniformDiscreteDistributionMoveParticleDirectionObject(DiscreteDistribution.begin(), DiscreteDistribution.end());

            UnsignedInt RandomMoveVectorIndex = UniformDiscreteDistributionMoveParticleDirectionObject(mt64R);
            MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleObject, MoveVectors[RandomMoveVectorIndex].X, MoveVectors[RandomMoveVectorIndex].Y, MoveVectors[RandomMoveVectorIndex].Z, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

            #ifdef SIMULATION_DETAILED_LOG
            LoggersManagerObject.Log(STREAM("Random Index = " << to_string(RandomMoveVectorIndex) << " " << to_string(MoveVectors[RandomMoveVectorIndex].X) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Y) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Z) << endl));
            #endif
        }
    }
    CATCH("generating one step of electric diffusion for one particle")
}


void CellEngineSimulationSpace::GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        ElectricChargeType NeighbourPoints[3][3][3];

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            GenerateOneStepOfElectricDiffusionForOneParticle(TypeOfLookingForParticles, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, ParticleIndex, &NeighbourPoints, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating one step of electric diffusion")
}

void CellEngineSimulationSpace::GenerateOneStepOfElectricDiffusionForSelectedSpace(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        ElectricChargeType NeighbourPoints[3][3][3];

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        auto ParticlesSortedByCapacityFoundInProximityCopy(LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity);

        for (auto& ParticleInProximityIndex : ParticlesSortedByCapacityFoundInProximityCopy)
            if (CellEngineUseful::IsDNA(GetParticleFromIndex(ParticleInProximityIndex).EntityId) == false)
                GenerateOneStepOfElectricDiffusionForOneParticle(TypeOfLookingForParticles, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, ParticleInProximityIndex, &NeighbourPoints, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating one step of electric diffusion")
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

        if (AllAreZero == true || (AllAreZero == false && (ReactionObject.ReactionIdNum == 30 || ReactionObject.ReactionIdNum == 80 || ReactionObject.ReactionIdNum == 70)))
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

void LogParticleData(const UniqueIdInt ParticleIndex, const UnsignedInt CenterIndex, const vector<vector3_16>& Centers, const ParticleKind& ParticleKindObjectForProduct, const vector3_16& ParticleKindVoxel)
{
    LoggersManagerObject.Log(STREAM(endl));
    LoggersManagerObject.Log(STREAM("I " << ParticleIndex << " " << Centers.size() << " " << CenterIndex << endl));
    LoggersManagerObject.Log(STREAM("C " << Centers.size() << " " << CenterIndex << " " << Centers[CenterIndex].X << " " << Centers[CenterIndex].Y << " " << Centers[CenterIndex].Z << endl));
    LoggersManagerObject.Log(STREAM("P " << ParticleKindObjectForProduct.XSizeDiv2 << " " << ParticleKindObjectForProduct.YSizeDiv2 << " " << ParticleKindObjectForProduct.ZSizeDiv2 << endl));
    LoggersManagerObject.Log(STREAM("K " << ParticleKindVoxel.X << " " << ParticleKindVoxel.Y << " " << ParticleKindVoxel.Z << endl));

    cout << endl;
    cout << "I " << ParticleIndex << " " << Centers.size() << " " << CenterIndex << endl;
    cout << "C " << Centers[CenterIndex].X << " " << Centers[CenterIndex].Y << " " << Centers[CenterIndex].Z << endl;
    cout << "P " << ParticleKindObjectForProduct.XSizeDiv2 << " " << ParticleKindObjectForProduct.YSizeDiv2 << " " << ParticleKindObjectForProduct.ZSizeDiv2 << endl;
    cout << "K " << ParticleKindVoxel.X << " " << ParticleKindVoxel.Y << " " << ParticleKindVoxel.Z << endl;
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

        vector<vector3_16> Centers;
        for (const auto& ParticleIndexChosenForReaction : ParticlesIndexesChosenForReaction)
            EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(ParticleIndexChosenForReaction.first, Centers);

        LoggersManagerObject.Log(STREAM("Reaction Step 2 - erasing particles chosen for reaction" << endl));

        LoggersManagerObject.Log(STREAM("Centers size = " << to_string(Centers.size()) << endl));

        UnsignedInt CenterIndex = 0;
        for (const auto& ReactionProduct : ReactionObject.Products)
        {
            UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ReactionProduct.EntityId, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));

            if (CurrentThreadIndex == 0 && Particles.contains(ParticleIndex) == false)
                cout << "PARTICLE INDEX ERROR = " << ParticleIndex << endl;
            else
            if (CurrentThreadIndex != 0 && ParticlesForThreads.contains(ParticleIndex) == false)
                cout << "PARTICLE INDEX ERROR IN THREADS = " << ParticleIndex << endl;

            GetParticleFromIndex(ParticleIndex).ListOfVoxels.clear();

            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

            if (CenterIndex < Centers.size())
            {
                bool ErrorInVoxels = false;

                for (const auto& ParticleKindVoxel : ParticleKindObjectForProduct.ListOfVoxels)
                {
                    vector3_64 NewVoxel(Centers[CenterIndex].X - ParticleKindObjectForProduct.XSizeDiv2 + ParticleKindVoxel.X, Centers[CenterIndex].Y - ParticleKindObjectForProduct.YSizeDiv2 + ParticleKindVoxel.Y, Centers[CenterIndex].Z - ParticleKindObjectForProduct.ZSizeDiv2 + ParticleKindVoxel.Z);
                    if (NewVoxel.X > CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension || NewVoxel.Y > CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension || NewVoxel.Z > CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                    {
                        if (PrintDetailsOfMakingReaction == true)
                            LogParticleData(ParticleIndex, CenterIndex, Centers, ParticleKindObjectForProduct, ParticleKindVoxel);
                        ErrorInVoxels = true;
                    }
                    else
                    {
                        if (CurrentThreadIndex == 0 && Particles.contains(ParticleIndex) == false)
                            cout << "PARTICLE INDEX WRONG = " << ParticleIndex << endl;
                        if (CurrentThreadIndex == 0 && Particles.contains(ParticleIndex) == true)
                            FillParticleElementInSpace(ParticleIndex, NewVoxel);
                        if (CurrentThreadIndex != 0 && ParticlesForThreads.contains(ParticleIndex) == false)
                            cout << "PARTICLE INDEX WRONG IN THREADS = " << ParticleIndex << endl;
                        if (CurrentThreadIndex != 0 && ParticlesForThreads.contains(ParticleIndex) == true)
                            FillParticleElementInSpace(ParticleIndex, NewVoxel);
                    }

                    LoggersManagerObject.Log(STREAM("New Centers From Product Added X = " << to_string(NewVoxel.X) << " Y = " << to_string(NewVoxel.Y) << " Z = " << to_string(NewVoxel.Z) << endl));
                }

                {
                    std::lock_guard LockGuardObject{ MainCountingErrorMutexObject };
                    if (ErrorInVoxels == true)
                        ErrorCounter++;
                    AddedParticlesInReactions++;
                }

                GetMinMaxCoordinatesForParticle(GetParticleFromIndex(ParticleIndex), false);
            }
            else
            {
            }

            CenterIndex++;
        }

        LoggersManagerObject.Log(STREAM("Reaction Step 3 - Reaction finished" << endl));

        if (SaveReactionsStatisticsBool == true)
            SaveReactionForStatistics(ReactionObject);
    }
    CATCH("making reaction")

    const auto stop_time = chrono::high_resolution_clock::now();

    CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactions += chrono::duration(stop_time - start_time);

    return true;
};

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticlesVersion3(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
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

void CellEngineSimulationSpace::FindAndExecuteRandomReactionVersion3(const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        set<UnsignedInt> PossibleReactionsIdNums;

        for (const auto& ParticleKindFoundInProximityObject : LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity)
            for (const auto& ReactionIdNum : ParticlesKindsManagerObject.GetParticleKind(ParticleKindFoundInProximityObject.first).AssociatedChemicalReactions)
                if (auto ReactionIter = ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.find(ReactionIdNum); ReactionIter != ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.end())
                    if (IsChemicalReactionPossible(ChemicalReactionsManagerObject.ChemicalReactions[ReactionIter->second]) == true)
                        PossibleReactionsIdNums.insert(ReactionIdNum);

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

void CellEngineSimulationSpace::GenerateOneRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const bool FindParticlesInProximityBool)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityBool == false || (FindParticlesInProximityBool == true && FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true))
            FindAndExecuteRandomReaction(min(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
    }
    CATCH("generating random reaction for selected voxel space")
}

void CellEngineSimulationSpace::GenerateOneRandomReactionForChosenParticle(const Particle& ParticleObject)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(ParticleObject, 20) == true)
            FindAndExecuteRandomReaction(min(LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
    }
    CATCH("generating random reaction for particle")
}

void CellEngineSimulationSpace::GenerateOneChosenReactionForSelectedSpace(const UnsignedInt ReactionId, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        if (FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            if (ReactionId != 0)
                FindAndExecuteChosenReaction(ReactionId);
    }
    CATCH("generating random reaction for particle")
}

void CellEngineSimulationSpace::GenerateOneStepOfRandomReactionsForAllParticles(const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.SelectedForReaction == false)
                GenerateOneRandomReactionForChosenParticle(ParticleObject.second);
    }
    CATCH("generating random reactions for all particles")
}

void CellEngineSimulationSpace::GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            if (auto ParticlesIterator = Particles.find(ParticleIndex); ParticlesIterator != Particles.end())
                GenerateOneRandomReactionForChosenParticle(ParticlesIterator->second);
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineSimulationSpace::GenerateOneStepOfRandomReactionsForOneParticleFromRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt ShiftIndexOfChosenParticle)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);
        GenerateOneRandomReactionForChosenParticle(GetParticleFromIndex(StartParticleIndexParam + ShiftIndexOfChosenParticle));
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineSimulationSpace::SaveNumberOfParticlesStatisticsToFile()
{
    try
    {
        GetNumberOfParticlesFromParticleKind(ParticlesKindsManagerObject.GetParticleKindFromStrId("M_glc__D_e")->EntityId);

        int CounterOfParticles = 0;
        for (const auto& ParticleData : Particles)
            if (ParticleData.first != 0)
            {
                if (ParticlesKindsManagerObject.GetParticleKind(GetParticleFromIndex(ParticleData.first).EntityId).IdStr == "M_glc__D_e")
                    CounterOfParticles++;
            }

        LoggersManagerObject.LogStatistics(STREAM("Particle Name = " << "D-Glucose" << " Number of Particles = " << CounterOfParticles));
    }
    CATCH("saving number of particles statistics to file")
}

void CellEngineSimulationSpace::SaveReactionsStatisticsToFile()
{
    try
    {
        for (const auto& ReactionData : SavedReactionsMap[SimulationStepNumber - 1])
        {
            LoggersManagerObject.LogStatistics(STREAM("REACTION ID = " << ReactionData.second.ReactionId << " REACTION NAME = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactionName << " REACTION ID_STR = #" << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactionIdStr << "# REACTION COUNTER = " << ReactionData.second.Counter));
            LoggersManagerObject.LogStatistics(STREAM("REACTANTS_STR = " << ChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactantsStr << endl));
        }
    }
    CATCH("saving reactions statistics to file")
}

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

void CellEngineSimulationSpace::GenerateNStepsOfOneChosenReactionForWholeCellSpace(const UnsignedInt ReactionId, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                        GenerateOneChosenReactionForSelectedSpace(ReactionId, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for whole cell space")
}

void CellEngineSimulationSpace::GenerateNStepsOfOneRandomReactionForWholeCellSpace(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForFindingParticles = std::chrono::seconds::zero();
        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForSavingFoundParticles = std::chrono::seconds::zero();
        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactions = std::chrono::seconds::zero();
        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions = std::chrono::seconds::zero();
        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions = std::chrono::seconds::zero();

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                        GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam, true);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating random reactions in whole cell space has taken time: ","Execution in threads")));

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForFindingParticles, "Execution of finding particles has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForSavingFoundParticles, "Execution of saving found particles has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactions, "Execution of making chemical reactions has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions, "Execution of making chemical reactions special functions has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions, "Execution of saving choosing particles for making chemical reactions has taken time: ","Execution in threads")));
    }
    CATCH("generating random reactions for whole cell space")
}

void CellEngineSimulationSpace::GenerateNStepsOfOneRandomReactionForBigPartOfCellSpace(const UnsignedInt SizeNMultiplyFactor, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const UnsignedInt NumberOfSimulationSteps)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt Step = 1; Step <= NumberOfSimulationSteps; Step++)
            for (UnsignedInt PosX = XStartParam - SizeNMultiplyFactor * XStepParam; PosX <= XStartParam + SizeNMultiplyFactor * XStepParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam - SizeNMultiplyFactor * YStepParam; PosY <= YStartParam + SizeNMultiplyFactor * YStepParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam - SizeNMultiplyFactor * ZStepParam; PosZ <= ZStartParam + SizeNMultiplyFactor * ZStepParam; PosZ += ZStepParam)
                        GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam, true);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for big part of cell")
}

void LogCenterOfParticleWithThreadIndex(const Particle& ParticleObject, const ThreadIdType ThreadXIndex, const ThreadIdType ThreadYIndex, const ThreadIdType ThreadZIndex)
{
    LoggersManagerObject.Log(STREAM("Center: " << ParticleObject.Center.X << ParticleObject.Center.Y << ParticleObject.Center.Z << endl));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex << endl));
    LoggersManagerObject.Log(STREAM(endl));

    cout << "Center: " << ParticleObject.Center.X << ParticleObject.Center.Y << ParticleObject.Center.Z << endl;
    cout << "THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex << endl;
    cout << endl;
}

void CellEngineSimulationSpace::FirstSendParticlesForThreads(const bool PrintCenterOfParticleWithThreadIndex, const bool PrintTime)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        for (auto& ThreadLocalParticlesInProximityXPos : CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer)
            for (auto& ThreadLocalParticlesInProximityYPos : ThreadLocalParticlesInProximityXPos)
                for (auto& ThreadLocalParticlesInProximityZPos : ThreadLocalParticlesInProximityYPos)
                    ThreadLocalParticlesInProximityZPos->ParticlesForThreads.clear();

        for (const auto& ParticleObject : Particles)
        {
            UnsignedInt ThreadXIndex = floor(ParticleObject.second.Center.X / CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace);
            UnsignedInt ThreadYIndex = floor(ParticleObject.second.Center.Y / CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace);
            UnsignedInt ThreadZIndex = floor(ParticleObject.second.Center.Z / CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace);

            if (PrintCenterOfParticleWithThreadIndex == true)
                LogCenterOfParticleWithThreadIndex(ParticleObject.second, ThreadXIndex, ThreadYIndex, ThreadZIndex);

            CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.insert({ ParticleObject.first, ParticleObject.second });
        }

        const auto stop_time = chrono::high_resolution_clock::now();

        if (PrintTime == true)
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "First sending particles to threads has taken time = ","Execution in threads")));

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    InitiateFreeParticleIndexes(CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads);

                    LoggersManagerObject.Log(STREAM("THREAD[" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << "] SIZE = " << CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.size()));
                }
    }
    CATCH("first sending particles for threads")
}

void PrintThreadIndexes(const ThreadIdType CurrentThreadIndex, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, const UnsignedInt ThreadXIndexNew, const UnsignedInt ThreadYIndexNew, const UnsignedInt ThreadZIndexNew)
{
    LoggersManagerObject.Log(STREAM("Current Thread: " << CurrentThreadIndex << endl));
    LoggersManagerObject.Log(STREAM("Center: " << ThreadXIndexNew << ThreadYIndexNew << ThreadZIndexNew << endl));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex << endl));
    LoggersManagerObject.Log(STREAM(endl));

    cout << "Current Thread: " << CurrentThreadIndex << endl;
    cout << "Center: " << ThreadXIndexNew << ThreadYIndexNew << ThreadZIndexNew << endl;
    cout << "THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex << endl;
    cout << endl;
}

void CellEngineSimulationSpace::ExchangeParticlesBetweenThreads(const UnsignedInt StepOutside) const
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    auto ParticleIter = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.begin();

                    while (ParticleIter != CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.end())
                    {
                        const UnsignedInt NumberOfXVoxelsInOneThreadInVoxelSimulationSpaceDiv2 = (StepOutside % 2 == 0 ? 0 : (CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace / 2));
                        const UnsignedInt NumberOfYVoxelsInOneThreadInVoxelSimulationSpaceDiv2 = (StepOutside % 2 == 0 ? 0 : (CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace / 2));
                        const UnsignedInt NumberOfZVoxelsInOneThreadInVoxelSimulationSpaceDiv2 = (StepOutside % 2 == 0 ? 0 : (CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace / 2));

                        if (StepOutside % 2 != 0)
                            if (ParticleIter->second.Center.X < NumberOfXVoxelsInOneThreadInVoxelSimulationSpaceDiv2 || ParticleIter->second.Center.Y < NumberOfYVoxelsInOneThreadInVoxelSimulationSpaceDiv2 || ParticleIter->second.Center.Z < NumberOfZVoxelsInOneThreadInVoxelSimulationSpaceDiv2)
                            {
                                ++ParticleIter;
                                continue;
                            }

                        UnsignedInt ThreadXIndexNew = floor((ParticleIter->second.Center.X - NumberOfXVoxelsInOneThreadInVoxelSimulationSpaceDiv2) / CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace);
                        UnsignedInt ThreadYIndexNew = floor((ParticleIter->second.Center.Y - NumberOfYVoxelsInOneThreadInVoxelSimulationSpaceDiv2) / CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace);
                        UnsignedInt ThreadZIndexNew = floor((ParticleIter->second.Center.Z - NumberOfZVoxelsInOneThreadInVoxelSimulationSpaceDiv2) / CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace);

                        PrintThreadIndexes(CurrentThreadIndex, ThreadXIndex, ThreadYIndex, ThreadZIndex, ThreadXIndexNew, ThreadYIndexNew, ThreadZIndexNew);

                        if ((ThreadXIndex - 1 != ThreadXIndexNew) || (ThreadYIndex - 1 != ThreadYIndexNew) || (ThreadZIndex - 1 != ThreadZIndexNew))
                            CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew]->ParticlesForThreads.insert(CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.extract(ParticleIter++));
                        else
                            ++ParticleIter;
                    }
                }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Next sending particles to threads v1 has taken time = ","Execution in threads")));
    }
    CATCH("sending particles for threads")
}

void CellEngineSimulationSpace::GatherParticlesForThreadsInMainParticles()
{
    try
    {
        std::lock_guard LockGuardObject{ CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderMenuAndVoxelSimulationSpaceMutexObject };

        const auto start_time = chrono::high_resolution_clock::now();

        Particles.clear();

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                    for (const auto& ParticleObject : CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads)
                        Particles.insert(ParticleObject);

        InitiateFreeParticleIndexes(Particles);

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Next sending particles to threads v2 has taken time = ","Execution in threads")));
    }
    CATCH("sending particles for threads")
}

void CellEngineSimulationSpace::GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    try
    {
        for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
        {
            LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));

            UnsignedInt XStartParam = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace;
            UnsignedInt XEndParam = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace + CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace;
            UnsignedInt YStartParam = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace;
            UnsignedInt YEndParam = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace + CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace;
            UnsignedInt ZStartParam = (ThreadZIndex - 1) * CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace;
            UnsignedInt ZEndParam = (ThreadZIndex - 1) * CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace + CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace;

            if (StepOutside % 2 == 0)
            {
                XStartParam += CellEngineConfigDataObject.NumberOfXVoxelsInOneThreadInVoxelSimulationSpace / 2;
                YStartParam += CellEngineConfigDataObject.NumberOfYVoxelsInOneThreadInVoxelSimulationSpace / 2;
                ZStartParam += CellEngineConfigDataObject.NumberOfZVoxelsInOneThreadInVoxelSimulationSpace / 2;
            }

            for (UnsignedInt PosX = XStartParam; PosX < XEndParam; PosX += CellEngineConfigDataObject.NumberOfXVoxelsInOneSectorInOneThreadInVoxelSimulationSpace)
                for (UnsignedInt PosY = YStartParam; PosY < YEndParam; PosY += CellEngineConfigDataObject.NumberOfYVoxelsInOneSectorInOneThreadInVoxelSimulationSpace)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZEndParam; PosZ += CellEngineConfigDataObject.NumberOfZVoxelsInOneSectorInOneThreadInVoxelSimulationSpace)
                    {
                        LoggersManagerObject.Log(STREAM("XStart = " << XStartParam << " YStart = " << YStartParam << " ZStart = " << ZStartParam << " XEnd = " << XEndParam << " YEnd = " << YEndParam << " ZEnd = " << ZEndParam << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ));

                        if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                            GenerateOneStepOfDiffusionForSelectedSpace(true, PosX, PosY, PosZ, CellEngineConfigDataObject.NumberOfXVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, CellEngineConfigDataObject.NumberOfYVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, CellEngineConfigDataObject.NumberOfZVoxelsInOneSectorInOneThreadInVoxelSimulationSpace);
                        if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                            GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.NumberOfXVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, CellEngineConfigDataObject.NumberOfYVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, CellEngineConfigDataObject.NumberOfZVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, false);
                        if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                            GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.NumberOfXVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, CellEngineConfigDataObject.NumberOfYVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, CellEngineConfigDataObject.NumberOfZVoxelsInOneSectorInOneThreadInVoxelSimulationSpace, true);
                    }
        }
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationSpace::GenerateNStepsOfSimulationForWholeCellSpaceInThreads(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside)
{
    try
    {
        LoggersManagerObject.Log(STREAM("MaxParticleIndex = " << MaxParticleIndex));

        ErrorCounter = 0;
        AddedParticlesInReactions = 0;

        std::barrier SyncPoint(CellEngineConfigDataObject.NumberOfXThreadsInSimulation * CellEngineConfigDataObject.NumberOfYThreadsInSimulation * CellEngineConfigDataObject.NumberOfZThreadsInSimulation);

        vector<vector<vector<thread*>>> Threads(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<thread*>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<thread*>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        auto GenerateNStepsOfSimulationForWholeCellSpaceInThreads = [NumberOfStepsOutside, &SyncPoint, this](const UnsignedInt NumberOfStepsInside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
        {
            for (UnsignedInt StepOutside = 1; StepOutside <= NumberOfStepsOutside; StepOutside++)
            {
                CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(NumberOfStepsInside, StepOutside, ThreadXIndex, ThreadYIndex, ThreadZIndex);

                SyncPoint.arrive_and_wait();

                if (CurrentThreadIndex == 1)
                    ExchangeParticlesBetweenThreads(StepOutside);

                SyncPoint.arrive_and_wait();
            }
        };

        LoggersManagerObject.Log(STREAM("START THREADS"));

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                    Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1] = new thread(GenerateNStepsOfSimulationForWholeCellSpaceInThreads, NumberOfStepsInside, ThreadXIndex, ThreadYIndex, ThreadZIndex);

        for (auto& ThreadX : Threads)
            for (auto& ThreadY : ThreadX)
                for (auto& ThreadZ : ThreadY)
                {
                    ThreadZ->join();
                    delete ThreadZ;
                }

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        string ResultText = "Execution in threads for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in threads")));

        LoggersManagerObject.Log(STREAM("END THREADS"));

        LoggersManagerObject.Log(STREAM("AddedParticlesInReactions =  " + to_string(AddedParticlesInReactions)));
        LoggersManagerObject.Log(STREAM("ErrorCounter = " + to_string(ErrorCounter)));
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationSpace::GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, bool PrintTime)
{
    try
    {
        FirstSendParticlesForThreads(false, true);
        GenerateNStepsOfSimulationForWholeCellSpaceInThreads(NumberOfStepsOutside, NumberOfStepsInside);
        GatherParticlesForThreadsInMainParticles();
    }
    CATCH("generate n steps of simulation with sending particles to threads and gathering particles to main threads for whole cell space")
}

void CellEngineSimulationSpace::CheckParticlesCenters() const
{
    try
    {
        UnsignedInt NumberOfZeroCenterParticles = 0;
        for (const auto& ParticleObject : Particles)
        {
            cout << "ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z << endl;

            if (ParticleObject.second.Center.X > CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension || ParticleObject.second.Center.Y > CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension || ParticleObject.second.Center.Z > CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                LoggersManagerObject.Log(STREAM("Wrong Particle Center -> ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));
            if (ParticleObject.second.Center.X == 0 || ParticleObject.second.Center.Y == 0 || ParticleObject.second.Center.Z == 0)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ZERO -> ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));
                NumberOfZeroCenterParticles++;
            }
        }
        LoggersManagerObject.Log(STREAM("All Particle Centers Checked. Number Of Zero Particles = " << NumberOfZeroCenterParticles));
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}