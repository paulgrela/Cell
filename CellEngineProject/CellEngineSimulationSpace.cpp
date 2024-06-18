
#include <set>
#include <map>
#include <algorithm>

#include "FileUtils.h"
#include "DoublyLinkedList.h"

#include "CellEngineConstants.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineChemicalReactionsManager.h"

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
    CATCH("generating one step of diffusion")
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

        for (const auto& NeighbourParticleIndexObjectToWrite : ParticlesSortedByCapacityFoundInProximity)
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

void CellEngineSimulationSpace::GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        ElectricChargeType NeighbourPoints[3][3][3];

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            if (GetParticleFromIndex(ParticleIndex).ElectricCharge != 0)
            {
                Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);

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

                UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(ParticleObject, &NeighbourPoints, MultiplyElectricChargeFactor);

                vector<vector3<SignedInt>> MoveVectors;
                UpdateNeighbourPointsForChosenVoxel([&MoveVectors](SignedInt X, SignedInt Y, SignedInt Z){ MoveVectors.emplace_back(X - 1, Y - 1, Z - 1); });

                vector<int> DiscreteDistribution;
                DiscreteDistribution.reserve(9);

                UpdateNeighbourPointsForChosenVoxel([&NeighbourPoints, &DiscreteDistribution](SignedInt X, SignedInt Y, SignedInt Z){ DiscreteDistribution.emplace_back(NeighbourPoints[X][Y][Z]); });
                #ifdef SIMULATION_DETAILED_LOG
                UnsignedInt NumberOfElement = 0;
                UpdateNeighbourPointsForChosenVoxel([&NeighbourPoints, &MoveVectors, &NumberOfElement](SignedInt X, SignedInt Y, SignedInt Z){ LoggersManagerObject.Log(STREAM("Element[" << NumberOfElement << "] = " << to_string(NeighbourPoints[X][Y][Z]) + " for (X,Y,Z) = (" << to_string(MoveVectors[NumberOfElement].X) << "," << to_string(MoveVectors[NumberOfElement].Y) << "," << to_string(MoveVectors[NumberOfElement].Z) << ")")); NumberOfElement++; });
                #endif

                discrete_distribution<int> UniformDiscreteDistributionMoveParticleDirectionObject(DiscreteDistribution.begin(), DiscreteDistribution.end());

                UnsignedInt RandomMoveVectorIndex = UniformDiscreteDistributionMoveParticleDirectionObject(mt64R);
                MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(ParticleObject, MoveVectors[RandomMoveVectorIndex].X, MoveVectors[RandomMoveVectorIndex].Y, MoveVectors[RandomMoveVectorIndex].Z, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

                #ifdef SIMULATION_DETAILED_LOG
                LoggersManagerObject.Log(STREAM("Random Index = " << to_string(RandomMoveVectorIndex) << " " << to_string(MoveVectors[RandomMoveVectorIndex].X) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Y) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Z) << endl));
                #endif
            }

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating one step of electric diffusion")
}

tuple<vector<pair<UniqueIdInt, UnsignedInt>>, bool> CellEngineSimulationSpace::ChooseParticlesForReactionFromAllParticlesInProximity(const ChemicalReaction& ReactionObject)
{
    bool AllAreZero = false;

    vector<pair<UniqueIdInt, UnsignedInt>> NucleotidesIndexesChosenForReaction, ParticlesIndexesChosenForReaction, AllParticlesIndexesChosenForReaction;

    vector<UnsignedInt> ReactantsCounters(ReactionObject.Reactants.size());

    try
    {
        for (UnsignedInt ReactantIndex = 0; ReactantIndex < ReactionObject.Reactants.size(); ReactantIndex++)
            ReactantsCounters[ReactantIndex] = ReactionObject.Reactants[ReactantIndex].Counter;

        for (const auto& ParticleObjectIndex : ParticlesSortedByCapacityFoundInProximity)
        {
            auto& ParticleObjectTestedForReaction = GetParticleFromIndex(ParticleObjectIndex);

            LoggersManagerObject.Log(STREAM("ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z)));

            vector<ParticleKindForChemicalReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return ParticleKindForReactionObjectParam.EntityId == ParticleObjectTestedForReaction.EntityId && CompareFitnessOfParticle(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
            else
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction, this](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsSpecialDNA(ParticleKindForReactionObjectParam.EntityId) && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByString, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

            auto PositionInReactants = ReactantIterator - ReactionObject.Reactants.begin();

            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == true)
                if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == false)
                    NucleotidesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == true)
                ParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0)
            {
                AllParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);
                LoggersManagerObject.Log(STREAM("CHOSEN ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z) << endl));
                ReactantsCounters[PositionInReactants]--;
            }

            AllAreZero = all_of(ReactantsCounters.begin(), ReactantsCounters.end(), [this](const UnsignedInt& Counter){ return Counter == 0; });
            if (AllAreZero == true)
            {
                LoggersManagerObject.Log(STREAM("ALL ARE ZERO"));
                break;
            }
            LoggersManagerObject.Log(STREAM(""));
        }

        if (AllAreZero == true || (AllAreZero == false && (ReactionObject.ReactionIdNum == 30 || ReactionObject.ReactionIdNum == 80 || ReactionObject.ReactionIdNum == 70)))
            if (ReactionObject.SpecialReactionFunction != nullptr)
                ReactionObject.SpecialReactionFunction(this, AllParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
    }
    CATCH("choosing particles for reaction from all particles in proximity")

    if (AllAreZero == true)
    {
        LoggersManagerObject.Log(STREAM("ALL ARE ZERO AT END = " << to_string(ParticlesIndexesChosenForReaction.size())));
        return { ParticlesIndexesChosenForReaction, true };
    }
    else
        return { vector<pair<UniqueIdInt, UnsignedInt>>(), false };
}

bool CellEngineSimulationSpace::MakeChemicalReaction(ChemicalReaction& ReactionObject)
{
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
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.clear();

            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

            if (CenterIndex < Centers.size())
            {
                for (const auto& ParticleKindVoxel : ParticleKindObjectForProduct.ListOfVoxels)
                {
                    vector3_64 NewVoxel(Centers[CenterIndex].X - ParticleKindObjectForProduct.XSizeDiv2 + ParticleKindVoxel.X, Centers[CenterIndex].Y - ParticleKindObjectForProduct.YSizeDiv2 + ParticleKindVoxel.Y, Centers[CenterIndex].Z - ParticleKindObjectForProduct.ZSizeDiv2 + ParticleKindVoxel.Z);

                    FillParticleElementInSpace(ParticleIndex, NewVoxel);

                    LoggersManagerObject.Log(STREAM("New Centers From Product Added X = " << to_string(NewVoxel.X) << " Y = " << to_string(NewVoxel.Y) << " Z = " << to_string(NewVoxel.Z) << endl));
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

    return true;
};

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticles(const UnsignedInt NumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, ParticlesKindsFoundInProximity.size() - 1);

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(std::next(std::begin(ParticlesKindsFoundInProximity), static_cast<int>(UniformDistributionObjectUint64t(mt64R)))->first);

            if (CellEngineUseful::IsSpecialDNA(RandomParticlesTypes.back()) == true)
                RandomParticlesTypes.back() = CellEngineConfigDataObject.DNAIdentifier;

            LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
        }
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

bool CellEngineSimulationSpace::IsChemicalReactionPossible(const ChemicalReaction& ReactionObject)
{
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForChemicalReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsFoundInProximity[CellEngineUseful::IfSpecialDNAThenReturnNormalDNACode(ReactionReactant.EntityId)]; });
}

void CellEngineSimulationSpace::PrepareRandomReaction()
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMainRandomCondition_Uint64t(0, 1);
        if (UniformDistributionObjectMainRandomCondition_Uint64t(mt64R) == 0)
            return;

        LoggersManagerObject.Log(STREAM("Looking for particles in proximity"));
    }
    CATCH("preparing random reaction")
}

void CellEngineSimulationSpace::FindAndExecuteRandomReaction()
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectNumberOfReactants_Uint64t(2, 2);

        UnsignedInt NumberOfTries = 0;
        while (NumberOfTries <= 100)
        {
            NumberOfTries++;

            UnsignedInt NumberOfReactants = UniformDistributionObjectNumberOfReactants_Uint64t(mt64R);

            if (TryToMakeRandomChemicalReaction(NumberOfReactants) == true)
                break;
        }
    }
    CATCH("finding and executing random reaction")
}

void CellEngineSimulationSpace::FindAndExecuteChosenReaction(const UnsignedInt ReactionId)
{
    try
    {
        auto ReactionIter = CellEngineChemicalReactionsManagerObject.ChemicalReactionsPosFromId.find(ReactionId);
        if (ReactionIter != CellEngineChemicalReactionsManagerObject.ChemicalReactionsPosFromId.end())
        {
            auto& ReactionObject = CellEngineChemicalReactionsManagerObject.ChemicalReactions[ReactionIter->second];

            bool IsPossible = IsChemicalReactionPossible(ReactionObject);
            if (IsPossible == true)
            {
                LoggersManagerObject.Log(STREAM("CHOSEN REACTION POSSIBLE" << endl));

                if (MakeChemicalReaction(ReactionObject) == false)
                    LoggersManagerObject.Log(STREAM("Chosen reaction not executed!"));
            }
            else
                LoggersManagerObject.Log(STREAM("Chosen reaction impossible!"));
        }
        else
            LoggersManagerObject.Log(STREAM("Chosen reaction Id not found!"));
    }
    CATCH("finding and executing chosen reaction")
}

void CellEngineSimulationSpace::GenerateRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            FindAndExecuteRandomReaction();
    }
    CATCH("generating random reaction for selected voxel space")
}

void CellEngineSimulationSpace::GenerateChosenReactionForSelectedSpace(UnsignedInt ReactionId, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        if (FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            if (ReactionId != 0)
                FindAndExecuteChosenReaction(ReactionId);
    }
    CATCH("generating random reaction for particle")
}

void CellEngineSimulationSpace::GenerateRandomReactionForParticle(Particle& ParticleObject)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(ParticleObject, 20) == true)
            FindAndExecuteRandomReaction();
    }
    CATCH("generating random reaction for particle")
}

void CellEngineSimulationSpace::GenerateRandomReactionsForAllParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.SelectedForReaction == false)
                GenerateRandomReactionForParticle(ParticleObject.second);
    }
    CATCH("generating random reactions for all particles")
}

void CellEngineSimulationSpace::GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto ParticlesIterator = Particles.find(ParticleIndex);
            if (ParticlesIterator != Particles.end())
                GenerateRandomReactionForParticle(ParticlesIterator->second);
        }
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineSimulationSpace::GenerateOneStepOfRandomReactionsForOneParticle(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);
        GenerateRandomReactionForParticle(GetParticleFromIndex(StartParticleIndexParam + 4));
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineSimulationSpace::SaveReactionsStatisticsToFile()
{
    try
    {
        for (const auto& ReactionData : SavedReactionsMap[SimulationStepNumber - 1])
            LoggersManagerObject.LogStatistics(STREAM("REACTION ID = " << ReactionData.second.ReactionId << " REACTION NAME = " << CellEngineChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactionName << " REACTION ID_STR = #" << CellEngineChemicalReactionsManagerObject.GetReactionFromNumId(ReactionData.second.ReactionId).ReactionIdStr << "# REACTION COUNTER = " << ReactionData.second.Counter));
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

void CellEngineSimulationSpace::GenerateChosenReactionsForWholeCellSpace(const UnsignedInt ReactionId, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                    GenerateChosenReactionForSelectedSpace(ReactionId, PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for whole cell space")
}

void CellEngineSimulationSpace::GenerateRandomReactionsForWholeCellSpace(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt PosX = XStartParam; PosX < XSizeParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam; PosY < YSizeParam; PosY += YStepParam)
                for (UnsignedInt PosZ = ZStartParam; PosZ < ZSizeParam; PosZ += ZStepParam)
                    GenerateRandomReactionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for whole cell space")
}