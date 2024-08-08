
#include <set>
#include <map>
#include <algorithm>

#include "FileUtils.h"
#include "Combinatorics.h"
#include "DoublyLinkedList.h"

#include "CellEngineUseful.h"
#include "CellEngineConstants.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineChemicalReactionsManager.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

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

void CellEngineSimulationSpace::GenerateOneStepOfDiffusionForSelectedSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-1, 1);

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
        for (auto& ParticleInProximityIndex : ParticlesSortedByCapacityFoundInProximity)
            if (CellEngineUseful::IsDNA(GetParticleFromIndex(ParticleInProximityIndex).EntityId) == false)
                MoveParticleByVectorIfSpaceIsEmpty(GetParticleFromIndex(ParticleInProximityIndex), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R));
    }
    CATCH("generating one step of diffusion for selected space")
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

        auto ParticlesSortedByCapacityFoundInProximityCopy(ParticlesSortedByCapacityFoundInProximity);

        for (auto& ParticleInProximityIndex : ParticlesSortedByCapacityFoundInProximityCopy)
            if (CellEngineUseful::IsDNA(GetParticleFromIndex(ParticleInProximityIndex).EntityId) == false)
                GenerateOneStepOfElectricDiffusionForOneParticle(TypeOfLookingForParticles, AdditionalSpaceBoundFactor, MultiplyElectricChargeFactor, ParticleInProximityIndex, &NeighbourPoints, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

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
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction, this](const ParticleKindForChemicalReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsDNA(ParticleKindForReactionObjectParam.EntityId) == true && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByString, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

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

std::vector<UnsignedInt> CellEngineSimulationSpace::GetRandomParticlesVersion3(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, ParticlesKindsFoundInProximity.size() - 1);
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
                RandomParticlesTypes.emplace_back(std::next(std::begin(ParticlesKindsFoundInProximity), static_cast<int>(ReactantNumberBitValuePos))->first);

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
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, ParticlesKindsFoundInProximity.size() - 1);

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(std::next(std::begin(ParticlesKindsFoundInProximity), static_cast<int>(UniformDistributionObjectUint64t(mt64R)))->first);

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
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForChemicalReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsFoundInProximity[ReactionReactant.EntityId]; });
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
        for (const auto& ParticleKindFoundInProximityObject : ParticlesKindsFoundInProximity)
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

        for (const auto& ParticleKindFoundInProximityObject : ParticlesKindsFoundInProximity)
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

void CellEngineSimulationSpace::GenerateRandomReactionForSelectedSpace(UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            FindAndExecuteRandomReaction(min(ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
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
            FindAndExecuteRandomReaction(min(ParticlesKindsFoundInProximity.size(), ChemicalReactionsManagerObject.MaxNumberOfReactants));
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
            if (auto ParticlesIterator = Particles.find(ParticleIndex); ParticlesIterator != Particles.end())
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

void CellEngineSimulationSpace::GenerateDiffusionForBigPartOfCellSpace(const UnsignedInt SizeNMultiplyFactor, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt PosX = XStartParam - SizeNMultiplyFactor * XStepParam; PosX <= XStartParam + SizeNMultiplyFactor * XStepParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam - SizeNMultiplyFactor * YStepParam; PosY <= YStartParam + SizeNMultiplyFactor * YStepParam; PosY += YStepParam)
                for (UnsignedInt PosZ = ZStartParam - SizeNMultiplyFactor * ZStepParam; PosZ <= ZStartParam + SizeNMultiplyFactor * ZStepParam; PosZ += ZStepParam)
                    GenerateOneStepOfDiffusionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating diffusion for big part of cell")
}

void CellEngineSimulationSpace::GenerateRandomReactionsForBigPartOfCellSpace(const UnsignedInt SizeNMultiplyFactor, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        for (UnsignedInt PosX = XStartParam - SizeNMultiplyFactor * XStepParam; PosX <= XStartParam + SizeNMultiplyFactor * XStepParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam - SizeNMultiplyFactor * YStepParam; PosY <= YStartParam + SizeNMultiplyFactor * YStepParam; PosY += YStepParam)
                for (UnsignedInt PosZ = ZStartParam - SizeNMultiplyFactor * ZStepParam; PosZ <= ZStartParam + SizeNMultiplyFactor * ZStepParam; PosZ += ZStepParam)
                    GenerateRandomReactionForSelectedSpace(PosX, PosY, PosZ, XStepParam, YStepParam, ZStepParam);

        CheckConditionsToIncSimulationStepNumberForStatistics();

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating random reactions for big part of cell")
}
