

#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"
#include "DestinationPlatform.h"

#include "CellEngineAtom.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"

#include "TerminalColorsUtils.h"

#include "CellEngineDataBuilderForVoxelSimulationSpace.h"

using namespace std;

inline void GetRangeOfParticlesForRandomParticles(UniqueIdInt& StartParticleIndexParam, UniqueIdInt& EndParticleIndexParam, UniqueIdInt MaxParticleIndex)
{
    if (EndParticleIndexParam == 0)
    {
        StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
        EndParticleIndexParam = MaxParticleIndex;
    }
}

SimulationSpaceVoxel CellEngineVoxelSimulationSpaceForOuterClass::GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z)
{
    return GetSpaceVoxel(X, Y, Z);
}

Particle& CellEngineVoxelSimulationSpaceForOuterClass::GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

[[nodiscard]] float CellEngineVoxelSimulationSpaceD1::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
{
    return static_cast<float>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace;
};

[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpaceD1::ConvertToSpaceCoordinate(double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace) + (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2);
};

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : Particles(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam), CellEngineAllTypesOfChemicalReactionsInVoxelSpace(ParticlesParam), CellEngineAllTypesOfParticlesGenerator(ParticlesParam)
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        SetValueToVoxelsForCuboidSelectedSpace(nullptr, GetZeroSimulationSpaceVoxel(), 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension);

        Genomes.resize(2);
        GenomesLines.resize(1);
    }
    CATCH("execution of constructor of voxel simulation space")
}

CellEngineVoxelSimulationSpace::~CellEngineVoxelSimulationSpace()
{
    try
    {
        free(SpacePointer);
    }
    CATCH("execution of destructor of voxel simulation space")
}

[[nodiscard]] stringstream CellEngineVoxelSimulationSpaceD1::PrintSpaceMinMaxValues() const
{
    stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMin) << " ][ Xmax = " << to_string(XMax) << " ][ Ymin = " << to_string(YMin) << " ][ Ymax = " << to_string(YMax) << " ][ Zmin = " << to_string(ZMin) << " ][ Zmax = " << to_string(XMax) << " ] " << endl;
    return ss;
}

void CellEngineVoxelSimulationSpaceD1::SetAtomInVoxelSimulationSpace(const UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom)
{
    try
    {
        UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates(PosX, PosY, PosZ, XMin, XMax, YMin, YMax, ZMin, ZMax);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
            SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, PosX, PosY, PosZ);
    }
    CATCH("setting atom in voxel simulation space")
}

void CellEngineVoxelSimulationSpaceD1::ClearSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        SetValueToVoxelsForCuboidSelectedSpace(nullptr, GetZeroSimulationSpaceVoxel(), StartXPosParam, StartYPosParam, StartZPosParam, StepXParam, StepYParam, StepZParam, SizeXParam, SizeYParam, SizeZParam);
    }
    CATCH("clearing selected space")
}

void CellEngineVoxelSimulationSpaceD1::ClearWholeVoxelSpace()
{
    try
    {
        SetValueToVoxelsForCuboidSelectedSpace(nullptr, GetZeroSimulationSpaceVoxel(), 0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension);
    }
    CATCH("clearing whole voxel space")
};

void CellEngineVoxelSimulationSpaceD1::ClearVoxelSpaceAndParticles()
{
    try
    {
        Particles.clear();
        ClearWholeVoxelSpace();
    }
    CATCH("clearing voxel space and particles")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(-1, 1);

        vector<vector3_16> NewVoxelsForParticle;

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
            MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(GetParticleFromIndex(ParticleIndex), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), UniformDistributionObjectMoveParticleDirection_int64t(mt64R), StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);
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

void CellEngineVoxelSimulationSpace::UpdateProbabilityOfMoveFromElectricInteractionForSelectedParticle(Particle& ParticleObject, ElectricChargeType (*NeighbourPoints)[3][3][3], const double MultiplyElectricChargeFactor)
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

void CellEngineVoxelSimulationSpace::GenerateOneStepOfElectricDiffusionForSelectedRangeOfParticles(const TypesOfLookingForParticlesInProximity TypeOfLookingForParticles, const UnsignedInt AdditionalSpaceBoundFactor, const double MultiplyElectricChargeFactor, UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        CellEngineUseful::SwitchOffLogs();

        vector<vector3_16> NewVoxelsForParticle;

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
                    case TypesOfLookingForParticlesInProximity::FromChosenParticleAsCenter : FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(false, ParticleObject.Center.X - ParticleKindObject.XSizeDiv2 - AdditionalSpaceBoundFactor, ParticleObject.Center.Y - ParticleKindObject.YSizeDiv2 - AdditionalSpaceBoundFactor, ParticleObject.Center.Z - ParticleKindObject.ZSizeDiv2 - AdditionalSpaceBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalSpaceBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalSpaceBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalSpaceBoundFactor); break;
                    case TypesOfLookingForParticlesInProximity::InChosenVoxelSpace : FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(false, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam); break;
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
                MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(ParticleObject, MoveVectors[RandomMoveVectorIndex].X, MoveVectors[RandomMoveVectorIndex].Y, MoveVectors[RandomMoveVectorIndex].Z, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

                #ifdef SIMULATION_DETAILED_LOG
                LoggersManagerObject.Log(STREAM("Random Index = " << to_string(RandomMoveVectorIndex) << " " << to_string(MoveVectors[RandomMoveVectorIndex].X) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Y) << " " << to_string(MoveVectors[RandomMoveVectorIndex].Z) << endl));
                #endif
            }

        CellEngineUseful::SwitchOnLogs();
    }
    CATCH("generating one step of electric diffusion")
}

tuple<vector<pair<UniqueIdInt, UnsignedInt>>, bool> CellEngineVoxelSimulationSpace::ChooseParticlesForReactionFromAllParticlesInProximity(const Reaction& ReactionObject)
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

            vector<ParticleKindForReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction](const ParticleKindForReaction& ParticleKindForReactionObjectParam){ return ParticleKindForReactionObjectParam.EntityId == ParticleObjectTestedForReaction.EntityId && CompareFitnessOfParticle(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
            else
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction, this](const ParticleKindForReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsSpecialDNA(ParticleKindForReactionObjectParam.EntityId) && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByVectorLoop, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

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

        if (AllAreZero == true || (AllAreZero == false && (ReactionObject.Id == 30 || ReactionObject.Id == 80 || ReactionObject.Id == 70)))
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

bool CellEngineVoxelSimulationSpace::MakeChemicalReaction(Reaction& ReactionObject)
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

                    SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(&GetParticleFromIndex(ParticleIndex).ListOfVoxels, ParticleIndex, NewVoxel.X, NewVoxel.Y, NewVoxel.Z);

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
    }
    CATCH("making reaction")

    return true;
};

std::vector<UnsignedInt> CellEngineVoxelSimulationSpace::GetRandomParticles(const UnsignedInt NumberOfReactants)
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

bool CellEngineVoxelSimulationSpace::IsChemicalReactionPossible(const Reaction& ReactionObject)
{
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsFoundInProximity[CellEngineUseful::IfSpecialDNAThenReturnNormalDNACode(ReactionReactant.EntityId)]; });
}

void CellEngineVoxelSimulationSpace::PrepareRandomReaction()
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

void CellEngineVoxelSimulationSpace::FindAndExecuteRandomReaction()
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

void CellEngineVoxelSimulationSpace::FindAndExecuteChosenReaction(const UnsignedInt ReactionId)
{
    try
    {
        auto ReactionIter = ReactionsPosFromId.find(ReactionId);
        if (ReactionIter != ReactionsPosFromId.end())
        {
            auto& ReactionObject = Reactions[ReactionIter->second];

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

void CellEngineVoxelSimulationSpace::GenerateRandomReactionForSelectedVoxelSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            FindAndExecuteRandomReaction();
    }
    CATCH("generating random reaction for selected voxel space")
}

void CellEngineVoxelSimulationSpace::GenerateChosenReactionForSelectedVoxelSpace(const UnsignedInt ReactionId, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        if (FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(true, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            if (ReactionId != 0)
                FindAndExecuteChosenReaction(ReactionId);
    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionForParticle(Particle& ParticleObject)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(ParticleObject, 20) == true)
            FindAndExecuteRandomReaction();
    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionsForAllParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.SelectedForReaction == false)
                GenerateRandomReactionForParticle(ParticleObject.second);
    }
    CATCH("generating random reactions for all particles")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
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

void CellEngineVoxelSimulationSpace::GenerateOneStepOfRandomReactionsForOneParticle(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);
        GenerateRandomReactionForParticle(GetParticleFromIndex(StartParticleIndexParam + 4));
    }
    CATCH("generating one step of random reactions for selected range of particles")
}