
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>

#include "FileUtils.h"
#include "DestinationPlatform.h"

#include "CellEngineAtom.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"

using namespace std;

[[nodiscard]] float CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
{
    return static_cast<float>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace;
};

[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace) + (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2);
};

inline SimulationSpaceVoxel& CellEngineVoxelSimulationSpace::GetSpaceVoxel(UnsignedInt x, UnsignedInt y, UnsignedInt z)
{
    return CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? (*static_cast<Space_2048_2048_2048*>(SpacePointer))[x][y][z] : (*static_cast<Space_1024_1024_1024*>(SpacePointer))[x][y][z];
}

inline Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndex(const UniqueIdInt ParticleIndex)
{
    return Particles[ParticleIndex];
}

inline SimulationSpaceVoxel GetZeroSimulationSpaceVoxel()
{
    return 0;
}

SimulationSpaceVoxel CellEngineVoxelSimulationSpace::GetSimulationSpaceVoxel(UnsignedInt X, UnsignedInt Y, UnsignedInt Z)
{
    return GetSpaceVoxel(X, Y, Z);
}

Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndexInSimulationSpaceVoxel(const UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                    GetSpaceVoxel(PosX, PosZ, PosY) = GetZeroSimulationSpaceVoxel();
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

void CellEngineVoxelSimulationSpace::InitiateFreeParticleIndexes()
{
    try
    {
        for (UnsignedInt FreeIndex = MaxParticleIndex + 100'000'000; FreeIndex >= MaxParticleIndex + 1; FreeIndex--)
            FreeIndexesOfParticles.push(FreeIndex);
    }
    CATCH("initiating free particle indexes")
}

void CellEngineVoxelSimulationSpace::GetMinMaxCoordinatesForAllParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            GetMinMaxCoordinatesForParticle(ParticleObject.second);
    }
    CATCH("getting min max coordinates for all particles")
}

void CellEngineVoxelSimulationSpace::GetMinMaxCoordinatesForParticle(Particle& ParticleObject)
{
    try
    {
        UnsignedInt ParticleXMin{}, ParticleXMax{}, ParticleYMin{}, ParticleYMax{}, ParticleZMin{}, ParticleZMax{};
        ParticleXMin = ParticleYMin = ParticleZMin = 10000;
        ParticleXMax = ParticleYMax = ParticleZMax = 0;

        for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
            GetMinMaxOfCoordinates(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);

        ParticleObject.SetCenterCoordinates(ParticleXMin + (ParticleXMax - ParticleXMin) / 2, ParticleYMin + (ParticleYMax - ParticleYMin) / 2, ParticleZMin + (ParticleZMax - ParticleZMin) / 2);

        auto& ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);
        if (ParticleKindObject.ListOfVoxels.empty() == true)
        {
            for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
                ParticleKindObject.ListOfVoxels.emplace_back(VoxelCoordinates.X - ParticleXMin, VoxelCoordinates.Y - ParticleYMin, VoxelCoordinates.Z - ParticleZMin);
            ParticleKindObject.XSizeDiv2 = (ParticleXMax - ParticleXMin) / 2;
            ParticleKindObject.YSizeDiv2 = (ParticleYMax - ParticleYMin) / 2;
            ParticleKindObject.ZSizeDiv2 = (ParticleZMax - ParticleZMin) / 2;
        }
    }
    CATCH("getting min max coordinates for one particle")
}

void CellEngineVoxelSimulationSpace::PreprocessData()
{
    try
    {
        InitiateFreeParticleIndexes();
        GetMinMaxCoordinatesForAllParticles();
    }
    CATCH("preprocessing data for voxel simulation space")
}

void CellEngineVoxelSimulationSpace::SetStartValuesForSpaceMinMax()
{
    XMin = YMin = ZMin = 10000;
    XMax = YMax = ZMax = 0;
}

inline void CellEngineVoxelSimulationSpace::GetMinMaxOfCoordinates(const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam)
{
    try
    {
        XMinParam = min(PosX, XMinParam);
        XMaxParam = max(PosX, XMaxParam);
        YMinParam = min(PosY, YMinParam);
        YMaxParam = max(PosY, YMaxParam);
        ZMinParam = min(PosZ, ZMinParam);
        ZMaxParam = max(PosZ, ZMaxParam);
    }
    CATCH("getting min max of coordinates")
}

[[nodiscard]] stringstream CellEngineVoxelSimulationSpace::PrintSpaceMinMaxValues() const
{
    stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMin) << " ][ Xmax = " << to_string(XMax) << " ][ Ymin = " << to_string(YMin) << " ][ Ymax = " << to_string(YMax) << " ][ Zmin = " << to_string(ZMin) << " ][ Zmax = " << to_string(XMax) << " ] " << endl;
    return ss;
}

void CellEngineVoxelSimulationSpace::CountStatisticsOfVoxelSimulationSpace()
{
    try
    {
        for(UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for(UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for(UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                {
                    if (GetParticleFromIndex(GetSpaceVoxel(PosX, PosY, PosZ)).EntityId != 0)
                        SumOfNotEmptyVoxels++;
                }
    }
    CATCH("counting statistics of voxel simulation space")
}

void CellEngineVoxelSimulationSpace::SetAtomInVoxelSimulationSpace(const UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom)
{
    try
    {
        UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates(PosX, PosY, PosZ, XMin, XMax, YMin, YMax, ZMin, ZMax);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
        {
            GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
        }
    }
    CATCH("setting atom in voxel simulation space")
}

UniqueIdInt CellEngineVoxelSimulationSpace::AddNewParticle(const Particle& ParticleParam)
{
    Particles[ParticleParam.Index] = ParticleParam;
    return MaxParticleIndex = ParticleParam.Index;
}

void CellEngineVoxelSimulationSpace::AddReaction(const Reaction& ReactionParam)
{
    try
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsIdByString.insert(make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
    CATCH("adding reaction")
}

UniqueIdInt CellEngineVoxelSimulationSpace::GetFreeIndexesOfParticleSize()
{
    return FreeIndexesOfParticles.size();
}

UniqueIdInt CellEngineVoxelSimulationSpace::GetNewFreeIndexOfParticle()
{
    if (FreeIndexesOfParticles.empty() == false)
    {
        UniqueIdInt FreeIndexOfParticle = FreeIndexesOfParticles.top();
        FreeIndexesOfParticles.pop();
        return FreeIndexOfParticle;
    }
    else
    {
        LoggersManagerObject.Log(STREAM("Lack of new free indexes of particles"));
        return MaxParticleIndex + 1;
    }
}

void CellEngineVoxelSimulationSpace::AddBasicParticlesKindsAndReactions()
{
    try
    {
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", 0});
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen6", "06", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 9, "Eten", "CH2CH2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 12, "DNA", "CGATATTAAATAGGGCCT", 0 });

        AddReaction(Reaction("C6H12O6 + O6 + ", { { 1, 1, true }, { 2, 2, true } }, { { 0, 6, true }, { 0, 6, true } }));
        AddReaction(Reaction("CH2CH2 + H2O + ", { { 4, 1, true }, { 0, 1, true } }, { { 5, 1, true } }));
        AddReaction(Reaction("CH3CHCH2 + HX + ", { { 6, 1, true }, { 7, 1, true } }, { { 8, 1, true } }));
        AddReaction(Reaction("CH2CH2 + O + ", { { 9, 1, true }, { 11, 1, true } }, { { 10, 1, true } }));
    }
    CATCH("adding particles kinds and reactions")
};

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        LoggersManagerObject.Log(STREAM("PARTICLES RANGE = " << MaxParticleIndex << " " << (MaxParticleIndex + NumberOfRandomParticles) << " " << NumberOfRandomParticles));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectObjectOfParticle_Uint64t(0, NumberOfRandomParticles);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);

        AddBasicParticlesKindsAndReactions();

        vector<UnsignedInt> LocalNewParticlesIndexes;

        for (UniqueIdInt ParticleNumber = 1; ParticleNumber <= NumberOfRandomParticles; ParticleNumber++)
            LocalNewParticlesIndexes.emplace_back(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), UniformDistributionObjectObjectOfParticle_Uint64t(mt64R), 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));

        vector<vector3_16> FilledVoxelsForRandomParticle;

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    GetSpaceVoxel(PosX, PosY, PosZ) = GetZeroSimulationSpaceVoxel();

        for (auto& LocalNewParticleIndex : LocalNewParticlesIndexes)
        {
            UnsignedInt RandomPosX = UniformDistributionObjectX_Uint64t(mt64R);
            UnsignedInt RandomPosY = UniformDistributionObjectY_Uint64t(mt64R);
            UnsignedInt RandomPosZ = UniformDistributionObjectZ_Uint64t(mt64R);

            UnsignedInt RandomSizeOfParticle = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);

            FilledVoxelsForRandomParticle.clear();

            if (RandomPosX + RandomSizeOfParticle < StartXPosParam + SizeXParam && RandomPosY + RandomSizeOfParticle < StartYPosParam + SizeYParam && RandomPosZ + RandomSizeOfParticle < StartZPosParam + SizeZParam)
                for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + RandomSizeOfParticle; PosX++)
                    for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + RandomSizeOfParticle; PosY++)
                        for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + RandomSizeOfParticle; PosZ++)
                        {
                            if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
                            {
                                FilledVoxelsForRandomParticle.emplace_back(PosX, PosY, PosZ);
                                GetSpaceVoxel(PosX, PosY, PosZ) = LocalNewParticleIndex;
                            }
                            else
                            {
                                for (auto& VoxelForRandomParticle : FilledVoxelsForRandomParticle)
                                    GetSpaceVoxel(VoxelForRandomParticle.X, VoxelForRandomParticle.Y, VoxelForRandomParticle.Z) = GetZeroSimulationSpaceVoxel();

                                FilledVoxelsForRandomParticle.clear();

                                goto NextRandomParticleOutsideLoopLabel;
                            }
                        }
            NextRandomParticleOutsideLoopLabel:;

            if (FilledVoxelsForRandomParticle.empty() == false)
            {
                LoggersManagerObject.Log(STREAM("LocalNewParticleIndex = " << LocalNewParticleIndex));
                Particles[LocalNewParticleIndex].ListOfVoxels = FilledVoxelsForRandomParticle;
            }
        }

        for (auto& LocalNewParticleIndex : LocalNewParticlesIndexes)
            GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
    }
    CATCH("generating random particles in selected space")
}

inline void GetRangeOfParticlesForRandomParticles(UniqueIdInt& StartParticleIndexParam, UniqueIdInt& EndParticleIndexParam, UniqueIdInt MaxParticleIndex)
{
    if (EndParticleIndexParam == 0)
    {
        StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
        EndParticleIndexParam = MaxParticleIndex;
    }
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectSizeOfParticle_int64t(-1, 1);

        vector<vector3_16> NewVoxelsForParticle;

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto &LocalParticleObject = Particles[ParticleIndex];

            for (auto &VoxelForParticle: LocalParticleObject.ListOfVoxels)
                GetSpaceVoxel(VoxelForParticle.X, VoxelForParticle.Y, VoxelForParticle.Z) = GetZeroSimulationSpaceVoxel();

            SignedInt ShiftX = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftY = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftZ = UniformDistributionObjectSizeOfParticle_int64t(mt64R);

            NewVoxelsForParticle.clear();
            bool Collision = false;

            for (auto &VoxelForParticle: LocalParticleObject.ListOfVoxels)
                if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY,VoxelForParticle.Z + ShiftZ) == 0 && VoxelForParticle.X + ShiftX >= StartXPosParam && VoxelForParticle.X + ShiftX < StartXPosParam + SizeXParam && VoxelForParticle.Y + ShiftY >= StartYPosParam && VoxelForParticle.Y + ShiftY < StartYPosParam + SizeYParam && VoxelForParticle.Z + ShiftZ >= StartZPosParam && VoxelForParticle.Z + ShiftZ < StartZPosParam + SizeZParam)
                {
                    GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) = ParticleIndex;
                    NewVoxelsForParticle.emplace_back(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ);
                }
                else
                {
                    for (auto &NewVoxelForParticle: NewVoxelsForParticle)
                        GetSpaceVoxel(NewVoxelForParticle.X, NewVoxelForParticle.Y, NewVoxelForParticle.Z) = GetZeroSimulationSpaceVoxel();

                    for (auto &OldVoxelForParticle: LocalParticleObject.ListOfVoxels)
                        GetSpaceVoxel(OldVoxelForParticle.X, OldVoxelForParticle.Y, OldVoxelForParticle.Z) = ParticleIndex;

                    Collision = true;
                    break;
                }

            if (Collision == false)
                LocalParticleObject.ListOfVoxels = NewVoxelsForParticle;
        }
    }
    CATCH("generating one step of diffusion")
}

std::vector<UnsignedInt> CellEngineVoxelSimulationSpace::GetRandomParticles(const UnsignedInt NumberOfReactants, map<EntityIdInt, UnsignedInt>& ParticlesKinds)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, ParticlesKinds.size() - 1);

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(std::next(std::begin(ParticlesKinds), static_cast<int>(UniformDistributionObjectUint64t(mt64R)))->first);
            LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
        }
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionForParticle(Particle& ParticleObject)
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMainRandomCondition_Uint64t(0, 1);
        if (UniformDistributionObjectMainRandomCondition_Uint64t(mt64R) == 0)
            return;

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectNumberOfReactants_Uint64t(2, 2);

        UnsignedInt AdditionalBoundFactor = 10;

        map<EntityIdInt, UnsignedInt> ParticlesKindsFoundInParticlesProximity;
        map<UnsignedInt, UniqueIdInt> ParticlesSortedByCapacityFoundInParticlesProximity;

        auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);
        for (UnsignedInt PosX = ParticleObject.XCenter - ParticleKindObject.XSizeDiv2 - AdditionalBoundFactor; PosX < ParticleObject.XCenter + ParticleKindObject.XSizeDiv2 + AdditionalBoundFactor; PosX++)
            for (UnsignedInt PosY = ParticleObject.YCenter - ParticleKindObject.YSizeDiv2 - AdditionalBoundFactor; PosX < ParticleObject.YCenter + ParticleKindObject.YSizeDiv2 + AdditionalBoundFactor; PosY++)
                for (UnsignedInt PosZ = ParticleObject.ZCenter - ParticleKindObject.ZSizeDiv2 - AdditionalBoundFactor; PosX < ParticleObject.ZCenter + ParticleKindObject.ZSizeDiv2 + AdditionalBoundFactor; PosZ++)
                    if (GetSpaceVoxel(PosX, PosZ, PosY) != 0 && GetSpaceVoxel(PosX, PosZ, PosY) != ParticleObject.EntityId)
                    {
                        UniqueIdInt ParticleIndex = GetSpaceVoxel(PosX, PosZ, PosY);
                        ParticlesSortedByCapacityFoundInParticlesProximity[GetParticleFromIndex(ParticleIndex).ListOfVoxels.size()] = ParticleIndex;

                        ParticlesKindsFoundInParticlesProximity[ParticlesKindsManagerObject.GetParticleKind(ParticleIndex).EntityId]++;
                    }

        UnsignedInt NumberOfTries = 0;
        while (NumberOfTries <= 10)
        {
            NumberOfTries++;

            UnsignedInt NumberOfReactants = UniformDistributionObjectNumberOfReactants_Uint64t(mt64R);
            auto RandomParticlesTypes = GetRandomParticles(NumberOfReactants, ParticlesKindsFoundInParticlesProximity);

            sort(begin(RandomParticlesTypes), end(RandomParticlesTypes));
            auto IteratorUnique = unique(begin(RandomParticlesTypes), end(RandomParticlesTypes));
            if (IteratorUnique == RandomParticlesTypes.end())
            {
                vector<string> ParticlesSymbolsForReactionToSort;
                ParticlesSymbolsForReactionToSort.reserve(10);
                for (auto& RandomParticleType : RandomParticlesTypes)
                    ParticlesSymbolsForReactionToSort.emplace_back(ParticlesKindsManagerObject.GetParticleKind(RandomParticleType).Symbol);
                std::sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());

                string ReactionSymbolsStr;
                for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
                    ReactionSymbolsStr += (ParticleSymbolForReaction + " + ");
                LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

                auto ReactionIter = ReactionsIdByString.find(ReactionSymbolsStr);
                if (ReactionIter != ReactionsIdByString.end())
                {
                    Reaction ReactionObject = Reactions[ReactionIter->second];

                    bool IsPossible = all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [&ParticlesKindsFoundInParticlesProximity](const ParticleKindForReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsFoundInParticlesProximity[ReactionReactant.EntityId]; });
                    if (IsPossible == true)
                    {
                        vector<UniqueIdInt> ParticlesChosenForReaction;

                        for (const auto& ParticleObjectIndexToBeErased : ParticlesSortedByCapacityFoundInParticlesProximity)
                        {
                            auto& ParticleObjectEntityId = GetParticleFromIndex(ParticleObjectIndexToBeErased.second).EntityId;

                            //CO JESLI WIELE ROZNYCH ODCINKOW DNA TO MA TO SAMO ENTITY_ID - to mozliwe tylko dla DNA
                            auto ReactantIterator = find_if(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [&ParticleObjectEntityId](ParticleKindForReaction& ParticleKindObject){ return ParticleKindObject.EntityId == ParticleObjectEntityId; });
                            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantIterator->Counter > 0 && ReactantIterator->ToRemoveInReaction == true)
                            {
                                if (CellEngineUseful::IsDNAorRNA(ParticleObjectEntityId) == true)
                                {
                                    auto& ParticleObjectForReaction = GetParticleFromIndex(ParticleObjectIndexToBeErased.second);

                                    bool FoundSequenceOuter = false;
                                    for (UnsignedInt GenomeIndex = ParticleObjectForReaction.GenomeIndex - 20; GenomeIndex < ParticleObjectForReaction.GenomeIndex + 20; GenomeIndex++)
                                    {
                                        bool FoundSequenceInner = true;
                                        auto ParticleKindForReactionObject = ParticlesKindsManagerObject.GetParticleKind(ReactantIterator->EntityId);
                                        for (UnsignedInt NucleotideNum = 0; NucleotideNum < ParticleKindForReactionObject.Sequence.size(); NucleotideNum++)
                                            if (ParticleKindForReactionObject.Sequence[NucleotideNum] != GetParticleFromIndex(Genome1[GenomeIndex + NucleotideNum]).ChainId)
                                            {
                                                FoundSequenceInner = false;
                                                break;
                                            }
                                        if (FoundSequenceInner == true)
                                        {
                                            FoundSequenceOuter = true;
                                            break;
                                        }
                                    }
                                    if (FoundSequenceOuter == true)
                                        ParticlesChosenForReaction.emplace_back(ParticleObjectIndexToBeErased.second);
                                }
                                else
                                    ParticlesChosenForReaction.emplace_back(ParticleObjectIndexToBeErased.second);
                            }

                            ReactantIterator->Counter--;
                        }

                        vector<vector3_16> Centers;

                        for (const auto& ParticleChosenForReaction : ParticlesChosenForReaction)
                        {
                            auto& ParticleObjectToBeErased = GetParticleFromIndex(ParticleChosenForReaction);

                            Centers.emplace_back(ParticleObjectToBeErased.XCenter, ParticleObjectToBeErased.YCenter, ParticleObjectToBeErased.ZCenter);

                            for (auto& VoxelCoordinates : ParticleObjectToBeErased.ListOfVoxels)
                                GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z) = GetZeroSimulationSpaceVoxel();

                            FreeIndexesOfParticles.push(ParticleChosenForReaction);
                            Particles.erase(ParticleChosenForReaction);
                        }

                        UnsignedInt CenterIndex = 0;
                        sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForReaction& PK1, ParticleKindForReaction& PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); } );
                        for (const auto& ReactionProduct : ReactionObject.Products)
                        {
                            UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ReactionProduct.EntityId, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));

                            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

                            for (auto& ParticleKindVoxel : ParticleKindObjectForProduct.ListOfVoxels)
                                if (CenterIndex < Centers.size())
                                {
                                    vector3_64 NewVoxel(Centers[CenterIndex].X - ParticleKindObjectForProduct.XSizeDiv2 + ParticleKindVoxel.X, Centers[CenterIndex].Y - ParticleKindObjectForProduct.YSizeDiv2 + ParticleKindVoxel.Y, Centers[CenterIndex].Z - ParticleKindObjectForProduct.ZSizeDiv2 + ParticleKindVoxel.Z);
                                    GetSpaceVoxel(NewVoxel.X, NewVoxel.Y, NewVoxel.Z) = ParticleIndex;
                                    GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(NewVoxel.X, NewVoxel.Y, NewVoxel.Z);
                                }
                                else
                                {
                                }

                            CenterIndex++;
                        }

                        break;
                    }
                    else
                        LoggersManagerObject.Log(STREAM("Particles types are the same!"));
                }
                else
                    LoggersManagerObject.Log(STREAM("Reaction for particles does not exist!"));
            }
            else
                LoggersManagerObject.Log(STREAM("Particles types for reaction are not unique!"));
        }
    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
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

void CellEngineVoxelSimulationSpace::GenerateRandomReactionsInWholeVoxelSimulationSpace()
{
    try
    {
        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                    if (GetSpaceVoxel(PosX, PosZ, PosY) != 0)
                    {
                        auto& ParticleObject = GetParticleFromIndex(GetSpaceVoxel(PosX, PosZ, PosY));
                        if (ParticleObject.SelectedForReaction == false)
                            GenerateRandomReactionForParticle(ParticleObject);
                    }
    }
    CATCH("generation random reactions in whole voxel simulation space")
}

void CellEngineVoxelSimulationSpace::EraseAllDNAParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
            {
                for (auto& VoxelCoordinates : ParticleObject.second.ListOfVoxels)
                    GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z) = GetZeroSimulationSpaceVoxel();

                FreeIndexesOfParticles.push(ParticleObject.first);
            }

        const auto RemovedDNAParticlesCounter = erase_if(Particles, [](const pair<UniqueIdInt, Particle>& item) { auto const& [key, value] = item; return (value.EntityId == CellEngineConfigDataObject.DNAIdentifier); });
        LoggersManagerObject.Log(STREAM("RemovedDNAParticlesCounter = " << RemovedDNAParticlesCounter));

        UnsignedInt DNAParticleCounter = 0;
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
                DNAParticleCounter++;
        LoggersManagerObject.Log(STREAM("DNAParticleCounter = " << DNAParticleCounter));

        Genome1.clear();
        Genome2.clear();
    }
    CATCH("erasing all dna particles")
}

void CellEngineVoxelSimulationSpace::UpdateRandomPositions(const UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, const UnsignedInt Size)
{
    switch (RandomMoveDirection)
    {
        case 1 : RandomPosX += Size; break;
        case 2 : RandomPosX -= Size; break;
        case 3 : RandomPosY += Size; break;
        case 4 : RandomPosY -= Size; break;
        case 5 : RandomPosZ += Size; break;
        case 6 : RandomPosZ -= Size; break;
        default : break;
    }
}

void CellEngineVoxelSimulationSpace::GenerateParticle(const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, vector<UniqueIdInt>& Genome, const vector3_16 UniqueColorParam)
{
    try
    {
        UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), EntityId, ChainId, GenomeIndex, UniqueColorParam));

        for (UnsignedInt PosX = StartPosX; PosX < StartPosX + ParticleSizeX; PosX++)
            for (UnsignedInt PosY = StartPosY; PosY < StartPosY + ParticleSizeY; PosY++)
                for (UnsignedInt PosZ = StartPosZ; PosZ < StartPosZ + ParticleSizeZ; PosZ++)
                {
                    GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
                    GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
                }

        Genome.emplace_back(ParticleIndex);
    }
    CATCH("generating particle")
}

void CellEngineVoxelSimulationSpace::GenerateTwoPairedNucleotides(const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt AddSizeX, const UnsignedInt AddSizeY, const UnsignedInt AddSizeZ, const vector3_16 UniqueColorParam)
{
    try
    {
        GenerateParticle(CellEngineConfigDataObject.DNAIdentifier, ChainId, Genome1.size(), StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, Genome1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
        GenerateParticle(CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainId(ChainId), Genome2.size(), StartPosX + AddSizeX, StartPosY + AddSizeY, StartPosZ + AddSizeZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, Genome2, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
    }
    CATCH("generating two paired nucleotides")
}

bool CellEngineVoxelSimulationSpace::TestFormerForbiddenPositions(unordered_set<string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, const UnsignedInt Size)
{
    UnsignedInt PosX = RandomPosX;
    UnsignedInt PosY = RandomPosY;
    UnsignedInt PosZ = RandomPosZ;

    UpdateRandomPositions(RandomMoveDirection, PosX, PosY, PosZ, Size);

    return (TestedFormerForbiddenPositions.find(to_string(PosX) + "|" + to_string(PosY) + "|" + to_string(PosZ)) != TestedFormerForbiddenPositions.end());
}

tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineVoxelSimulationSpace::EraseLastRandomParticle(vector<UniqueIdInt>& Genome)
{
    UnsignedInt LocalRandomPosX, LocalRandomPosY, LocalRandomPosZ;

    try
    {
        UnsignedInt PreviousParticleIndex = Genome.back();
        Genome.pop_back();

        LocalRandomPosX = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].X;
        LocalRandomPosY = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Y;
        LocalRandomPosZ = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Z;

        for (auto& VoxelCoordinates : GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels)
            GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z) = GetZeroSimulationSpaceVoxel();

        Particles.erase(PreviousParticleIndex);

        FreeIndexesOfParticles.push(PreviousParticleIndex);

        LoggersManagerObject.Log(STREAM("ERASED PARTICLE PreviousParticleIndex = " << PreviousParticleIndex << " RandomPosX = " << LocalRandomPosX << " RandomPosY = " << LocalRandomPosY << " RandomPosZ = " << LocalRandomPosZ));
    }
    CATCH("erasing last random particle")

    return { LocalRandomPosX, LocalRandomPosY, LocalRandomPosZ };
}

UnsignedInt Sqr(UnsignedInt Value)
{
    return Value * Value;
}

void CellEngineVoxelSimulationSpace::GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5)
{
    try
    {
        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        LoggersManagerObject.InitializePrintingParameters(false, false, false, false, false, false, false, false, false, false, false, false, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);

        UnsignedInt NumberOfGeneratedNucleotides = 0;

        unordered_set<string> TestedFormerForbiddenPositions;

        vector<UnsignedInt> RandomMovesDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = [](vector<UnsignedInt>& RandomMovesDirections) { return all_of(RandomMovesDirections.cbegin(), RandomMovesDirections.cend(), [] (const UnsignedInt Element) { return Element == 1; }); };

        auto timeStart = clock();

        while (NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
        {
            if ((clock() - timeStart) / CLOCKS_PER_SEC >= 15)
                break;

            UnsignedInt RandomMoveDirection = 0;

            do
            {
                RandomMoveDirection = UniformDistributionObjectMoveOfParticle_Uint64t(mt64R);

                if (TestFormerForbiddenPositions(TestedFormerForbiddenPositions, RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize1) == true)
                    RandomMovesDirections[RandomMoveDirection - 1] = 1;

                LoggersManagerObject.Log(STREAM("RandomMoveDirection = " << RandomMoveDirection << " " << RandomMovesDirections[RandomMoveDirection - 1] << " " << CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections)));
            }
            while (RandomMovesDirections[RandomMoveDirection - 1] == 1 && CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections) == false);

            if (RandomMovesDirections[RandomMoveDirection - 1] == 0)
            {
                UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize2);

                RandomMovesDirections[RandomMoveDirection - 1] = 1;
            }

            bool EmptyVoxelSpaceForNewNucleotideBool = true;
            while (EmptyVoxelSpaceForNewNucleotideBool == true && NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
            {
                for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + ParticleSizeX; PosX++)
                    for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + ParticleSizeY; PosY++)
                        for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + ParticleSizeZ; PosZ++)
                            if (GetSpaceVoxel(PosX, PosY, PosZ) != 0)
                            {
                                LoggersManagerObject.Log(STREAM("NOT EMPTY A POS = " << PosX << " " << PosY << " " << PosZ << " " << GetSpaceVoxel(PosX, PosY, PosZ)));

                                EmptyVoxelSpaceForNewNucleotideBool = false;

                                goto BreakOutOfLoop;
                            }

                BreakOutOfLoop:;

                if (EmptyVoxelSpaceForNewNucleotideBool == false)
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize3);

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
                {
                    EmptyVoxelSpaceForNewNucleotideBool = false;
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize4);
                }

                LoggersManagerObject.Log(STREAM("EmptyVoxelSpaceForNewNucleotideBool = " << EmptyVoxelSpaceForNewNucleotideBool << " " << CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections)));

                if (EmptyVoxelSpaceForNewNucleotideBool == true)
                {
                    NumberOfGeneratedNucleotides++;

                    ChainIdInt ChainId = UniformDistributionObjectChainOfParticle_Uint64t(mt64R);

                    if (RandomMoveDirection == 1 || RandomMoveDirection == 2)
                        GenerateTwoPairedNucleotides(CellEngineConfigDataObject.DNAIdentifier, ChainId, Genome1.size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, 1, ParticleSizeZ, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
                    else
                    if (RandomMoveDirection == 3 || RandomMoveDirection == 4 || RandomMoveDirection == 5 || RandomMoveDirection == 6)
                        GenerateTwoPairedNucleotides(CellEngineConfigDataObject.DNAIdentifier, ChainId, Genome1.size(), RandomPosX, RandomPosY, RandomPosZ, 1, ParticleSizeY, ParticleSizeZ, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

                    fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);

                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize5);

                    LoggersManagerObject.Log(STREAM("ADDED PARTICLE NumberOfGeneratedNucleotides = " << NumberOfGeneratedNucleotides << " RX = " << RandomPosX << " RY = " << RandomPosY << " RZ = " << RandomPosZ));
                }
            }

            if (EmptyVoxelSpaceForNewNucleotideBool == false && CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections) == true)
            {
                NumberOfGeneratedNucleotides--;

                tie(RandomPosX, RandomPosY, RandomPosZ) = EraseLastRandomParticle(Genome1);

                TestedFormerForbiddenPositions.insert(to_string(RandomPosX) + "|" + to_string(RandomPosY) + "|" + to_string(RandomPosZ));

                EraseLastRandomParticle(Genome2);

                fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);
            }

            LoggersManagerObject.Log(STREAM("END OF GOING IN ONE DIRECTION"));
        }

        GetMinMaxCoordinatesForDNA();

        LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}

void CellEngineVoxelSimulationSpace::GetMinMaxCoordinatesForDNA()
{
    try
    {
        for (auto& LocalNewParticleIndex : Genome1)
            GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
        for (auto& LocalNewParticleIndex : Genome2)
            GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
    }
    CATCH("getting min max of coordinates for DNA")
}

void CellEngineVoxelSimulationSpace::SaveGenomeDataToFile(UnsignedInt ParticleSize)
{
    try
    {
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_POSITIONS.DAT");
        ofstream FileToWriteGenome;
        FileToWriteGenome.open(FileName, ios_base::out | ios_base::trunc);
        FileToWriteGenome << to_string(ParticleSize) << endl;
        FileToWriteGenome << to_string(Genome1.size()) << endl;
        for (const auto& Nucleotide : Genome1)
        {
            EntityIdInt EntityId = GetParticleFromIndex(Nucleotide).EntityId;
            ChainIdInt ChainId = GetParticleFromIndex(Nucleotide).ChainId;
            UniqueIdInt GenomeIndex = GetParticleFromIndex(Nucleotide).GenomeIndex;
            UnsignedInt PosX = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].X;
            UnsignedInt PosY = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Y;
            UnsignedInt PosZ = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Z;
            FileToWriteGenome << to_string(EntityId) << "," << to_string(ChainId)  << "," << to_string(GenomeIndex) << "," << to_string(PosX) << "," << to_string(PosY) << "," << to_string(PosZ) << endl;
        }
        FileToWriteGenome.close();
    }
    CATCH("saving genome data to file")
}

void CellEngineVoxelSimulationSpace::ReadGenomeDataFromFile(bool Paired)
{
    try
    {
        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        string Line, Word;
        vector<string> Row;

        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_POSITIONS.DAT");
        fstream FileToReadGenome(FileName, ios::in);

        getline(FileToReadGenome, Line);
        UnsignedInt ParticleSize = stoi(Line);
        getline(FileToReadGenome, Line);
        UnsignedInt GenomeSize = stoi(Line);

        UnsignedInt PrevStartPosX = 0;
        UnsignedInt GenomeIndex = 0;

        while(getline(FileToReadGenome, Line))
        {
            Row.clear();
            stringstream Str(Line);

            while(getline(Str, Word, ','))
                Row.push_back(Word);

            UnsignedInt StartPosX = stoi(Row[3]);
            UnsignedInt StartPosY = stoi(Row[4]);
            UnsignedInt StartPosZ = stoi(Row[5]);

            if (Paired == true)
            {
                if (abs(static_cast<long>(PrevStartPosX - StartPosX)) > 0)
                    GenerateTwoPairedNucleotides(stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, ParticleSize, 1, ParticleSize, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
                else
                    GenerateTwoPairedNucleotides(stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, 1, ParticleSize, ParticleSize, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
            }
            else
                GenerateParticle(stoi(Row[0]), stoi(Row[1]), stoi(Row[2]), stoi(Row[3]), stoi(Row[4]), stoi(Row[5]), ParticleSize, ParticleSize, ParticleSize, Genome1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

            PrevStartPosX = StartPosX;

            GenomeIndex++;
        }

        FileToReadGenome.close();

        GetMinMaxCoordinatesForDNA();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("reading genome data from file")
}

void CellEngineVoxelSimulationSpace::ReadGenomeSequenceFromFile()
{
    try
    {
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_SEQUENCE.DAT");

        fstream FileToReadGenome(FileName, ios::in);

        getline(FileToReadGenome, GenomeLine);

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genome1.size(); GenomeIndex++)
        {
            GetParticleFromIndex(Genome1[GenomeIndex]).ChainId = CellEngineUseful::GetChainIdFromLetter(GenomeLine[GenomeIndex]);
            GetParticleFromIndex(Genome2[GenomeIndex]).ChainId = CellEngineUseful::GetPairedChainId(CellEngineUseful::GetChainIdFromLetter(GenomeLine[GenomeIndex]));
        }
    }
    CATCH("reading real genome data from file")
}

void CellEngineVoxelSimulationSpace::TestGeneratedGenomeCorrectness(const UnsignedInt ParticleSize)
{
    try
    {
        bool FoundBreak = false;

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genome1.size(); GenomeIndex++)
            if (GenomeIndex > 3)
            {
                if ((abs(GetParticleFromIndex(Genome1[GenomeIndex]).ListOfVoxels[0].X - GetParticleFromIndex(Genome1[GenomeIndex - 1]).ListOfVoxels[0].X) != ParticleSize) &&
                    (abs(GetParticleFromIndex(Genome1[GenomeIndex]).ListOfVoxels[0].Y - GetParticleFromIndex(Genome1[GenomeIndex - 1]).ListOfVoxels[0].Y) != ParticleSize) &&
                    (abs(GetParticleFromIndex(Genome1[GenomeIndex]).ListOfVoxels[0].Z - GetParticleFromIndex(Genome1[GenomeIndex - 1]).ListOfVoxels[0].Z) != ParticleSize))
                {
                    LoggersManagerObject.Log(STREAM("GenomeIndex = " << GenomeIndex << " ChainId = " << GetParticleFromIndex(Genome1[GenomeIndex]).ChainId << " Letter = " << CellEngineUseful::GetLetterForDNAChainId(GetParticleFromIndex(Genome1[GenomeIndex]).ChainId)));
                    LoggersManagerObject.Log(STREAM("DIFF X = " << GetParticleFromIndex(Genome1[GenomeIndex]).ListOfVoxels[0].X << " " << GetParticleFromIndex(Genome1[GenomeIndex - ParticleSize]).ListOfVoxels[0].X));
                    LoggersManagerObject.Log(STREAM("DIFF Y = " << GetParticleFromIndex(Genome1[GenomeIndex]).ListOfVoxels[0].Y << " " << GetParticleFromIndex(Genome1[GenomeIndex - ParticleSize]).ListOfVoxels[0].Y));
                    LoggersManagerObject.Log(STREAM("DIFF Z = " << GetParticleFromIndex(Genome1[GenomeIndex]).ListOfVoxels[0].Z << " " << GetParticleFromIndex(Genome1[GenomeIndex - ParticleSize]).ListOfVoxels[0].Z));
                    FoundBreak = true;
                    break;
                }
            }

        if (FoundBreak == false)
            LoggersManagerObject.Log(STREAM("Genome is continuous and correctly generated - OK!"));
    }
    CATCH("testing generated genome correctness")
}