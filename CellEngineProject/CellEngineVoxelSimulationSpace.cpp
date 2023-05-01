
#include "FileUtils.h"

#include "CellEngineAtom.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"

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

inline SimulationSpaceVoxel GetZerSimulationSpaceVoxel()
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
                    GetSpaceVoxel(PosX, PosZ, PosY) = GetZerSimulationSpaceVoxel();
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

void CellEngineVoxelSimulationSpace::SetStartValuesForSpaceMinMax()
{
    XMin = YMin = ZMin = 10000;
    XMax = YMax = ZMax = 0;
}

void CellEngineVoxelSimulationSpace::GetMinMaxOfCoordinates(const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ)
{
    try
    {
        XMin = min(PosX, XMin);
        XMax = max(PosX, XMax);
        YMin = min(PosY, YMin);
        YMax = max(PosY, YMax);
        ZMin = min(PosZ, ZMin);
        ZMax = max(PosZ, ZMax);
    }
    CATCH("getting min max of coordinates")
}

[[nodiscard]] std::stringstream CellEngineVoxelSimulationSpace::PrintSpaceMinMaxValues() const
{
    std::stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << std::to_string(XMin) << " ][ Xmax = " << std::to_string(XMax) << " ][ Ymin = " << std::to_string(YMin) << " ][ Ymax = " << std::to_string(YMax) << " ][ Zmin = " << std::to_string(ZMin) << " ][ Zmax = " << std::to_string(XMax) << " ] " << std::endl;
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

        GetMinMaxOfCoordinates(PosX, PosY, PosZ);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
        {
            GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
        }
    }
    CATCH("setting atom in voxel simulation space")
}

void CellEngineVoxelSimulationSpace::AddParticleKind(const ParticleKind& ParticleParam)
{
    ParticlesKinds.emplace_back(ParticleParam);
}

UniqueIdInt CellEngineVoxelSimulationSpace::AddNewParticle(UniqueIdInt ParticleIndex, const Particle& ParticleParam)
{
    #ifdef PARTICLES_IN_VECTOR
    Particles.emplace_back(ParticleParam);
    return MaxParticleIndex = Particles.size() - 1;
    #else
    Particles[ParticleIndex] = ParticleParam;
    return MaxParticleIndex = ParticleIndex;
    #endif
}

void CellEngineVoxelSimulationSpace::AddReaction(const Reaction& ReactionParam)
{
    try
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsIdByString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
    CATCH("adding reaction")
}

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        vector<UnsignedInt> LocalNewParticlesIndexes;

        LoggersManagerObject.Log(STREAM("PARTICLES RANGE = " << MaxParticleIndex << " " << (MaxParticleIndex + NumberOfRandomParticles) << " " << NumberOfRandomParticles));

        UnsignedInt LocalMaxParticleIndex = MaxParticleIndex;
        for (UniqueIdInt ParticleIndex = LocalMaxParticleIndex + 1; ParticleIndex <= LocalMaxParticleIndex + NumberOfRandomParticles; ParticleIndex++)
            LocalNewParticlesIndexes.emplace_back(AddNewParticle(ParticleIndex, Particle(CellEngineConfigDataObject.DNAIdentifier, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));

        vector<vector3_64> FilledVoxelsForRandomParticle;

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectObjectOfParticle_Uint64t(1, NumberOfRandomParticles);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    GetSpaceVoxel(PosX, PosY, PosZ) = GetZerSimulationSpaceVoxel();

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
                                    GetSpaceVoxel(VoxelForRandomParticle.X, VoxelForRandomParticle.Y, VoxelForRandomParticle.Z) = GetZerSimulationSpaceVoxel();

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
    }
    CATCH("generating random particles in selected space")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusion(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        std::uniform_int_distribution<SignedInt> UniformDistributionObjectSizeOfParticle_int64t(-1, 1);

        vector<vector3_64> NewVoxelsForParticle;

        if (EndParticleIndexParam == 0)
        {
            StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
            EndParticleIndexParam = MaxParticleIndex;
        }

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto& LocalParticleObject = Particles[ParticleIndex];

            for (auto& VoxelForParticle : LocalParticleObject.ListOfVoxels)
                GetSpaceVoxel(VoxelForParticle.X, VoxelForParticle.Y, VoxelForParticle.Z) = GetZerSimulationSpaceVoxel();

            SignedInt ShiftX = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftY = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftZ = UniformDistributionObjectSizeOfParticle_int64t(mt64R);

            NewVoxelsForParticle.clear();
            bool Collision = false;

            for (auto& VoxelForParticle : LocalParticleObject.ListOfVoxels)
                if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) == 0 && VoxelForParticle.X + ShiftX >= StartXPosParam && VoxelForParticle.X + ShiftX < StartXPosParam + SizeXParam && VoxelForParticle.Y + ShiftY >= StartYPosParam && VoxelForParticle.Y + ShiftY < StartYPosParam + SizeYParam && VoxelForParticle.Z + ShiftZ >= StartZPosParam && VoxelForParticle.Z + ShiftZ < StartZPosParam + SizeZParam)
                {
                    GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) = ParticleIndex;
                    NewVoxelsForParticle.emplace_back(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ);
                }
                else
                {
                    for (auto& NewVoxelForParticle : NewVoxelsForParticle)
                        GetSpaceVoxel(NewVoxelForParticle.X, NewVoxelForParticle.Y, NewVoxelForParticle.Z) = GetZerSimulationSpaceVoxel();

                    for (auto& OldVoxelForParticle : LocalParticleObject.ListOfVoxels)
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

void CellEngineVoxelSimulationSpace::GetDNASequenceFromNucleotides()
{
    try
    {
//        vector<char> Genome;
//        string GenomeStr0;
//        string GenomeStr1;
//
//        UnsignedInt NumberOfNucleotides = 0;
//
//        Particle* PrevParticlePtr = nullptr;
//
//        LoggersManagerObject.Log(STREAM("LIST SIZE " << ParticlesKinds[1].ParticlesObjects.size() << " " << ParticlesKinds[2].ParticlesObjects.size() << " " << ParticlesKinds[3].ParticlesObjects.size() << " " << ParticlesKinds[4].ParticlesObjects.size()));
//        auto& ParticleObject = ParticlesKinds[1].ParticlesObjects.front();
//
//        for (UnsignedInt Counter = 1; Counter <= 4; Counter++)
//        {
//            UnsignedInt FoundEmpty = 0;
//            for (const auto& ParticleObjectLoop : ParticlesKinds[Counter].ParticlesObjects)
//                if (ParticleObjectLoop.ListOfVoxels.empty() == true)
//                    FoundEmpty++;
//            LoggersManagerObject.Log(STREAM("EMPTY " << Counter << " " << FoundEmpty));
//        }
//
//        while (NumberOfNucleotides < 580079)
//        {
//            unordered_map<Particle* , UnsignedInt> NeighbourNucleotideVoxelCounter;
//
//            LoggersManagerObject.Log(STREAM("1")); getchar();
//
//            if (ParticleObject.ListOfVoxels.empty() == true)
//                LoggersManagerObject.Log(STREAM("EMPTY VECTOR"));
//
//            bool WRITE = false;
//
//            for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
//            {
//                auto Voxel = GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z);
//
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA W " << VoxelCoordinates.X << " " << VoxelCoordinates.Y << " " << VoxelCoordinates.Z << " ParticlePtr " << Voxel.ParticlePtr << " EntityId = " << Voxel.EntityId << " ChainId = " << Voxel.ChainId));
//
//                auto VoxelNeighbourX1P = GetSpaceVoxel(VoxelCoordinates.X + 1, VoxelCoordinates.Y, VoxelCoordinates.Z);
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA X1P " << VoxelCoordinates.X + 1 << " " << VoxelCoordinates.Y << " " << VoxelCoordinates.Z << " ParticlePtr " << VoxelNeighbourX1P.ParticlePtr << " EntityId = " << VoxelNeighbourX1P.EntityId << " ChainId = " << VoxelNeighbourX1P.ChainId));
//                if (VoxelNeighbourX1P.EntityId != 0 && CellEngineUseful::IsDNAorRNA(VoxelNeighbourX1P.EntityId) == true && VoxelNeighbourX1P.ParticlePtr != PrevParticlePtr && VoxelNeighbourX1P.ParticlePtr != Voxel.ParticlePtr)
//                {
//                    if (WRITE) LoggersManagerObject.Log(STREAM("XP"));
//                    NeighbourNucleotideVoxelCounter[VoxelNeighbourX1P.ParticlePtr]++;
//                }
//
//                auto VoxelNeighbourX1M = GetSpaceVoxel(VoxelCoordinates.X - 1, VoxelCoordinates.Y, VoxelCoordinates.Z);
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA X1M " << VoxelCoordinates.X - 1 << " " << VoxelCoordinates.Y << " " << VoxelCoordinates.Z << " ParticlePtr " << VoxelNeighbourX1M.ParticlePtr << " EntityId = " << VoxelNeighbourX1M.EntityId << " ChainId = " << VoxelNeighbourX1M.ChainId));
//                if (VoxelNeighbourX1M.EntityId != 0 && CellEngineUseful::IsDNAorRNA(VoxelNeighbourX1M.EntityId) == true && VoxelNeighbourX1M.ParticlePtr != PrevParticlePtr && VoxelNeighbourX1M.ParticlePtr != Voxel.ParticlePtr)
//                {
//                    if (WRITE) LoggersManagerObject.Log(STREAM("XM"));
//                    NeighbourNucleotideVoxelCounter[VoxelNeighbourX1M.ParticlePtr]++;
//                }
//
//                auto VoxelNeighbourY1P = GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y + 1, VoxelCoordinates.Z);
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA Y1P " << VoxelCoordinates.X << " " << VoxelCoordinates.Y + 1 << " " << VoxelCoordinates.Z << " ParticlePtr " << VoxelNeighbourY1P.ParticlePtr << " EntityId = " << VoxelNeighbourY1P.EntityId << " ChainId = " << VoxelNeighbourY1P.ChainId));
//                if (VoxelNeighbourY1P.EntityId != 0 && CellEngineUseful::IsDNAorRNA(VoxelNeighbourY1P.EntityId) == true && VoxelNeighbourY1P.ParticlePtr != PrevParticlePtr && VoxelNeighbourY1P.ParticlePtr != Voxel.ParticlePtr)
//                {
//                    if (WRITE) LoggersManagerObject.Log(STREAM("YP"));
//                    NeighbourNucleotideVoxelCounter[VoxelNeighbourY1P.ParticlePtr]++;
//                }
//
//                auto VoxelNeighbourY1M = GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y - 1, VoxelCoordinates.Z);
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA Y1M " << VoxelCoordinates.X << " " << VoxelCoordinates.Y - 1 << " " << VoxelCoordinates.Z << " ParticlePtr " << VoxelNeighbourY1M.ParticlePtr << " EntityId = " << VoxelNeighbourY1M.EntityId << " ChainId = " << VoxelNeighbourY1M.ChainId));
//                if (VoxelNeighbourY1M.EntityId != 0 && CellEngineUseful::IsDNAorRNA(VoxelNeighbourY1M.EntityId) == true && VoxelNeighbourY1M.ParticlePtr != PrevParticlePtr && VoxelNeighbourY1M.ParticlePtr != Voxel.ParticlePtr)
//                {
//                    if (WRITE) LoggersManagerObject.Log(STREAM("YM"));
//                    NeighbourNucleotideVoxelCounter[VoxelNeighbourY1M.ParticlePtr]++;
//                }
//
//                auto VoxelNeighbourZ1P = GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z + 1);
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA Z1P " << VoxelCoordinates.X << " " << VoxelCoordinates.Y << " " << VoxelCoordinates.Z + 1 << " ParticlePtr " << VoxelNeighbourZ1P.ParticlePtr << " EntityId = " << VoxelNeighbourZ1P.EntityId << " ChainId = " << VoxelNeighbourZ1P.ChainId));
//                if (VoxelNeighbourZ1P.EntityId != 0 && CellEngineUseful::IsDNAorRNA(VoxelNeighbourZ1P.EntityId) == true && VoxelNeighbourZ1P.ParticlePtr != PrevParticlePtr && VoxelNeighbourZ1P.ParticlePtr != Voxel.ParticlePtr)
//                {
//                    if (WRITE) LoggersManagerObject.Log(STREAM("ZP"));
//                    NeighbourNucleotideVoxelCounter[VoxelNeighbourZ1P.ParticlePtr]++;
//                }
//
//                auto VoxelNeighbourZ1M = GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z - 1);
//                if (WRITE) LoggersManagerObject.Log(STREAM("VOXEL DATA Z1P " << VoxelCoordinates.X << " " << VoxelCoordinates.Y << " " << VoxelCoordinates.Z - 1 << " ParticlePtr " << VoxelNeighbourZ1P.ParticlePtr << " EntityId = " << VoxelNeighbourZ1P.EntityId << " ChainId = " << VoxelNeighbourZ1P.ChainId));
//                if (VoxelNeighbourZ1M.EntityId != 0 && CellEngineUseful::IsDNAorRNA(VoxelNeighbourZ1M.EntityId) == true && VoxelNeighbourZ1M.ParticlePtr != PrevParticlePtr && VoxelNeighbourZ1M.ParticlePtr != Voxel.ParticlePtr)
//                {
//                    if (WRITE) LoggersManagerObject.Log(STREAM("ZM"));
//                    NeighbourNucleotideVoxelCounter[VoxelNeighbourZ1M.ParticlePtr]++;
//                }
//            }
//
//            LoggersManagerObject.Log(STREAM("2")); getchar();
//
//            auto NeighbourNucleotideWithMaximalNumberOfBorderVoxel = std::max_element(NeighbourNucleotideVoxelCounter.begin(), NeighbourNucleotideVoxelCounter.end(), [](const auto &x, const auto &y) { return x.second < y.second; });
//
//            if (NeighbourNucleotideVoxelCounter.empty() == false)
//                for (const auto& NeighbourNucleotideVoxelCounterOneObject : NeighbourNucleotideVoxelCounter)
//                    LoggersManagerObject.Log(STREAM("MAP = " << NeighbourNucleotideVoxelCounterOneObject.first << " " << NeighbourNucleotideVoxelCounterOneObject.second));
//            else
//                LoggersManagerObject.Log(STREAM("EMPTY MAP"));
//
//            LoggersManagerObject.Log(STREAM("3")); getchar();
//
//            PrevParticlePtr = &ParticleObject;
//
//            LoggersManagerObject.Log(STREAM("4")); getchar();
//
//            ParticleObject = *NeighbourNucleotideWithMaximalNumberOfBorderVoxel->first;
//
//            LoggersManagerObject.Log(STREAM("5")); getchar();
//
//            NumberOfNucleotides++;
//
//            Genome.emplace_back(ParticleObject.ChainId);
//
//            LoggersManagerObject.Log(STREAM("6")); getchar();
//
//            GenomeStr0 += CellEngineUseful::GetLetterForDNAChainId0(ParticleObject.ChainId);
//            GenomeStr1 += CellEngineUseful::GetLetterForDNAChainId1(ParticleObject.ChainId);
//
//            LoggersManagerObject.Log(STREAM("NumberOfNucleotides = " << NumberOfNucleotides << " PrevParticlePtr = " << PrevParticlePtr << " ParticleObjectPtr = " << &ParticleObject << " PrevParticlePtr->UniqueId = " << PrevParticlePtr->UniqueId << " ParticleObject.UniqueId = " << ParticleObject.UniqueId << " ParticleObject.ChainId " << ParticleObject.ChainId));
//            getchar();
//        }
//
//        FileUtils::RewriteTextToFile("Genome0.txt", GenomeStr0);
//        FileUtils::RewriteTextToFile("Genome1.txt", GenomeStr1);
    }
    CATCH("getting dna sequence from nucleotides")
}
