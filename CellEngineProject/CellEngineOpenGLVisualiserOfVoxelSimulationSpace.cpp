
#include <omp.h>
#include "../Common/Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include <string>

#include <sb7.h>
#include <tuple>

#include "CellEngineDataFile.h"

#include "CellEngineUseful.h"
#include "CellEngineConstants.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

constexpr bool DrawError = false;

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetStartPositions()
{
    return { SelectionStartXPos, SelectionStartYPos, SelectionStartZPos };
}

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetSteps()
{
    return { SelectionStepX, SelectionStepY, SelectionStepZ };
}

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetSizes()
{
    return { SelectionSizeX, SelectionSizeY, SelectionSizeZ };
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SetVoxelSpaceSelection(const UnsignedInt SelectionStartXParam, const UnsignedInt SelectionStartYParam, const UnsignedInt SelectionStartZParam, const UnsignedInt SelectionStepXParam, const UnsignedInt SelectionStepYParam, const UnsignedInt SelectionStepZParam, const UnsignedInt SelectionSizeXParam, const UnsignedInt SelectionSizeYParam, const UnsignedInt SelectionSizeZParam)
{
    SelectionStartXPos = SelectionStartXParam, SelectionStartYPos = SelectionStartYParam, SelectionStartZPos = SelectionStartZParam;
    SelectionStepX = SelectionStepXParam, SelectionStepY = SelectionStepYParam, SelectionStepZ = SelectionStepZParam;
    SelectionSizeX = SelectionSizeXParam, SelectionSizeY = SelectionSizeYParam, SelectionSizeZ = SelectionSizeZParam;
}

CellEngineOpenGLVisualiserOfVoxelSimulationSpace::CellEngineOpenGLVisualiserOfVoxelSimulationSpace()
{
    SetVoxelSpaceSelection(CellEngineConfigDataObject.SimulationSpaceSelectionStartXPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.SimulationSpaceSelectionStartZPos, CellEngineConfigDataObject.SimulationSpaceSelectionStepX, CellEngineConfigDataObject.SimulationSpaceSelectionStepY, CellEngineConfigDataObject.SimulationSpaceSelectionStepZ, CellEngineConfigDataObject.SimulationSpaceSelectionSizeX, CellEngineConfigDataObject.SimulationSpaceSelectionSizeY, CellEngineConfigDataObject.SimulationSpaceSelectionSizeZ);
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SaveVoxelPositionChosenByMouse()
{
    SelectionStartXPos = SaveXMousePosition;
    SelectionStartYPos = SaveYMousePosition;
    SelectionStartZPos = SaveZMousePosition;
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SetSaveXYZPositions(const UnsignedInt SaveXParam, const UnsignedInt SaveYParam, const UnsignedInt SaveZParam)
{
    SaveXMousePosition = SaveXParam;
    SaveYMousePosition = SaveYParam;
    SaveZMousePosition = SaveZParam;
}

static inline float CovertToGraphicsCoordinateSelected(const UnsignedInt StartParam, const UnsignedInt SpaceParam, const UnsignedInt SizeParam)
{
    return ((static_cast<float>((StartParam + SizeParam) - SpaceParam) - static_cast<float>(SizeParam) / 2) * 4);
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::ConvertAtomPosToGraphicCoordinate(CellEngineAtom& CellEngineAtomObjectParam, const UnsignedInt StartXParam, const UnsignedInt StartYParam, const UnsignedInt StartZParam, const UnsignedInt SpaceXParam, const UnsignedInt SpaceYParam, const UnsignedInt SpaceZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam) const
{
    if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
        CellEngineAtomObjectParam.SetAtomPositionsData(CovertToGraphicsCoordinateSelected(StartXParam, SpaceXParam, SizeXParam), CovertToGraphicsCoordinateSelected(StartYParam, SpaceYParam, SizeYParam), CovertToGraphicsCoordinateSelected(StartZParam, SpaceZParam, SizeZParam));
    else
    if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
        CellEngineAtomObjectParam.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceXParam), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceYParam), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceZParam));
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::SetParticleParametersToDraw(CellEngineAtom& TempAtomObject, Particle& ParticleObject)
{
    try
    {
        TempAtomObject.EntityId = ParticleObject.EntityId;
        auto ParticleKindObject = ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId);
        TempAtomObject.AtomColor = ParticleKindObject.AtomColor;
        ParticleObject.ParticleColor = ParticleKindObject.ParticleColor;
        ParticleObject.RandomParticleKindColor = (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true ? CellEngineColorsObject.GetDNAorRNAColor(TempAtomObject.EntityId, ParticleObject.ChainId) : ParticleKindObject.RandomParticleColor);

        if (CellEngineConfigDataObject.DNAPaired == true)
            if (CellEngineUseful::IsDNA(TempAtomObject.EntityId) == true && ((CellEngineConfigDataObject.GenomeReadFromFile == false && ParticleObject.GenomeIndex == 0) || (CellEngineConfigDataObject.GenomeReadFromFile == true && (ParticleObject.Prev != nullptr && ParticleObject.Next != nullptr && ParticleObject.PairedNucleotidePtr == nullptr))))
                ParticleObject.RandomParticleKindColor = ParticleObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Purple));

        if (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true && ((CellEngineConfigDataObject.GenomeReadFromFile == false && ParticleObject.GenomeIndex == 0) || (CellEngineConfigDataObject.GenomeReadFromFile == true && (ParticleObject.Prev == nullptr || ParticleObject.Next == nullptr))))
            ParticleObject.RandomParticleKindColor = ParticleObject.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::Blue));

        if (CellEngineUseful::IsDNAorRNA(TempAtomObject.EntityId) == true && (CellEngineConfigDataObject.GenomeReadFromFile == true))
        {
            ParticleObject.GenomeIndexPrev = (ParticleObject.Prev != nullptr ? ParticleObject.Prev->GenomeIndex : 0);
            ParticleObject.GenomeIndexNext = (ParticleObject.Next != nullptr ? ParticleObject.Next->GenomeIndex : 0);
        }
        ParticleObject.Nucleotide = ((CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.ChainId) : '0');
    }
    CATCH("setting particle parameters to draw")
};

uint32_t AtomTotalIndex = 0;
constexpr UnsignedInt NumberOfSectors = 40;
constexpr UnsignedInt MaxNumberOfSectors = 16;

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GenerateVoxelsForGPUParallel(const UnsignedInt MainPosX, const UnsignedInt MainPosY, const UnsignedInt MainPosZ, const SimulationSpaceVoxel SimulationSpaceVoxelObject, const SimulationSpaceVoxel LastSimulationSpaceVoxel, UnsignedInt& AtomsCounter, const CellEngineAtom& TempAtomObject, const Particle& ParticleObject)
{
    try
    {
        //cout << "K = " << MainPosX << " " << MainPosY << " " << MainPosZ << endl;

        if (LastSimulationSpaceVoxel != SimulationSpaceVoxelObject)
        {
            GPUParticle GPUParticleObject;

            GPUParticleObject.AtomOffset = AtomOffsetInSectors[MainPosX][MainPosY][MainPosZ];
            GPUParticleObject.AtomCount = AtomsCounter + 1;

            GPUParticlesInSectors[MainPosX][MainPosY][MainPosZ].emplace_back(std::move(GPUParticleObject));

            AtomOffsetInSectors[MainPosX][MainPosY][MainPosZ] += (AtomsCounter + 1);

            AtomsCounter = 0;
        }
        else
        {
            if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
            {
                GPUAtomLocal GPUAtomLocalObject;
                GPUAtomLocalObject.ParticleSectorXIndex = 0;
                GPUAtomLocalObject.ParticleSectorYIndex = 0;
                GPUAtomLocalObject.ParticleSectorZIndex = 0;
                GPUAtomLocalObject.Index = ParticleObject.Index;
                GPUAtomLocalObject.AtomOffset = AtomTotalIndex;
                GPUAtomsLocalInSectors[MainPosX][MainPosY][MainPosZ].emplace_back(std::move(GPUAtomLocalObject));
            }

            GPUAtom GPUAtomObject;

            GPUAtomObject.X = TempAtomObject.X;
            GPUAtomObject.Y = TempAtomObject.Y;
            GPUAtomObject.Z = TempAtomObject.Z;

            const auto ParticleColor = CellEngineConfigDataObject.RenderTheWholePickedParticleInOnePickingColorForVoxelSpace == false ? CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(TempAtomObject, ParticleObject, ChosenAtomObjectIndex == AtomTotalIndex)) : CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(TempAtomObject, ParticleObject, ChosenParticleObject.Index == ParticleObject.Index));

            GPUAtomObject.ColorR = ParticleColor.X();
            GPUAtomObject.ColorG = ParticleColor.Y();
            GPUAtomObject.ColorB = ParticleColor.Z();

            GPUAtomsInSectors[MainPosX][MainPosY][MainPosZ].emplace_back(std::move(GPUAtomObject));

            AtomsCounter++;
            AtomTotalIndex++;
        }
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GenerateVoxelsForGPU(const SimulationSpaceVoxel SimulationSpaceVoxelObject, const SimulationSpaceVoxel LastSimulationSpaceVoxel, UnsignedInt& AtomsCounter, const CellEngineAtom& TempAtomObject, const Particle& ParticleObject)
{
    try
    {
        if (LastSimulationSpaceVoxel != SimulationSpaceVoxelObject)
        {
            GPUParticle GPUParticleObject;

            GPUParticleObject.AtomOffset = AtomOffsetTotal;
            GPUParticleObject.AtomCount = AtomsCounter;

            GPUParticles[ParticlesOffsetTotal] = std::move(GPUParticleObject);

            ParticlesOffsetTotal++;

            AtomOffsetTotal += AtomsCounter;

            AtomsCounter = 0;
        }
        else
        {
            if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
            {
                GPUAtomLocal GPUAtomLocalObject;
                GPUAtomLocalObject.ParticleSectorXIndex = 0;
                GPUAtomLocalObject.ParticleSectorYIndex = 0;
                GPUAtomLocalObject.ParticleSectorZIndex = 0;
                GPUAtomLocalObject.Index = ParticleObject.Index;
                GPUAtomLocalObject.AtomOffset = AtomTotalIndex;
                GPUAtomsLocal[AtomLocalOffsetTotal] = std::move(GPUAtomLocalObject);
                AtomLocalOffsetTotal++;
            }

            GPUAtom GPUAtomObject;

            GPUAtomObject.X = TempAtomObject.X;
            GPUAtomObject.Y = TempAtomObject.Y;
            GPUAtomObject.Z = TempAtomObject.Z;

            const auto ParticleColor = CellEngineConfigDataObject.RenderTheWholePickedParticleInOnePickingColorForVoxelSpace == false ? CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(TempAtomObject, ParticleObject, ChosenAtomObjectIndex == AtomTotalIndex)) : CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(TempAtomObject, ParticleObject, ChosenParticleObject.Index == ParticleObject.Index));

            GPUAtomObject.ColorR = ParticleColor.X();
            GPUAtomObject.ColorG = ParticleColor.Y();
            GPUAtomObject.ColorB = ParticleColor.Z();

            GPUAtoms[AtomTotalIndex] = std::move(GPUAtomObject);

            AtomsCounter++;
            AtomTotalIndex++;
        }
    }
    CATCH("");
}

typedef SimulationSpaceVoxel (*PointerToSpace_2048_2048_2048)[2048][2048][2048];

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderSelectedSpace(const vmath::mat4& ViewMatrix, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, CellEngineAtom& TempAtomObject)
{
    try
    {
        const auto SpacePointer = static_cast<PointerToSpace_2048_2048_2048>(CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SpacePointer);

        UnsignedInt AtomsCounter = 0;

        for (UnsignedInt PosX = XStartParam; PosX < XStartParam + XSizeParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam; PosY < YStartParam + YSizeParam; PosY += YStepParam)
            {
                //SimulationSpaceVoxel LastSimulationSpaceVoxel = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSpaceVoxelForOuterClass(PosX, PosY, ZStartParam);
                SimulationSpaceVoxel LastSimulationSpaceVoxel = (*SpacePointer)[PosX][PosY][ZStartParam];

                for (UnsignedInt PosZ = ZStartParam; PosZ < ZStartParam + ZSizeParam; PosZ += ZStepParam)
                    if (PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                    {
                        //const SimulationSpaceVoxel SimulationSpaceVoxelObject = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSpaceVoxelForOuterClass(PosX, PosY, PosZ);
                        const SimulationSpaceVoxel SimulationSpaceVoxelObject = (*SpacePointer)[PosX][PosY][PosZ];

                        if (SimulationSpaceVoxelObject != CellEngineParticlesVoxelsOperations::GetZeroSimulationSpaceVoxel())
                        {
                            if (auto FoundParticleIter = CellEngineDataFileObjectPointer->GetParticleIteratorFromIndex(SimulationSpaceVoxelObject); FoundParticleIter != CellEngineDataFileObjectPointer->GetParticleEnd())
                            {
                                Particle& ParticleObject = FoundParticleIter->second;

                                //if (RenderObject(TempAtomObject, ParticleObject, ViewMatrix, true, false, true, false) == true)
                                {
                                    ConvertAtomPosToGraphicCoordinate(TempAtomObject, XStartParam, YStartParam, ZStartParam, PosX, PosY, PosZ, XSizeParam, YSizeParam, ZSizeParam);

                                    if (DrawEmptyVoxels == false || SimulationSpaceVoxelObject != 0)
                                        SetParticleParametersToDraw(TempAtomObject, ParticleObject);

                                    if (DrawEmptyVoxels == true || (DrawEmptyVoxels == false && SimulationSpaceVoxelObject != 0 && ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true))
                                    {
                                        if ((SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull && CellEngineConfigDataObject.ViewPositionZ <= CellEngineConfigDataObject.Distance + 700) || SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
                                            GenerateVoxelsForGPU(SimulationSpaceVoxelObject, LastSimulationSpaceVoxel, AtomsCounter, TempAtomObject, ParticleObject);
                                        else
                                        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull && CellEngineConfigDataObject.ViewPositionZ > CellEngineConfigDataObject.Distance + 700)
                                            GenerateVoxelsForGPUParallel(XStartParam >> 6, YStartParam >> 6, ZStartParam >> 6, SimulationSpaceVoxelObject, LastSimulationSpaceVoxel, AtomsCounter, TempAtomObject, ParticleObject);
                                            //GenerateVoxelsForGPUParallel(XStartParam / SelectionStepX, YStartParam / SelectionStepY, ZStartParam / SelectionStepZ, SimulationSpaceVoxelObject, LastSimulationSpaceVoxel, AtomsCounter, TempAtomObject, ParticleObject);

                                        LastSimulationSpaceVoxel = SimulationSpaceVoxelObject;
                                    }
                                }
                            }
                            else
                            if (DrawError ==  true)
                                LoggersManagerObject.LogError(STREAM("Try to draw the particle from not existing index = " << SimulationSpaceVoxelObject));
                        }
                        else
                        if (DrawEmptyVoxels == true)
                        {
                            Particle ParticleObject;

                            TempAtomObject.AtomColor = ParticleObject.ParticleColor = ParticleObject.UniqueParticleColor = ParticleObject.RandomParticleKindColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(vmath::FromVec4ToVec3(sb7::color::DeepSkyBlue));

                            ConvertAtomPosToGraphicCoordinate(TempAtomObject, XStartParam, YStartParam, ZStartParam, PosX, PosY, PosZ, XSizeParam, YSizeParam, ZSizeParam);

                            GenerateVoxelsForGPU(SimulationSpaceVoxelObject, LastSimulationSpaceVoxel, AtomsCounter, TempAtomObject, ParticleObject);

                            LastSimulationSpaceVoxel = SimulationSpaceVoxelObject;
                        }
                    }
            }
    }
    CATCH("rendering selected voxel simulation space")
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderSpace(const vmath::mat4& ViewMatrix)
{
    //SpacePointer = static_cast<Space_1024_1024_1024>(CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->SpacePointer);

    try
    {
        AtomTotalIndex = 0;

        CellEngineAtom TempAtomObject;

        //CellEngineConfigDataObject.LoadOfAtomsStep > 10 ? CellEngineConfigDataObject.LoadOfAtomsStep = 4 : 1;

        UnsignedInt NumberOfRenderedSelectedSpaces = 0;

        cout << "P = " << SelectionStartXPos << " " << SelectionStartYPos << " " << SelectionStartZPos << " " << SelectionStepX << " " << SelectionStepY << " " << SelectionStepZ << " " << SelectionSizeX << " " << SelectionSizeY << " " << SelectionSizeZ << endl;

        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull && CellEngineConfigDataObject.ViewPositionZ <= CellEngineConfigDataObject.Distance + 700)
        {
            const auto start_time111 = chrono::high_resolution_clock::now();

            for (UnsignedInt PosX = SelectionStartXPos; PosX < SelectionSizeX; PosX += SelectionStepX)
                for (UnsignedInt PosY = SelectionStartYPos; PosY < SelectionSizeY; PosY += SelectionStepY)
                    for (UnsignedInt PosZ = SelectionStartZPos; PosZ < SelectionSizeZ; PosZ += SelectionStepZ)
                    {
                        TempAtomObject.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosX), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosY), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosZ));

                        if (RenderObject(TempAtomObject, Particle(), ViewMatrix, true, false, true, false) == true)
                        {
                            RenderSelectedSpace(ViewMatrix, PosX, PosY, PosZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, TempAtomObject);
                            NumberOfRenderedSelectedSpaces++;
                        }
                    }

            const auto stop_time111 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);
        }
        else
        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull && CellEngineConfigDataObject.ViewPositionZ > CellEngineConfigDataObject.Distance + 700)
        {
            for (UnsignedInt PosX = 0; PosX < NumberOfSectors; PosX++)
                for (UnsignedInt PosY = 0; PosY < NumberOfSectors; PosY++)
                    for (UnsignedInt PosZ = 0; PosZ < NumberOfSectors; PosZ++)
                    {
                        GPUParticlesInSectors[PosX][PosY][PosZ].clear();
                        GPUAtomsInSectors[PosX][PosY][PosZ].clear();
                        GPUAtomsLocalInSectors[PosX][PosY][PosZ].clear();
                        AtomOffsetInSectors[PosX][PosY][PosZ] = 0;
                        TempAtomObjectInSectors[PosX][PosY][PosZ] = {};
                    }

            const auto start_time111 = chrono::high_resolution_clock::now();

            #pragma omp parallel for collapse(3) num_threads(256) default(none) shared(CellEngineConfigDataObject, ViewMatrix, TempAtomObject, NumberOfRenderedSelectedSpaces, AtomOffsetInSectors, GPUAtomsInSectors, GPUParticlesInSectors, GPUAtomsLocalInSectors, TempAtomObjectInSectors) schedule(dynamic)
            for (UnsignedInt PosX = 0; PosX < MaxNumberOfSectors; PosX++)
                for (UnsignedInt PosY = 0; PosY < MaxNumberOfSectors; PosY++)
                    for (UnsignedInt PosZ = 0; PosZ < MaxNumberOfSectors; PosZ++)
                    {
                        TempAtomObjectInSectors[PosX][PosY][PosZ].SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosX * SelectionStepX), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosY * SelectionStepY), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosZ * SelectionStepZ));

                        if (RenderObject(TempAtomObjectInSectors[PosX][PosY][PosZ], Particle(), ViewMatrix, true, false, true, false) == true)
                        {
                            //RenderSelectedSpace(PosX * SelectionStepX, PosY * SelectionStepY, PosZ * SelectionStepZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, TempAtomObject);
                            RenderSelectedSpace(ViewMatrix, PosX * SelectionStepX, PosY * SelectionStepY, PosZ * SelectionStepZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, TempAtomObjectInSectors[PosX][PosY][PosZ]);
                            NumberOfRenderedSelectedSpaces++;
                        }
                    }

            const auto stop_time111 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);

            LoggersManagerObject.LogOnlyToConsole(STREAM("END OF PARALLEL -> ParticlesOffsetTotal = " << ParticlesOffsetTotal << " AtomsOffsetTotal = " << AtomOffsetTotal));

            const auto start_time114 = chrono::high_resolution_clock::now();

            for (UnsignedInt PosX = 0; PosX < MaxNumberOfSectors; PosX++)
                for (UnsignedInt PosY = 0; PosY < MaxNumberOfSectors; PosY++)
                    for (UnsignedInt PosZ = 0; PosZ < MaxNumberOfSectors; PosZ++)
                    {
                        for (auto& GPUParticleInSectors : GPUParticlesInSectors[PosX][PosY][PosZ])
                        {
                            GPUParticleInSectors.AtomOffset += AtomOffsetTotal;
                            GPUParticles[ParticlesOffsetTotal++] = std::move(GPUParticleInSectors);
                        }

                        memcpy(&GPUAtoms[AtomOffsetTotal], GPUAtomsInSectors[PosX][PosY][PosZ].data(), GPUAtomsInSectors[PosX][PosY][PosZ].size() * sizeof(GPUAtom));

                        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                            memcpy(&GPUAtomsLocal[AtomLocalOffsetTotal], GPUAtomsLocalInSectors[PosX][PosY][PosZ].data(), GPUAtomsLocalInSectors[PosX][PosY][PosZ].size() * sizeof(GPUAtomLocal));

                        AtomOffsetTotal += GPUAtomsInSectors[PosX][PosY][PosZ].size();
                        AtomLocalOffsetTotal += GPUAtomsLocalInSectors[PosX][PosY][PosZ].size();
                    }

            const auto stop_time114 = chrono::high_resolution_clock::now();

            ExecutionDurationTimeForCopyingParticlesToGraphicMemory3 += chrono::duration(stop_time114 - start_time114);
        }
        else
        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
        {
            UnsignedInt SubStartPos = 0;
            if (CellEngineConfigDataObject.SelectedSpaceStartParametersDrawTypesObject == CellEngineConfigData::SelectedSpaceStartParametersDrawTypes::DrawFromCenter)
                SubStartPos = SelectionSizeX / 2;

            RenderSelectedSpace(ViewMatrix, SelectionStartXPos - SubStartPos, SelectionStartYPos - SubStartPos, SelectionStartZPos - SubStartPos, SelectionStepX, SelectionStepY, SelectionStepZ, SelectionSizeX, SelectionSizeY, SelectionSizeY, TempAtomObject);
        }

        LoggersManagerObject.Log(STREAM("NumberOfRenderedSelectedSpaces = " << NumberOfRenderedSelectedSpaces));
    }
    CATCH("rendering voxel simulation space");
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::DrawChosenAtomUsingStencilBuffer1(const GLuint ChosenParticleCenterIndex)
{
    try
    {
        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
        {
            if (ChosenParticleCenterIndex > 0)
            {
                if (ChosenParticleCenterIndex < GPUAtomsLocal.size())
                {
                    if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticleIteratorFromIndex(GPUAtomsLocal[ChosenParticleCenterIndex].Index); ParticleIter != CellEngineDataFileObjectPointer->GetParticleEnd())
                    {
                        cout << "PICKING OBJECT = " << std::to_string(ChosenParticleCenterIndex) + " " + to_string(GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset) + " " + to_string(ParticleIter->second.ListOfAtoms.size()) << endl;

                        ChosenParticleObject = ParticleIter->second;
                        ChosenAtomObjectIndex = GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset;
                    }
                    else
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + std::to_string(ChosenParticleCenterIndex));
                }

                PrintAtomDescriptionOnScreen(ChosenAtomObject, ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using buffer")
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetStartCenterPoint()
{
    Center = { 0.0f, 0.0f, 0.0f };
}

#endif