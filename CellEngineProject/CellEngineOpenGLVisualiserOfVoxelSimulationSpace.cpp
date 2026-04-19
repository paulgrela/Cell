
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

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GenerateVoxelsForGPU(const SimulationSpaceVoxel SimulationSpaceVoxelObject, const SimulationSpaceVoxel LastSimulationSpaceVoxel, UnsignedInt& AtomsCounter, UnsignedInt& AtomsLocalCounter, const CellEngineAtom& TempAtomObject, const Particle& ParticleObject, const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XSizeParam, const UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
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

            AtomsLocalCounter = 0;
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

                //GPUAtomLocalObject.AtomOffset = 0;
                GPUAtomLocalObject.AtomOffset = AtomsLocalCounter + 1;
                //GPUAtomLocalObject.AtomOffset = AtomsLocalCounter;
                GPUAtomsLocal[AtomLocalOffsetTotal] = std::move(GPUAtomLocalObject);

                AtomLocalOffsetTotal++;
                AtomsLocalCounter++;
            }

            GPUAtom GPUAtomObject;

            GPUAtomObject.X = TempAtomObject.X;
            GPUAtomObject.Y = TempAtomObject.Y;
            GPUAtomObject.Z = TempAtomObject.Z;

            const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(TempAtomObject, ParticleObject, ChosenParticleObject.Index == ParticleObject.Index && ChosenAtomObjectIndex == AtomsLocalCounter));
            //const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(TempAtomObject, ParticleObject, false));

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

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderSelectedSpace(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, const vmath::mat4& ViewMatrix, CellEngineAtom& TempAtomObject, std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList, UnsignedInt StencilBufferLoopCounter)
{
    try
    {
        UnsignedInt AtomsCounter = 0;
        UnsignedInt AtomsLocalCounter = 0;

        for (UnsignedInt PosX = XStartParam; PosX < XStartParam + XSizeParam; PosX += XStepParam)
            for (UnsignedInt PosY = YStartParam; PosY < YStartParam + YSizeParam; PosY += YStepParam)
            {
                SimulationSpaceVoxel LastSimulationSpaceVoxel = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSpaceVoxelForOuterClass(PosX, PosY, ZStartParam);

                for (UnsignedInt PosZ = ZStartParam; PosZ < ZStartParam + ZSizeParam; PosZ += ZStepParam)
                    if (PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                    {
                        const SimulationSpaceVoxel SimulationSpaceVoxelObject = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSpaceVoxelForOuterClass(PosX, PosY, PosZ);

                        if (SimulationSpaceVoxelObject != CellEngineParticlesVoxelsOperations::GetZeroSimulationSpaceVoxel())
                        {
                            if (auto FoundParticleIter = CellEngineDataFileObjectPointer->GetParticleIteratorFromIndex(SimulationSpaceVoxelObject); FoundParticleIter != CellEngineDataFileObjectPointer->GetParticleEnd())
                            {
                                Particle& ParticleObject = FoundParticleIter->second;

                                if (DrawEmptyVoxels == true || (DrawEmptyVoxels == false && SimulationSpaceVoxelObject != 0 && ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true))
                                {
                                    ConvertAtomPosToGraphicCoordinate(TempAtomObject, XStartParam, YStartParam, ZStartParam, PosX, PosY, PosZ, XSizeParam, YSizeParam, ZSizeParam);

                                    if (DrawEmptyVoxels == false || SimulationSpaceVoxelObject != 0)
                                        SetParticleParametersToDraw(TempAtomObject, ParticleObject);

                                    GenerateVoxelsForGPU(SimulationSpaceVoxelObject, LastSimulationSpaceVoxel, AtomsCounter, AtomsLocalCounter, TempAtomObject, ParticleObject, PosX, PosY, PosZ, XStartParam, YStartParam, ZStartParam, XSizeParam, YSizeParam, ZSizeParam);

                                    LastSimulationSpaceVoxel = SimulationSpaceVoxelObject;

                                    // if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                    // {
                                    //     glStencilFunc(GL_ALWAYS, static_cast<uint8_t>((TemporaryRenderedVoxelsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                    //     TemporaryRenderedVoxelsList.emplace_back(TemporaryRenderedVoxel{ TempAtomObject, FoundParticleIter, PosX, PosY, PosZ });
                                    // }
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

                            GenerateVoxelsForGPU(SimulationSpaceVoxelObject, LastSimulationSpaceVoxel, AtomsCounter, AtomsLocalCounter, TempAtomObject, ParticleObject, PosX, PosY, PosZ, XStartParam, YStartParam, ZStartParam, XSizeParam, YSizeParam, ZSizeParam);

                            LastSimulationSpaceVoxel = SimulationSpaceVoxelObject;
                        }
                    }
            }
    }
    CATCH("rendering selected voxel simulation space")
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderSpace(const vmath::mat4& ViewMatrix)
{
    try
    {
        AtomTotalIndex = 0;

        GLuint PartOfStencilBufferIndex[3];

        CellEngineAtom TempAtomObject;

        std::vector<TemporaryRenderedVoxel> TemporaryRenderedVoxelsList;

        CellEngineConfigDataObject.LoadOfAtomsStep > 10 ? CellEngineConfigDataObject.LoadOfAtomsStep = 4 : 1;

        UnsignedInt NumberOfRenderedSelectedSpaces = 0;

        TemporaryRenderedVoxelsList.clear();

        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
            for (UnsignedInt PosX = SelectionStartXPos; PosX < SelectionSizeX; PosX += SelectionStepX)
                for (UnsignedInt PosY = SelectionStartYPos; PosY < SelectionSizeY; PosY += SelectionStepY)
                    for (UnsignedInt PosZ = SelectionStartZPos; PosZ < SelectionSizeZ; PosZ += SelectionStepZ)
                    {
                        TempAtomObject.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosX), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosY), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosZ));

                        if (RenderObject(TempAtomObject, Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                        {
                            RenderSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, 0);
                            NumberOfRenderedSelectedSpaces++;
                        }
                    }
        else
        {
            UnsignedInt SubStartPos = 0;
            if (CellEngineConfigDataObject.SelectedSpaceStartParametersDrawTypesObject == CellEngineConfigData::SelectedSpaceStartParametersDrawTypes::DrawFromCenter)
                SubStartPos = SelectionSizeX / 2;

            RenderSelectedSpace(SelectionStartXPos - SubStartPos, SelectionStartYPos - SubStartPos, SelectionStartZPos - SubStartPos, SelectionStepX, SelectionStepY, SelectionStepZ, SelectionSizeX, SelectionSizeY, SelectionSizeY, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, 0);
        }

        // if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
        //     glReadPixels(static_cast<GLint>(MousePositionLocal.s.X), static_cast<GLint>(static_cast<float>(Info.WindowHeight) - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);

        LoggersManagerObject.Log(STREAM("NumberOfRenderedSelectedSpaces = " << NumberOfRenderedSelectedSpaces));

        // if (PressedRightMouseButton != 1)
        //     DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedVoxelsList);
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
                    if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticles()[GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorXIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorYIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorZIndex].Particles.find(GPUAtomsLocal[ChosenParticleCenterIndex].Index); ParticleIter != CellEngineDataFileObjectPointer->GetParticles()[GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorXIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorYIndex][GPUAtomsLocal[ChosenParticleCenterIndex].ParticleSectorZIndex].Particles.end())
                    //if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticles()[0][0][0].Particles.find(GPUAtomsLocal[ChosenParticleCenterIndex].Index); ParticleIter != CellEngineDataFileObjectPointer->GetParticles()[0][0][0].Particles.end())
                    //if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticleIteratorFromIndex(GPUAtomsLocal[ChosenParticleCenterIndex].Index); ParticleIter != CellEngineDataFileObjectPointer->GetParticleEnd())
                    {
                        // if (GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset > ParticleIter->second.ListOfVoxels.size())
                        //     throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + std::to_string(ChosenParticleCenterIndex) + " " + to_string(GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset) + " " + to_string(ParticleIter->second.ListOfAtoms.size()));
                        // else
                        {
                            cout << "PICKING OBJECT = " << std::to_string(ChosenParticleCenterIndex) + " " + to_string(GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset) + " " + to_string(ParticleIter->second.ListOfAtoms.size()) << endl;
                            ChosenParticleObject = ParticleIter->second;
                            //cout << "END1" << endl;
                            //ChosenAtomObject = ParticleIter->second.ListOfAtoms[GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset];
                            //ChosenAtomObject = ParticleIter->second.ListOfAtoms[0];
                            //cout << "END2" << endl;
                            //ChosenAtomObjectIndex = GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset;
                            ChosenAtomObjectIndex = 0;
                            //cout << "END" << endl;
                        }
                    }
                    else
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + std::to_string(ChosenParticleCenterIndex));
                }

                CellEngineAtom ChosenAtomObject{};
                PrintAtomDescriptionOnScreen(ChosenAtomObject, ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using buffer")
}

inline void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList)
{
    try
    {
        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
        {
            UnsignedInt ChosenVoxelIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

            if (ChosenVoxelIndex > 0)
            {
                TemporaryRenderedVoxel ChosenObject{};

                if (ChosenVoxelIndex > TemporaryRenderedVoxelsList.size())
                    throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + std::to_string(ChosenVoxelIndex) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(TemporaryRenderedVoxelsList.size()));
                else
                    ChosenObject = TemporaryRenderedVoxelsList[ChosenVoxelIndex];

                SetSaveXYZPositions(ChosenObject.X, ChosenObject.Y, ChosenObject.Z);

                RenderObject(ChosenObject.CellEngineAtomObject, ChosenObject.CellEngineParticleObject->second, ViewMatrix, false, false, false, true, RenderObjectsBool);

                PrintAtomDescriptionOnScreen(ChosenObject.CellEngineAtomObject, ChosenObject.CellEngineParticleObject->second);
            }
        }
    }
    CATCH("choosing atom using stencil buffer")
}

void CellEngineOpenGLVisualiserOfVoxelSimulationSpace::GetStartCenterPoint()
{
    Center = { 0.0f, 0.0f, 0.0f };
}

#endif