#pragma once

#ifndef CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H

#include <string>

#include <sb7.h>
#include <tuple>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineVoxelSimulationSpace.h"

class CellEngineOpenGLVisualiserOfVoxelSimulationSpace : public CellEngineOpenGLVisualiser
{
public:
    enum class VoxelSpaceDrawingTypes : UnsignedInt
    {
        DrawVoxelSpaceFull = 1,
        DrawVoxelSpaceSelected = 2
    };
public:
    VoxelSpaceDrawingTypes SpaceDrawingType = VoxelSpaceDrawingTypes::DrawVoxelSpaceFull;
    bool DrawEmptyVoxels = false;
private:
    UnsignedInt SelectionStartXPos{}, SelectionStartYPos{}, SelectionStartZPos{};
    UnsignedInt SelectionStepX{}, SelectionStepY{}, SelectionStepZ{};
    UnsignedInt SelectionSizeX{}, SelectionSizeY{}, SelectionSizeZ{};
public:
    inline std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetStartPositions()
    {
        return {SelectionStartXPos, SelectionStartYPos, SelectionStartZPos };
    }
    inline std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetSteps()
    {
        return { SelectionStepX, SelectionStepY, SelectionStepZ };
    }
    inline std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetSizes()
    {
        return { SelectionSizeX, SelectionSizeY, SelectionSizeZ };
    }
    inline void SetVoxelSpaceSelection(const UnsignedInt SelectionStartXParam, const UnsignedInt SelectionStartYParam, const UnsignedInt SelectionStartZParam, const UnsignedInt SelectionStepXParam, const UnsignedInt SelectionStepYParam, const UnsignedInt SelectionStepZParam, const UnsignedInt SelectionSizeXParam, const UnsignedInt SelectionSizeYParam, const UnsignedInt SelectionSizeZParam)
    {
        SelectionStartXPos = SelectionStartXParam, SelectionStartYPos = SelectionStartYParam, SelectionStartZPos = SelectionStartZParam;
        SelectionStepX = SelectionStepXParam, SelectionStepY = SelectionStepYParam, SelectionStepZ = SelectionStepZParam;
        SelectionSizeX = SelectionSizeXParam, SelectionSizeY = SelectionSizeYParam, SelectionSizeZ = SelectionSizeZParam;
    }
public:
    CellEngineOpenGLVisualiserOfVoxelSimulationSpace()
    {
        SetVoxelSpaceSelection(CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartXPos, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartZPos, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStepX, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStepY, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStepZ, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionSizeX, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionSizeY, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionSizeZ);
    }
private:
    struct TemporaryRenderedVoxel
    {
        CellEngineAtom CellEngineAtomObject;
        UnsignedInt X, Y, Z;
    };
    UnsignedInt SaveXMousePosition{}, SaveYMousePosition{}, SaveZMousePosition{};
public:
    inline void SaveVoxelPositionChosenByMouse()
    {
        SelectionStartXPos = SaveXMousePosition;
        SelectionStartYPos = SaveYMousePosition;
        SelectionStartZPos = SaveZMousePosition;
    }
private:
    inline void SetSaveXYZPositions(const UnsignedInt SaveXParam, const UnsignedInt SaveYParam, const UnsignedInt SaveZParam)
    {
        SaveXMousePosition = SaveXParam;
        SaveYMousePosition = SaveYParam;
        SaveZMousePosition = SaveZParam;
    }
private:
    static inline float CovertToGraphicsCoordinateSelected(const UnsignedInt StartParam, const UnsignedInt SpaceP, const UnsignedInt SizeParam)
    {
        return ((static_cast<float>((StartParam + SizeParam) - SpaceP) - static_cast<float>(SizeParam) / 2) * 4);
    }
private:
    inline void ConvertAtomPosToGraphicCoordinate(CellEngineAtom& CellEngineAtomObjectParam, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt SpaceXP, const UnsignedInt SpaceYP, const UnsignedInt SpaceZP, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam) const
    {
        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceSelected)
            CellEngineAtomObjectParam.SetAtomPositionsData(CovertToGraphicsCoordinateSelected(XStartParam, SpaceXP, XSizeParam), CovertToGraphicsCoordinateSelected(YStartParam, SpaceYP, YSizeParam), CovertToGraphicsCoordinateSelected(ZStartParam, SpaceZP, ZSizeParam));
        else
        if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
            CellEngineAtomObjectParam.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceXP), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceYP), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceZP));
    }
private:
    void RenderSelectedSpace(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam, UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, CellEngineAtom& TempAtomObject, std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList, UnsignedInt StencilBufferLoopCounter)
    {
        try
        {
            for (UnsignedInt PosX = XStartParam; PosX < XStartParam + XSizeParam; PosX += XStepParam)
                for (UnsignedInt PosY = YStartParam; PosY < YStartParam + YSizeParam; PosY += YStepParam)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZStartParam + ZSizeParam; PosZ += ZStepParam)
                        if (PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                        {
                            SimulationSpaceVoxel SimulationSpaceVoxelObject = CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GetSimulationSpaceVoxel(PosX, PosY, PosZ);

                            if (DrawEmptyVoxels == true || (DrawEmptyVoxels == false && SimulationSpaceVoxelObject.EntityId != 0 && CheckVisibilityOfParticles(SimulationSpaceVoxelObject.EntityId) == true))
                            {
                                ConvertAtomPosToGraphicCoordinate(TempAtomObject, XStartParam, YStartParam, ZStartParam, PosX, PosY, PosZ, XSizeParam, YSizeParam, ZSizeParam);

                                if (DrawEmptyVoxels == false || (DrawEmptyVoxels == true && SimulationSpaceVoxelObject.EntityId != 0))
                                {
                                    TempAtomObject.EntityId = SimulationSpaceVoxelObject.EntityId;
                                    auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds[CellEngineConfigDataObject.ParticlesKindsPos.find(SimulationSpaceVoxelObject.EntityId)->second];
                                    TempAtomObject.AtomColor = ParticleKindObject.AtomColor;
                                    TempAtomObject.ParticleColor = ParticleKindObject.ParticleColor;
                                    TempAtomObject.RandomParticleColor = (CellEngineConfigDataObject.IsDNAorRNA(TempAtomObject.EntityId) ? CellEngineConfigDataObject.GetDNAorRNAColor(TempAtomObject.EntityId, SimulationSpaceVoxelObject.ChainId) : ParticleKindObject.RandomParticleColor);
                                }
                                else
                                if (DrawEmptyVoxels == true)
                                    TempAtomObject.AtomColor = TempAtomObject.ParticleColor = TempAtomObject.RandomParticleColor = GetVector3FormVMathVec3(sb7::FromVec4ToVec3(sb7::color::DeepSkyBlue));

                                if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                {
                                    glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedVoxelsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                    TemporaryRenderedVoxelsList.emplace_back(TemporaryRenderedVoxel{TempAtomObject, PosX, PosY, PosZ });
                                }

                                RenderObject(TempAtomObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                            }
                        }
        }
        CATCH("rendering selected voxel simulation space")
    }
private:
    void RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix) override
    {
        try
        {
            std::lock_guard<std::mutex> LockGuardObject{RenderMenuAndVoxelSimulationSpaceMutexObject};

            GLuint PartOfStencilBufferIndex[3];

            CellEngineAtom TempAtomObject;
            TempAtomObject.Visible = true;

            NumberOfAllRenderedAtoms = 0;

            std::vector<TemporaryRenderedVoxel> TemporaryRenderedVoxelsList;

            CellEngineConfigDataObject.LoadOfAtomsStep > 10 ? CellEngineConfigDataObject.LoadOfAtomsStep = 4 : 1;

            for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
            {
                NumberOfAllRenderedAtoms = 0;

                TemporaryRenderedVoxelsList.clear();

                if (SpaceDrawingType == VoxelSpaceDrawingTypes::DrawVoxelSpaceFull)
                    for (UnsignedInt PosX = SelectionStartXPos; PosX < SelectionSizeX; PosX += SelectionStepX)
                        for (UnsignedInt PosY = SelectionStartYPos; PosY < SelectionSizeY; PosY += SelectionStepY)
                            for (UnsignedInt PosZ = SelectionStartZPos; PosZ < SelectionSizeZ; PosZ += SelectionStepZ)
                            {
                                TempAtomObject.SetAtomPositionsData(CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosX), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosY), CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(PosZ));

                                if (RenderObject(TempAtomObject, ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                                {
                                    NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;
                                    RenderSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);
                                }
                            }
                else
                    RenderSelectedSpace(SelectionStartXPos, SelectionStartYPos, SelectionStartZPos, SelectionStepX, SelectionStepY, SelectionStepZ, SelectionSizeX, SelectionSizeY, SelectionSizeY, NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);

                if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                    glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
            }

            if (PressedRightMouseButton != 1)
                DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedVoxelsList);
        }
        CATCH("rendering voxel simulation space");
    }
private:
    inline void DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList)
    {
        try
        {
            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
            {
                UnsignedInt ChosenVoxelIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

                if (ChosenVoxelIndex > 0)
                {
                    TemporaryRenderedVoxel ChosenParticleObject{};

                    if (ChosenVoxelIndex > TemporaryRenderedVoxelsList.size())
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + std::to_string(ChosenVoxelIndex) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(TemporaryRenderedVoxelsList.size()));
                    else
                        ChosenParticleObject = TemporaryRenderedVoxelsList[ChosenVoxelIndex];

                    SetSaveXYZPositions(ChosenParticleObject.X, ChosenParticleObject.Y, ChosenParticleObject.Z);

                    RenderObject(ChosenParticleObject.CellEngineAtomObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool);

                    PrintAtomDescriptionOnScreen(ChosenParticleObject.CellEngineAtomObject);
                }
            }
        }
        CATCH("choosing atom using stencil buffer")
    }
public:
    void GetStartCenterPoint() override
    {
        Center = { 0.0f, 0.0f, 0.0f };
    }
public:
    void GetMemoryForBondsBetweenAtomsToDraw() override
    {
    }
public:
    void DrawBondsForParticlesCenters(std::vector<std::pair<UnsignedInt, UnsignedInt>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix) override
    {
    }
public:
    static inline std::mutex RenderMenuAndVoxelSimulationSpaceMutexObject;
};

#endif
