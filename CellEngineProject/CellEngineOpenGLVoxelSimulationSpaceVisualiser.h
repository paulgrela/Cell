#pragma once

#ifndef CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H

#include <string>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"

class CellEngineOpenGLVoxelSimulationSpaceVisualiser : public CellEngineOpenGLVisualiser
{
public:
    UnsignedIntType XStart = 0, YStart = 0, ZStart = 0;
    UnsignedIntType XStep = 64, YStep = 64, ZStep = 64;
    UnsignedIntType XEnd = 1024, YEnd = 1024, ZEnd = 1024;
private:
    void SetVoxelSpaceSelection(UnsignedIntType XStartParam, UnsignedIntType YStartParam, UnsignedIntType ZStartParam, UnsignedIntType XStepParam, UnsignedIntType YStepParam, UnsignedIntType ZStepParam, UnsignedIntType XEndParam, UnsignedIntType YEndParam, UnsignedIntType ZEndParam)
    {
        XStart = XStartParam, YStart = YStartParam, ZStart = ZStartParam;
        XStep = XStepParam, YStep = YStepParam, ZStep = ZStepParam;
        XEnd = XEndParam, YEnd = YEndParam, ZEnd = ZEndParam;
    }
private:
    void RenderSelectedSpace(UnsignedIntType& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, CellEngineAtom& TempAtomObject, std::vector<CellEngineAtom>& TemporaryRenderedVoxelsList, UnsignedIntType StencilBufferLoopCounter)
    {
        try
        {
            for (UnsignedIntType SpaceXP = XStart; SpaceXP < XEnd; SpaceXP += XStep)
                for (UnsignedIntType SpaceYP = YStart; SpaceYP < YEnd; SpaceYP += YStep)
                    for (UnsignedIntType SpaceZP = ZStart; SpaceZP < ZEnd; SpaceZP += ZStep)
                        if (CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP] != 0)
                        {
                            TempAtomObject.X = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceXP));
                            TempAtomObject.Y = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceYP));
                            TempAtomObject.Z = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceZP));
                            TempAtomObject.EntityId = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP];
                            auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds[CellEngineConfigDataObject.ParticlesKindsPos.find(CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP])->second];
                            TempAtomObject.AtomColor = GetColor(ParticleKindObject, false);
                            TempAtomObject.ParticleColor = GetColor(ParticleKindObject, false);
                            TempAtomObject.RandomParticleColor = GetColor(ParticleKindObject, false);

                            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                            {
                                glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedVoxelsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                TemporaryRenderedVoxelsList.emplace_back(TempAtomObject);
                            }

                            RenderObject(TempAtomObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                        }
        }
        CATCH("rendering selected voxel simulation space")
    }
private:
    void RenderSpace(UnsignedIntType& NumberOfAllRenderedAtoms, UnsignedIntType& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix, const Point2fT& MousePositionLocal) override
    {
        try
        {
            GLuint PartOfStencilBufferIndex[3];

            CellEngineAtom TempAtomObject;
            TempAtomObject.Visible = true;

            NumberOfAllRenderedAtoms = 0;

            std::vector<CellEngineAtom> TemporaryRenderedVoxelsList;

            CellEngineConfigDataObject.LoadOfAtomsStep > 10 ? CellEngineConfigDataObject.LoadOfAtomsStep = 4 : 1;

            for (UnsignedIntType StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
            {
                NumberOfAllRenderedAtoms = 0;

                TemporaryRenderedVoxelsList.clear();

                for (IntType SpaceX = 0; SpaceX < 1024; SpaceX += 64)
                    for (IntType SpaceY = 0; SpaceY < 1024; SpaceY += 64)
                        for (IntType SpaceZ = 0; SpaceZ < 1024; SpaceZ += 64)
                        {
                            TempAtomObject.X = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceX);
                            TempAtomObject.Y = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceY);
                            TempAtomObject.Z = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceZ);
                            if (RenderObject(TempAtomObject, ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                            {
                                NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                                SetVoxelSpaceSelection(SpaceX, SpaceY, SpaceZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, SpaceX + 64, SpaceY + 64, SpaceZ + 64);

                                RenderSelectedSpace(NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);

                            }
                        }

                if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                    glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
            }

            DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedVoxelsList);
        }
        CATCH("rendering voxel simulation space");
    }
public:
    inline void DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedIntType& NumberOfAllRenderedAtoms, const std::vector<CellEngineAtom>& TemporaryRenderedVoxelsList)
    {
        try
        {
            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
            {
                UnsignedIntType ChosenVoxelIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

                if (ChosenVoxelIndex > 0)
                {
                    CellEngineAtom ChosenParticleObject{};

                    if (ChosenVoxelIndex > TemporaryRenderedVoxelsList.size())
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + std::to_string(ChosenVoxelIndex) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(TemporaryRenderedVoxelsList.size()));
                    else
                        ChosenParticleObject = TemporaryRenderedVoxelsList[ChosenVoxelIndex];

                    RenderObject(ChosenParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool);

                    PrintAtomDescriptionOnScreen(ChosenParticleObject);
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
    void DrawBondsForParticlesCenters(std::vector<std::pair<UnsignedIntType, UnsignedIntType>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix) override
    {
    }
};

#endif
