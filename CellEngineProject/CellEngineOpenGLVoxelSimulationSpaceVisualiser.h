#pragma once

#ifndef CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H

#include <string>

#include <sb7.h>
#include <tuple>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"

class CellEngineOpenGLVoxelSimulationSpaceVisualiser : public CellEngineOpenGLVisualiser
{
private:
    UnsignedIntType XStart = 512, YStart = 512, ZStart = 512;
    UnsignedIntType XStep = 1, YStep = 1, ZStep = 1;
    UnsignedIntType XSize = 16, YSize = 16, ZSize = 16;
public:
    bool DrawEmptyVoxels = false;
    bool OnlySelectedSpace = false;
public:
    std::tuple<UnsignedIntType, UnsignedIntType, UnsignedIntType> GetStartPositions()
    {
        return { XStart, YStart, ZStart };
    }
    std::tuple<UnsignedIntType, UnsignedIntType, UnsignedIntType> GetSteps()
    {
        return { XStep, YStep, ZStep };
    }
    std::tuple<UnsignedIntType, UnsignedIntType, UnsignedIntType> GetSizes()
    {
        return { XSize, YSize, ZSize };
    }
    void SetVoxelSpaceSelection(UnsignedIntType XStartParam, UnsignedIntType YStartParam, UnsignedIntType ZStartParam, UnsignedIntType XStepParam, UnsignedIntType YStepParam, UnsignedIntType ZStepParam, UnsignedIntType XSizeParam, UnsignedIntType YSizeParam, UnsignedIntType ZSizeParam)
    {
        XStart = XStartParam, YStart = YStartParam, ZStart = ZStartParam;
        XStep = XStepParam, YStep = YStepParam, ZStep = ZStepParam;
        XSize = XSizeParam, YSize = YSizeParam, ZSize = ZSizeParam;
    }
private:
    void RenderSelectedSpace(UnsignedIntType XStartParam, UnsignedIntType YStartParam, UnsignedIntType ZStartParam, UnsignedIntType XStepParam, UnsignedIntType YStepParam, UnsignedIntType ZStepParam, UnsignedIntType XSizeParam, UnsignedIntType YSizeParam, UnsignedIntType ZSizeParam, UnsignedIntType& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, CellEngineAtom& TempAtomObject, std::vector<CellEngineAtom>& TemporaryRenderedVoxelsList, UnsignedIntType StencilBufferLoopCounter)
    {
        try
        {
            for (UnsignedIntType SpaceXP = XStartParam; SpaceXP < XStartParam + XSizeParam; SpaceXP += XStepParam)
                for (UnsignedIntType SpaceYP = YStartParam; SpaceYP < YStartParam + YSizeParam; SpaceYP += YStepParam)
                    for (UnsignedIntType SpaceZP = ZStartParam; SpaceZP < ZStartParam + ZSizeParam; SpaceZP += ZStepParam)
                        if (DrawEmptyVoxels == true || (DrawEmptyVoxels == false && CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP] != 0))
                        {
                            TempAtomObject.X = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceXP));
                            TempAtomObject.Y = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceYP));
                            TempAtomObject.Z = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceZP));

                            if (DrawEmptyVoxels == false || (DrawEmptyVoxels == true && CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP] != 0))
                            {
                                TempAtomObject.EntityId = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP];
                                auto ParticleKindObject = CellEngineConfigDataObject.ParticlesKinds[CellEngineConfigDataObject.ParticlesKindsPos.find(CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP])->second];
                                TempAtomObject.AtomColor = GetColor(ParticleKindObject, false);
                                TempAtomObject.ParticleColor = GetColor(ParticleKindObject, false);
                                TempAtomObject.RandomParticleColor = GetColor(ParticleKindObject, false);
                            }
                            else
                            if (DrawEmptyVoxels == true)
                            {
                                TempAtomObject.AtomColor = GetVector3FormVMathVec3(sb7::FromVec4ToVec3(sb7::color::DeepSkyBlue));
                                TempAtomObject.ParticleColor = GetVector3FormVMathVec3(sb7::FromVec4ToVec3(sb7::color::DeepSkyBlue));
                                TempAtomObject.RandomParticleColor = GetVector3FormVMathVec3(sb7::FromVec4ToVec3(sb7::color::DeepSkyBlue));
                            }

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

                if (OnlySelectedSpace == false)
                    for (IntType SpaceX = XStart; SpaceX < XSize; SpaceX += XStep)
                        for (IntType SpaceY = YStart; SpaceY < YSize; SpaceY += YStep)
                            for (IntType SpaceZ = ZStart; SpaceZ < ZSize; SpaceZ += ZStep)
                            {
                                TempAtomObject.X = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceX);
                                TempAtomObject.Y = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceY);
                                TempAtomObject.Z = CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(SpaceZ);
                                if (RenderObject(TempAtomObject, ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                                {
                                    NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                                    RenderSelectedSpace(SpaceX, SpaceY, SpaceZ, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, CellEngineConfigDataObject.LoadOfAtomsStep, 64, 64, 64, NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);
                                }
                            }
                else
                    RenderSelectedSpace(XStart, YStart, ZStart, XStep, YStep, ZStep, XSize, YSize, YSize, NumberOfAllRenderedAtoms, ViewMatrix, TempAtomObject, TemporaryRenderedVoxelsList, StencilBufferLoopCounter);

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
