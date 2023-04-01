#pragma once

#ifndef CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H

#include <string>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"

class CellEngineOpenGLVoxelSimulationSpaceVisualiser : public CellEngineOpenGLVisualiser
{
private:
    void RenderSpace(UnsignedIntType& NumberOfAllRenderedAtoms, UnsignedIntType& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix, const Point2fT& MousePositionLocal) override
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
                        TempAtomObject.X = CellEngineSimulationSpace::ConvertToGraphicsCoordinate(SpaceX);
                        TempAtomObject.Y = CellEngineSimulationSpace::ConvertToGraphicsCoordinate(SpaceY);
                        TempAtomObject.Z = CellEngineSimulationSpace::ConvertToGraphicsCoordinate(SpaceZ);
                        if (RenderObject(TempAtomObject, ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                        {
                            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                            for (UnsignedIntType SpaceXP = SpaceX; SpaceXP < SpaceX + 64; SpaceXP += CellEngineConfigDataObject.LoadOfAtomsStep)
                                for (UnsignedIntType SpaceYP = SpaceY; SpaceYP < SpaceY + 64; SpaceYP += CellEngineConfigDataObject.LoadOfAtomsStep)
                                    for (UnsignedIntType SpaceZP = SpaceZ; SpaceZP < SpaceZ + 64; SpaceZP += CellEngineConfigDataObject.LoadOfAtomsStep)
                                        if (CellEngineDataFileObjectPointer->CellEngineSimulationSpaceObjectPointer->Space[SpaceXP][SpaceYP][SpaceZP] != 0)
                                        {
                                            TempAtomObject.X = CellEngineSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceXP));
                                            TempAtomObject.Y = CellEngineSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceYP));
                                            TempAtomObject.Z = CellEngineSimulationSpace::ConvertToGraphicsCoordinate(static_cast<IntType>(SpaceZP));
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
                    }

            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
        }

        DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedVoxelsList);
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
