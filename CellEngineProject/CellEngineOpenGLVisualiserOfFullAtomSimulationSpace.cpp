
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h"

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix)
{
    try
    {
        auto ParticlesCenters = CellEngineDataFileObjectPointer->GetParticlesCenters();

        GLuint PartOfStencilBufferIndex[3];
        vector<pair<UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = 0;
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedAtomsList.clear();

            lock_guard LockGuardObject{ CellEngineDataFileObjectPointer->ChosenStructureMutexObject };

            for (auto ParticlesCenterIterator = ParticlesCenters.begin(); ParticlesCenterIterator != ParticlesCenters.end(); ++ParticlesCenterIterator)
            {
                if (CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject == CellEngineConfigData::StencilForDrawingObjectsTypes::StencilForDrawingOnlyParticlesCenters)
                    glStencilFunc(GL_ALWAYS, uint8_t((NumberOfAllRenderedAtoms >> (8 * StencilBufferLoopCounter))), -1);

                auto ParticlesCenterObject = *ParticlesCenterIterator;

                bool FinalVisibilityInModelWorld = RenderObject(ParticlesCenterObject, ViewMatrix, true, ParticlesCenterIterator == ParticlesCenters.end() - 1, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale);

                if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                    if (FinalVisibilityInModelWorld == true)
                        if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticlesCenterObject.EntityId).Visible == true)
                        {
                            NumberOfFoundParticlesCenterToBeRenderedInAtomDetails++;

                            DrawBonds(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex], BondsBetweenAtomsToDraw[ParticlesCenterObject.AtomIndex], CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);

                            UnsignedInt AtomObjectIndex;
                            for (AtomObjectIndex = 0; AtomObjectIndex < CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex].size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                            {
                                if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                                    if (CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject == CellEngineConfigData::StencilForDrawingObjectsTypes::StencilForDrawingOnlyInAtomScale)
                                    {
                                        glStencilFunc(GL_ALWAYS, uint8_t((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                        TemporaryRenderedAtomsList.emplace_back(ParticlesCenterObject.AtomIndex, AtomObjectIndex);
                                    }
                                RenderObject(CellEngineDataFileObjectPointer->GetAllAtoms()[ParticlesCenterObject.AtomIndex][AtomObjectIndex], ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, false, RenderObjectsBool);
                            }
                        }
            }

            if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
                glReadPixels(GLint(MousePositionLocal.s.X), GLint((float)Info.WindowHeight - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
        }

        if (PressedRightMouseButton != 1)
            DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedAtomsList);
    }
    CATCH("rendering full atom simulation space");
}

inline void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const std::vector<std::pair<UnsignedInt, UnsignedInt>>& TemporaryRenderedAtomsList)
{
    try
    {
        if (CellEngineConfigDataObject.NumberOfStencilBufferLoops > 1)
        {
            UnsignedInt ChosenParticleCenterIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16);

            if (ChosenParticleCenterIndex > 0)
            {
                CellEngineAtom ChosenParticleObject{};
                if (CellEngineConfigDataObject.StencilForDrawingObjectsTypesObject == CellEngineConfigData::StencilForDrawingObjectsTypes::StencilForDrawingOnlyParticlesCenters)
                {
                    if (ChosenParticleCenterIndex > CellEngineDataFileObjectPointer->GetParticlesCenters().size())
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG = " + std::to_string(ChosenParticleCenterIndex) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(CellEngineDataFileObjectPointer->GetParticlesCenters().size()));
                    else
                        ChosenParticleObject = CellEngineDataFileObjectPointer->GetParticlesCenters()[ChosenParticleCenterIndex];
                }
                else
                if (ChosenParticleCenterIndex < TemporaryRenderedAtomsList.size())
                {
                    if (TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first > CellEngineDataFileObjectPointer->GetAllAtoms().size())
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + std::to_string(TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(CellEngineDataFileObjectPointer->GetAllAtoms().size()));
                    else
                    {
                        if (TemporaryRenderedAtomsList[ChosenParticleCenterIndex].second > CellEngineDataFileObjectPointer->GetAllAtoms()[TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first].size())
                            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + std::to_string(TemporaryRenderedAtomsList[ChosenParticleCenterIndex].second) + " MAXIMAL NUMBER OF OBJECTS = " + std::to_string(CellEngineDataFileObjectPointer->GetAllAtoms()[TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first].size()));
                        else
                            ChosenParticleObject = CellEngineDataFileObjectPointer->GetAllAtoms()[TemporaryRenderedAtomsList[ChosenParticleCenterIndex].first][TemporaryRenderedAtomsList[ChosenParticleCenterIndex].second];
                    }
                }

                RenderObject(ChosenParticleObject, ViewMatrix, false, false, false, NumberOfAllRenderedAtoms, true, RenderObjectsBool);

                PrintAtomDescriptionOnScreen(ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using stencil buffer")
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetStartCenterPoint()
{
    Center = CellEngineDataFileObjectPointer->GetCenter(CellEngineDataFileObjectPointer->GetParticlesCenters());
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetMemoryForBondsBetweenAtomsToDraw()
{
    BondsBetweenAtomsToDraw.resize(CellEngineDataFileObjectPointer->GetParticlesCenters().size());
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawBondsForParticlesCenters(std::vector<std::pair<UnsignedInt, UnsignedInt>>& BondsToDraw, const bool DrawBonds, const vmath::mat4& ViewMatrix)
{
    CellEngineOpenGLVisualiser::DrawBonds(CellEngineDataFileObjectPointer->GetParticlesCenters(), BondsToDraw, CellEngineConfigDataObject.DrawBondsBetweenParticlesCenters, ViewMatrix);
}