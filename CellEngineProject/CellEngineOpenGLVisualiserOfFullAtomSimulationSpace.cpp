
#include <omp.h>
#include <execution>

#include "../Common/Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineOpenGLVisualiserOfFullAtomSimulationSpace.h"

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetStartPositions()
{
    return { SelectionStartXPos, SelectionStartYPos, SelectionStartZPos };
}

std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetSizes()
{
    return { SelectionSizeX, SelectionSizeY, SelectionSizeZ };
}

bool CheckVisibility(const bool Visibility)
{
    if (CellEngineConfigDataObject.TypeOfFileToRead == CellEngineConfigData::TypesOfFileToRead::PDBFile)
        return true;
    else
        return Visibility;
}

//void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace1(const vmath::mat4& ViewMatrix, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms, vector<GPUAtomLocal>& GPUAtomsLocal)
void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace1(const vmath::mat4& ViewMatrix, uint32_t& AtomOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms, vector<GPUAtomLocal>& GPUAtomsLocal)
{
    try
    {
        const auto start_time111 = chrono::high_resolution_clock::now();

        //uint32_t AtomOffsetTotal = 0;
                                                                                                                        uint32_t AtomTotalIndex = 0;//CCC

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles | views::values)
                            if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.EntityId).Visible == true)
                            {
                                GPUParticle GPUParticleObject;

                                GPUParticleObject.AtomOffset = AtomOffsetTotal;
                                GPUParticleObject.AtomCount = ParticleObject.ListOfAtoms.size();

                                GPUParticles.emplace_back(GPUParticleObject);

                                UnsignedInt AtomIndex = 0;
                                for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                {
                                    const auto& AtomObject = ParticleObject.ListOfAtoms[AtomObjectIndex];

                                    if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                    {
                                        GPUAtomLocal GPUAtomLocalObject;
                                        GPUAtomLocalObject.ParticleSectorXIndex = ParticleSectorXIndex;
                                        GPUAtomLocalObject.ParticleSectorYIndex = ParticleSectorYIndex;
                                        GPUAtomLocalObject.ParticleSectorZIndex = ParticleSectorZIndex;
                                        GPUAtomLocalObject.Index = ParticleObject.Index;
                                        GPUAtomLocalObject.AtomOffset = AtomIndex++;
                                        GPUAtomsLocal.emplace_back(GPUAtomLocalObject);
                                    }

                                    GPUAtom GPUAtomObject;

                                    GPUAtomObject.X = AtomObject.X;
                                    GPUAtomObject.Y = AtomObject.Y;
                                    GPUAtomObject.Z = AtomObject.Z;

                                    const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(AtomObject, ParticleObject, ChosenParticleObject.Index == ParticleObject.Index && ChosenAtomObjectIndex == AtomIndex));

                                    GPUAtomObject.ColorR = ParticleColor.X();
                                    GPUAtomObject.ColorG = ParticleColor.Y();
                                    GPUAtomObject.ColorB = ParticleColor.Z();

                                    //GPUAtoms.emplace_back(GPUAtomObject);
                                                                                                                        GPUAtoms[AtomTotalIndex] = std::move(GPUAtomObject);//CCC
                                                                                                                        //GPUAtoms[AtomTotalIndex] = GPUAtomObject;//CCC
                                                                                                                        //memcpy(&GPUAtoms[AtomTotalIndex], &GPUAtomObject, sizeof(GPUAtomObject));//CCC
                                                                                                                        AtomTotalIndex++;//CCC
                                }

                                AtomOffsetTotal += ParticleObject.ListOfAtoms.size();
                            }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);
    }
    CATCH("");
}

//void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace2(const vmath::mat4& ViewMatrix, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms1, vector<GPUAtomLocal>& GPUAtomsLocal)
void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace2(const vmath::mat4& ViewMatrix, uint32_t& AtomOffsetTotal, uint32_t& ParticlesOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms1, vector<GPUAtomLocal>& GPUAtomsLocal)
{
    try
    {
        const auto start_time111 = chrono::high_resolution_clock::now();

        //uint32_t AtomOffsetTotal = 0;

        uint32_t AtomOffsetInSectors[40][40][40];
        vector<GPUParticle> GPUParticlesInSectors[40][40][40];
        vector<GPUAtom> GPUAtomsInSectors[40][40][40];

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = 0;

        omp_set_nested(1);
        omp_set_max_active_levels(2);
        omp_set_dynamic(0);

        #pragma omp parallel for collapse(3) num_threads(128) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, AtomOffsetInSectors, GPUAtomsInSectors, GPUParticlesInSectors, GPUAtomsLocal) schedule(dynamic)
        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale) == true)
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                        {
                            GPUParticle GPUParticleObject;

                            GPUParticleObject.AtomOffset = AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex];
                            GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();

                            GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUParticleObject);

                            UnsignedInt AtomIndex = 0;
                            for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                            {
                                const auto& AtomObject = ParticleObject.second.ListOfAtoms[AtomObjectIndex];

                                if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                {
                                    GPUAtomLocal GPUAtomLocalObject;
                                    GPUAtomLocalObject.ParticleSectorXIndex = ParticleSectorXIndex;
                                    GPUAtomLocalObject.ParticleSectorYIndex = ParticleSectorYIndex;
                                    GPUAtomLocalObject.ParticleSectorZIndex = ParticleSectorZIndex;
                                    GPUAtomLocalObject.Index = ParticleObject.second.Index;
                                    GPUAtomLocalObject.AtomOffset = AtomIndex++;
                                    GPUAtomsLocal.emplace_back(GPUAtomLocalObject);
                                }

                                GPUAtom GPUAtomObject;

                                GPUAtomObject.X = AtomObject.X;
                                GPUAtomObject.Y = AtomObject.Y;
                                GPUAtomObject.Z = AtomObject.Z;

                                const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(AtomObject, ParticleObject.second, ChosenParticleObject.Index == ParticleObject.second.Index && ChosenAtomObjectIndex == AtomIndex));

                                GPUAtomObject.ColorR = ParticleColor.X();
                                GPUAtomObject.ColorG = ParticleColor.Y();
                                GPUAtomObject.ColorB = ParticleColor.Z();

                                GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUAtomObject);
                            }

                            AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] += ParticleObject.second.ListOfAtoms.size();
                        }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);


        const auto start_time114 = chrono::high_resolution_clock::now();

                                                                                                                        UnsignedInt AtomIndex = 0;
        //AtomOffsetTotal = 0;
        //czy moze jak rownolegle wstawia to nie wstawia po kolei i przez to ze sa nie po kolei
        //albo nie wstarcza 128 watkow aby GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].emplace_back(GPUAtomObject); wstawiala niezaleznie dobre

        FOR_EACH_SECTOR_IN_XYZ_ONLY
            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
            {
                AtomOffsetTotal += AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex];
                for (auto& GPUParticleInSectors : GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex])
                {
                    GPUParticleInSectors.AtomOffset += AtomOffsetTotal;
                    GPUParticles.emplace_back(GPUParticleInSectors);
                }
                //GPUAtoms1.insert(GPUAtoms1.end(), make_move_iterator(GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].begin()), make_move_iterator(GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].end()));
                                                                                                                        //copy(GPUAtoms, GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].data(), GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size());
                                                                                                                        memcpy(&GPUAtoms[AtomIndex], GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].data(), GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size() * sizeof(GPUAtom));
                                                                                                                        AtomIndex += GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].size();
            }

        const auto stop_time114 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory3 += chrono::duration(stop_time114 - start_time114);
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace3(const vmath::mat4& ViewMatrix, uint32_t& AtomOffsetTotal, uint32_t& ParticlesOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms1)
{
    try
    {
        const auto start_time114 = chrono::high_resolution_clock::now();

        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt>> ParticlesSectorsToBeRendered;

        //uint32_t AtomOffsetTotal = 0;
        //uint32_t ParticlesOffsetTotal = 0;

        //NIECH SPAMIETUJE POCZATEK I KONIEC KAZDEGO PRZEDZIALU ATOMOW DO KOPIOWANIA
        uint32_t ParticleOffsetInSectors[40][40][40][3000];
        uint32_t AtomOffsetInSectors[40][40][40][3000];
        //uint32_t AtomOffsetInSectors[40][40][40];
        FOR_EACH_SECTOR_IN_XYZ_ONLY
        {
            // ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = 0;
            // AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = 0;

            if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale))
                    if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                    {
                        ParticlesSectorsToBeRendered.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex);

                        UnsignedInt ParticleIndex = 0;
                        for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                        {
                            ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex] = ParticlesOffsetTotal;
                            ParticlesOffsetTotal++;

                            //LoggersManagerObject.Log(STREAM("P=" << ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex]));

                            AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex] = AtomOffsetTotal;
                            //AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] = AtomOffsetTotal;
                            AtomOffsetTotal += ParticleObject.second.ListOfAtoms.size();

                            ParticleIndex++;

                            //LoggersManagerObject.Log(STREAM("A=" << AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] << " S=" << ParticleObject.second.ListOfAtoms.size()));
                        }
                    }
        }

        const auto stop_time114 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory3 += chrono::duration(stop_time114 - start_time114);

        //LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time114, stop_time114, "Execution B has taken time = ","Execution in threads")));


        LoggersManagerObject.Log(STREAM("PARTICLES SIZE = " << ParticlesOffsetTotal << " ATOMS SIZE = " << AtomOffsetTotal));

        const auto start_time111 = chrono::high_resolution_clock::now();

        GPUParticles.resize(ParticlesOffsetTotal);
        //GPUAtoms.resize(AtomOffsetTotal);

         //LoggersManagerObject.Log(STREAM("A=" << ParticlesSectorsToBeRendered.size()));
         //cout << "MAX THREAD " << omp_get_max_threads() << " Thread " << omp_get_thread_num() << " at level " << omp_get_level() << endl;


         // #pragma omp parallel for collapse(3) num_threads(64) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, NumberOfAllRenderedAtoms) schedule(dynamic)
         // FOR_EACH_SECTOR_IN_XYZ_ONLY

         UnsignedInt ParticlesSectorToBeRenderedIndex;
         // omp_set_nested(1);
         // omp_set_max_active_levels(2);
         // omp_set_dynamic(0);
         //#pragma omp parallel for num_threads(64) default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered,NumberOfAllRenderedAtoms) private(ParticlesSectorToBeRenderedIndex) schedule(static)

         //#pragma omp parallel for default(none) shared(CellEngineConfigDataObject, CellEngineDataFileObjectPointer, ViewMatrix, ParticlesKindsManagerObject, RenderObjectsBool, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered,NumberOfAllRenderedAtoms) private(ParticlesSectorToBeRenderedIndex) schedule(static)

         //#pragma omp parallel for num_threads(2) default(none) shared(CellEngineDataFileObjectPointer, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) schedule(dynamic)

         //#pragma omp parallel for default(none) shared(CellEngineDataFileObjectPointer, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) schedule(dynamic)
         //#pragma omp parallel for default(none) shared(CellEngineDataFileObjectPointer, GPUParticles, GPUAtoms, ParticleOffsetInSectors, AtomOffsetInSectors, AtomOffsetTotal, ParticlesSectorsToBeRendered) private(ParticlesSectorToBeRenderedIndex) schedule(dynamic)
         //#pragma omp parallel for

         // size_t totalSize = ParticlesSectorsToBeRendered.size();
         // size_t halfSize = (totalSize + 1) / 2;
         // #pragma omp parallel for num_threads(2) schedule(static, halfSize)



         #pragma omp parallel for num_threads(4) schedule(static)
         for (ParticlesSectorToBeRenderedIndex = 0; ParticlesSectorToBeRenderedIndex < ParticlesSectorsToBeRendered.size(); ParticlesSectorToBeRenderedIndex++)
         {
             const UnsignedInt ParticleSectorXIndex = get<0>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);
             const UnsignedInt ParticleSectorYIndex = get<1>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);
             const UnsignedInt ParticleSectorZIndex = get<2>(ParticlesSectorsToBeRendered[ParticlesSectorToBeRenderedIndex]);

                 //cout << "MAX THREAD S " << omp_get_max_threads() << " Thread S " << omp_get_thread_num() << " at level S " << omp_get_level() << endl;

                 // if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                 //     if (RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, NumberOfAllRenderedAtoms, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale, ParticleSectorXIndex))
                 //         if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                 UnsignedInt ParticleIndex = 0;
                 for (const auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                 {
                     GPUParticle GPUParticleObject;

                     //GPUParticleObject.AtomOffset = AtomOffsetTotal;
                     //TU JEDEN ATOM_OFFSET DLA CZASTKI ATOM - A POWINIEN BYC JESZCZE ZAPAMIETANY DLA KAZDEJ CZASTKI Z OSOBNA
                     //BO W JEDNYM SEKTORZE MUSI BYC ZAPAMIETANYCH WIELE PRZESUNIEC DLA ATOMOW
                     GPUParticleObject.AtomOffset = AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex];
                     GPUParticleObject.AtomCount = ParticleObject.second.ListOfAtoms.size();

                     //LoggersManagerObject.Log(STREAM("P=" << ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex]));

                     //GPUParticles[ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] + ParticleIndex] = GPUParticleObject;
                     GPUParticles[ParticleOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex]] = GPUParticleObject;

                     UnsignedInt AtomIndex = 0;
                     for (const auto& Atom : ParticleObject.second.ListOfAtoms)
                     {
                         GPUAtom GPUAtomObject;

                         GPUAtomObject.X = Atom.X;
                         GPUAtomObject.Y = Atom.Y;
                         GPUAtomObject.Z = Atom.Z;

                         const auto ParticleColor = CellEngineUseful::GetVMathVec3FromVector3ForColor(GetColor<CellEngineAtom>(Atom, ParticleObject.second, ChosenParticleObject.Index == ParticleObject.second.Index && ChosenAtomObjectIndex == AtomIndex));

                         GPUAtomObject.ColorR = ParticleColor.X();
                         GPUAtomObject.ColorG = ParticleColor.Y();
                         GPUAtomObject.ColorB = ParticleColor.Z();

                         //LoggersManagerObject.Log(STREAM("A=" << AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] + AtomIndex << " S=" << ParticleObject.second.ListOfAtoms.size()));

                         GPUAtoms[AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex][ParticleIndex] + AtomIndex] = GPUAtomObject;
                         //GPUAtoms[AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex] + AtomIndex] = GPUAtomObject;
                         AtomIndex++;
                     }
                     ParticleIndex++;
                 }
         }

        const auto stop_time111 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory0 += chrono::duration(stop_time111 - start_time111);

        // LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time111, stop_time111, "Execution Z has taken time = ","Execution in threads")));
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace0(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix)
{
    try
    {
        GLuint PartOfStencilBufferIndex[3];
        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        for (UnsignedInt StencilBufferLoopCounter = 0; StencilBufferLoopCounter < CellEngineConfigDataObject.NumberOfStencilBufferLoops; StencilBufferLoopCounter++)
        {
            NumberOfAllRenderedAtoms = 0;

            TemporaryRenderedAtomsList.clear();

            lock_guard LockGuardObject{ CellEngineDataFile::ChosenStructureMutexObject };

            FOR_EACH_SECTOR_IN_XYZ_ONLY
            {
                if (CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.empty() == false)
                {
                    bool FinalVisibilityInModelWorld = RenderObject(CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.begin()->second.ListOfAtoms.back(), Particle(), ViewMatrix, true, false, true, false, !CellEngineConfigDataObject.ShowDetailsInAtomScale);

                    FinalVisibilityInModelWorld = CheckVisibility(FinalVisibilityInModelWorld);

                    if (FinalVisibilityInModelWorld == true)
                        if (CellEngineConfigDataObject.ShowDetailsInAtomScale == true)
                            for (auto& ParticleObject : CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles)
                                if (ParticlesKindsManagerObject.GetGraphicParticleKind(ParticleObject.second.EntityId).Visible == true)
                                {
                                    DrawBonds(ParticleObject.second, ParticleObject.second.BondsBetweenAtomsToDraw, CellEngineConfigDataObject.DrawBondsBetweenAtoms, ViewMatrix);

                                    ParticleObject.second.ParticleSectorPos = SectorPosType{ static_cast<SignedInt>(ParticleSectorXIndex), static_cast<SignedInt>(ParticleSectorYIndex), static_cast<SignedInt>(ParticleSectorZIndex) };

                                    for (UnsignedInt AtomObjectIndex = 0; AtomObjectIndex < ParticleObject.second.ListOfAtoms.size(); AtomObjectIndex += CellEngineConfigDataObject.LoadOfAtomsStep)
                                    {
                                        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                        {
                                            glStencilFunc(GL_ALWAYS, static_cast<uint8_t>((TemporaryRenderedAtomsList.size()) >> (8 * StencilBufferLoopCounter)), -1);
                                            TemporaryRenderedAtomsList.emplace_back(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, ParticleObject.first, AtomObjectIndex);
                                        }

                                        RenderObject(ParticleObject.second.ListOfAtoms[AtomObjectIndex], ParticleObject.second, ViewMatrix, false, false, false, false, RenderObjectsBool);
                                        NumberOfAllRenderedAtoms++;
                                    }
                                }
                }
            }

            if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                glReadPixels(static_cast<GLint>(MousePositionLocal.s.X), static_cast<GLint>(static_cast<float>(Info.WindowHeight) - MousePositionLocal.s.Y - 1), 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &PartOfStencilBufferIndex[StencilBufferLoopCounter]);
        }

        if (PressedRightMouseButton != 1)
            DrawChosenAtomUsingStencilBuffer(ViewMatrix, PartOfStencilBufferIndex, NumberOfAllRenderedAtoms, TemporaryRenderedAtomsList);
    }
    CATCH("");
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix)
{
    try
    {
        //RenderSpace0(ViewMatrix);

        vector<GPUParticle> GPUParticles;
        GPUParticles.reserve(10'000'000);
        vector<GPUAtom> GPUAtoms1;
                                                                                                                        //GPUAtoms1.reserve(1000'000'000);
        vector<GPUAtomLocal> GPUAtomsLocal;
        GPUAtomsLocal.reserve(1000'000'000);

        // if (CellEngineConfigDataObject.ViewPositionZ < 2300)
        //     RenderSpace1(ViewMatrix, GPUParticles, GPUAtoms, GPUAtomsLocal);
        // else
        //   RenderSpace2(ViewMatrix, GPUParticles, GPUAtoms, GPUAtomsLocal);
        //   RenderSpace3(ViewMatrix, GPUParticles, GPUAtoms1);
        uint32_t AtomOffsetTotal = 0;//CCC
        uint32_t ParticlesOffsetTotal = 0;//CCC

            //RenderSpace1(ViewMatrix, AtomOffsetTotal, GPUParticles, GPUAtoms, GPUAtomsLocal);//CCC
            RenderSpace2(ViewMatrix, AtomOffsetTotal, ParticlesOffsetTotal, GPUParticles, GPUAtoms1, GPUAtomsLocal);//CCC
            //RenderSpace3(ViewMatrix, AtomOffsetTotal, ParticlesOffsetTotal, GPUParticles, GPUAtoms1);//CCC

        //LoggersManagerObject.Log(STREAM("P=" << GPUParticles.size() << " " << GPUAtoms.size()));

        //NumberOfAllRenderedAtoms = GPUAtoms1.size();
        NumberOfAllRenderedAtoms = AtomOffsetTotal;//CCC
        NumberOfFoundParticlesCenterToBeRenderedInAtomDetails = GPUParticles.size();


        glUseProgram(ComputeShaderProgramPhong);

        glUniform3fv(glGetUniformLocation(ComputeShaderProgramPhong, "Center"), 1, Center);
        glUniformMatrix4fv(glGetUniformLocation(ComputeShaderProgramPhong, "ViewMatrix"), 1, GL_FALSE, ViewMatrix);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ParticleSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, AtomSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

        const auto start_time112 = chrono::high_resolution_clock::now();

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ParticleSSBO);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUParticles.size() * sizeof(GPUParticle), GPUParticles.data());

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, AtomSSBO);
        //glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, GPUAtoms1.size() * sizeof(GPUAtom), GPUAtoms1.data());
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (AtomOffsetTotal - 1) * sizeof(GPUAtom), GPUAtoms.data());//CCC

        const auto stop_time112 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory1 += chrono::duration(stop_time112 - start_time112);

        const auto start_time113 = chrono::high_resolution_clock::now();

        glDispatchCompute((GPUParticles.size() + 255) / 256, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        const auto stop_time113 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory2 += chrono::duration(stop_time113 - start_time113);

                                                                                                                        glBindFramebuffer(GL_FRAMEBUFFER, FrameBufferObject);


        const vmath::vec3 BackgroundColor = CellEngineConfigDataObject.BackgroundColors[CellEngineConfigDataObject.ChosenBackgroundColor];
        glClearColor(BackgroundColor.data[0], BackgroundColor.data[1], BackgroundColor.data[2], 0.0f);

                                                                                                                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                                                                                                                        constexpr GLuint ClearValue = 0xFFFFFFFF;
                                                                                                                        glClearTexImage(ScreenBufferInstanceTexture, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &ClearValue);

                                                                                                                        glEnable(GL_DEPTH_TEST);
                                                                                                                        glDepthFunc(GL_LESS);

        glUseProgram(ShaderProgramPhong);

        glUniformMatrix4fv(glGetUniformLocation(ShaderProgramPhong, "ProjectionMatrix"), 1, GL_FALSE, ProjectionMatrixGlobal);

                                                                                                                        // glUniform3fv(glGetUniformLocation(ShaderProgramPhong, "cameraPos"), 1, &CellEngineConfigDataObject.ViewPositionZ);//???
                                                                                                                        // glUniform1f(glGetUniformLocation(ShaderProgramPhong, "billboardDistance"), 2300.0f);//???
                                                                                                                        // glEnable(GL_PROGRAM_POINT_SIZE);//???
                                                                                                                        // glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);//???



        const auto start_time1 = chrono::high_resolution_clock::now();

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ParticlesAtomsBufferSharedBetweenComputeShaderAndVertexShaderSSBO);

        vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>> TemporaryRenderedAtomsList;

        //AtomGraphicsObject.RenderSubGraphicObject(0, GPUAtoms1.size(), 0);
        AtomGraphicsObject.RenderSubGraphicObject(0, AtomOffsetTotal, 0);//CCC

                                                                                                                        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

                                                                                                                        glReadBuffer(GL_COLOR_ATTACHMENT0);

                                                                                                                        glBlitFramebuffer(0, 0, Info.WindowWidth, Info.WindowHeight, 0, 0, Info.WindowWidth, Info.WindowHeight, GL_COLOR_BUFFER_BIT,GL_LINEAR);

        const auto stop_time1 = chrono::high_resolution_clock::now();

        ExecutionDurationTimeForCopyingParticlesToGraphicMemory += chrono::duration(stop_time1 - start_time1);




                                                                                                                        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
                                                                                                                        {
                                                                                                                            uint32_t ClickedObjectID = 0xFFFFFFFF;

                                                                                                                            glBindFramebuffer(GL_READ_FRAMEBUFFER, FrameBufferObject);
                                                                                                                            glReadBuffer(GL_COLOR_ATTACHMENT1);

                                                                                                                            glReadPixels(static_cast<GLint>(MousePositionLocal.s.X), static_cast<GLint>(static_cast<float>(Info.WindowHeight) - MousePositionLocal.s.Y - 1), 1, 1, GL_RED_INTEGER, GL_UNSIGNED_INT, &ClickedObjectID);

                                                                                                                            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

                                                                                                                            if (PressedRightMouseButton != 1)
                                                                                                                                DrawChosenAtomUsingStencilBuffer1(ClickedObjectID, GPUAtomsLocal);

                                                                                                                            //LoggersManagerObject.Log(STREAM("C=" << CellEngineConfigDataObject.NumberOfStencilBufferLoops << " " << MousePositionLocal.s.X << " " << MousePositionLocal.s.Y << " " << Info.WindowWidth << " " << Info.WindowHeight << " " << ClickedObjectID));
                                                                                                                        }

    }
    CATCH("rendering full atom simulation space");
}

inline void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawChosenAtomUsingStencilBuffer1(const GLuint ChosenParticleCenterIndex, const vector<GPUAtomLocal>& GPUAtomsLocal)
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
                    {
                        if (GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset > ParticleIter->second.ListOfAtoms.size())
                            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + std::to_string(ChosenParticleCenterIndex));
                        else
                        {
                            ChosenParticleObject = ParticleIter->second;
                            ChosenAtomObject = ParticleIter->second.ListOfAtoms[GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset];
                            ChosenAtomObjectIndex = GPUAtomsLocal[ChosenParticleCenterIndex].AtomOffset + 1;
                        }
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

inline void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const vector<tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>>& TemporaryRenderedAtomsList)
{
    try
    {
        if (CellEngineConfigDataObject.ShowDetailsOfPickedAtomParticle == true)
        {
            if (const UnsignedInt ChosenParticleCenterIndex = PartOfStencilBufferIndex[0] | (PartOfStencilBufferIndex[1] << 8) | (PartOfStencilBufferIndex[2] << 16); ChosenParticleCenterIndex > 0)
            {
                if (ChosenParticleCenterIndex < TemporaryRenderedAtomsList.size())
                {
                    const UnsignedInt ParticleSectorXIndex = get<0>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]);
                    const UnsignedInt ParticleSectorYIndex = get<1>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]);
                    const UnsignedInt ParticleSectorZIndex = get<2>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]);

                    if (const auto ParticleIter = CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.find(get<3>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])); ParticleIter != CellEngineDataFileObjectPointer->GetParticles()[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].Particles.end())
                    {
                        if (get<4>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex]) > ParticleIter->second.ListOfAtoms.size())
                            throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 2 = " + std::to_string(get<4>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])));
                        else
                        {
                            ChosenParticleObject = ParticleIter->second;
                            ChosenAtomObject = ParticleIter->second.ListOfAtoms[get<4>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])];
                        }
                    }
                    else
                        throw std::runtime_error("ERROR STENCIL INDEX TOO BIG IN INNER 1 = " + std::to_string(get<3>(TemporaryRenderedAtomsList[ChosenParticleCenterIndex])));
                }

                RenderObject(ChosenAtomObject, ChosenParticleObject, ViewMatrix, false, false, false, true, RenderObjectsBool);

                PrintAtomDescriptionOnScreen(ChosenAtomObject, ChosenParticleObject);
            }
        }
    }
    CATCH("choosing atom using stencil buffer")
}

void CellEngineOpenGLVisualiserOfFullAtomSimulationSpace::GetStartCenterPoint()
{
    Center = CellEngineDataFileObjectPointer->GetCenterForAllParticles();

    LoggersManagerObject.Log(STREAM("CENTER OF CELL = " << Center.X() << " " << Center.Y() << " " << Center.Z()));
}

#endif