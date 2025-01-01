
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include <memory>

#include "CellEngineAtom.h"
#include "CellEngineConfigData.h"
#include "CellEngineFilmOfStructures.h"
#include "CellEngineVoxelSimulationSpace.h"

class CellEngineDataFile : public CellEngineFilmOfStructures
{
public:
    CellEngineDataFile() = default;
public:
    virtual ~CellEngineDataFile() = default;
public:
    std::vector<std::vector<std::vector<std::unique_ptr<CellEngineSimulationSpace>>>> CellEngineSimulationSpaceForThreadsObjectsPointer;
    std::unique_ptr<CellEngineVoxelSimulationSpace> CellEngineVoxelSimulationSpaceObjectPointer;
public:
    bool FilmOfStructuresActive = false;
protected:
    std::unordered_map<UniqueIdInt, Particle> Particles;
protected:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return Particles[ParticleIndex];
    }
protected:
    std::vector<std::vector<CellEngineAtom>> ParticlesCenters;
    std::vector<std::vector<CellEngineAtom>> AllAtoms;
public:
    virtual void ReadDataFromFile(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type) = 0;
public:
    virtual void SaveDataToFile()
    {
    }
public:
    static vmath::vec3 GetCenter(const std::vector<CellEngineAtom>& AtomsParam)
    {
        vmath::vec3 Center(0.0, 0.0, 0.0);

        try
        {
            for (const CellEngineAtom& AtomObject : AtomsParam)
                Center += AtomObject.Position();
            Center /= static_cast<float>(AtomsParam.size());
        }
        CATCH_AND_THROW("counting mass center")

        return Center;
    }
public:
    [[nodiscard]] std::vector<std::vector<CellEngineAtom>>& GetAllAtoms()
    {
        return AllAtoms;
    }
    [[nodiscard]] std::vector<CellEngineAtom>& GetParticlesCenters()
    {
        return ParticlesCenters[CellEngineConfigDataObject.ChosenStructureIndex];
    }
    [[nodiscard]] UnsignedInt GetNumberOfStructures() override
    {
        return ParticlesCenters.size();
    }
};

inline std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;

#endif