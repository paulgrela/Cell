
#ifndef CELL_ENGINE_RANDOM_DEVICE_ENGINE_H
#define CELL_ENGINE_RANDOM_DEVICE_ENGINE_H

#include <random>

struct SavedGeneratedRandomValue
{
    std::uint16_t Type{ 0 };
    UnsignedInt ValueUnsignedInt{ 0 };
    SignedInt ValueSignedInt{ 0 };
    float ValueFloat{ 0 };
};

class CellEngineRandomDeviceEngine
{
protected:
    std::mt19937_64 mt64R{ std::random_device{}() };
public:
    void RandomGeneratorSetSeedByRandomDevice();
    void RandomGeneratorSetSeedByTime();
public:
    template <template <class T> class U, class T>
    T GetRandomValue(U<T>& UniformDistributionObject);
public:
    bool SaveGeneratedRandomValuesInVectorToFileBool = false;
    bool GetRandomValuesFromSavedGeneratedRandomValuesInVector = false;
public:
    UnsignedInt GetRandomValueIndex = 0;
protected:
    std::vector<SavedGeneratedRandomValue> SavedGeneratedRandomValuesVector;
public:
    void ClearSavedGeneratedRandomValuesVector()
    {
        SavedGeneratedRandomValuesVector.clear();
    }
    void WriteSavedGeneratedRandomValuesVectorToFile();
    void ReadSavedGeneratedRandomValuesVectorFromFile();
};

#endif

