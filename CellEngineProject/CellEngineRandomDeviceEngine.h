
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

    template <class T>
    void SaveRandomValue(T RandomValue);

    template <template <class T> class U, class T>
    T GetRandomValueInside1(U<T>& UniformDistributionObject);

    template <template <class T> class U, class T>
    T GetRandomValueInside2(U<T>& UniformDistributionObject);
public:
    void SetGeneratingRandomValuesParameters(const bool SaveGeneratedRandomValuesInVectorToFileBoolParam, const bool GetRandomValuesFromSavedGeneratedRandomValuesInVectorBoolParam)
    {
        SaveGeneratedRandomValuesInVectorToFileBool = SaveGeneratedRandomValuesInVectorToFileBoolParam;
        GetRandomValuesFromSavedGeneratedRandomValuesInVectorBool = GetRandomValuesFromSavedGeneratedRandomValuesInVectorBoolParam;
    }
protected:
    bool SaveGeneratedRandomValuesInVectorToFileBool = false;
    bool GetRandomValuesFromSavedGeneratedRandomValuesInVectorBool = false;
public:
    UnsignedInt GetRandomValueIndex = 0;
public:
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

