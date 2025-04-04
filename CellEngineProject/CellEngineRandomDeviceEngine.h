
#ifndef CELL_ENGINE_RANDOM_DEVICE_ENGINE_H
#define CELL_ENGINE_RANDOM_DEVICE_ENGINE_H

#include <random>

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
    std::vector<std::variant<UnsignedInt, SignedInt, float>> SavedGeneratedRandomValuesVector;
public:
    void ClearSavedGeneratedRandomValuesVector()
    {
        SavedGeneratedRandomValuesVector.clear();
    }
};

#endif

