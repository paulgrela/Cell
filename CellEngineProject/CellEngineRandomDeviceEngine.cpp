
#include <ctime>
#include "CellEngineTypes.h"
#include "CellEngineConfigData.h"
#include "CellEngineRandomDeviceEngine.h"

void CellEngineRandomDeviceEngine::RandomGeneratorSetSeedByRandomDevice()
{
    mt64R.seed(std::random_device{}());
}

void CellEngineRandomDeviceEngine::RandomGeneratorSetSeedByTime()
{
    mt64R.seed(time(nullptr));
}

template <template <class T> class U, class T>
T CellEngineRandomDeviceEngine::GetRandomValue(U<T>& UniformDistributionObject)
{
    T RandomValue = UniformDistributionObject(mt64R);

    if (CellEngineConfigDataObject.SaveGeneratedRandomValuesInVectorBool == true)
        SavedGeneratedRandomValuesVector.emplace_back(RandomValue);

    return RandomValue;
}

template ElectricChargeType CellEngineRandomDeviceEngine::GetRandomValue<std::uniform_int_distribution, ElectricChargeType>(std::uniform_int_distribution<ElectricChargeType>& UniformDistributionObject);
template UnsignedInt CellEngineRandomDeviceEngine::GetRandomValue<std::uniform_int_distribution, UnsignedInt>(std::uniform_int_distribution<UnsignedInt>& UniformDistributionObject);
template SignedInt CellEngineRandomDeviceEngine::GetRandomValue<std::uniform_int_distribution, SignedInt>(std::uniform_int_distribution<SignedInt>& UniformDistributionObject);
template float CellEngineRandomDeviceEngine::GetRandomValue<std::uniform_real_distribution, float>(std::uniform_real_distribution<float>& UniformDistributionObject);

