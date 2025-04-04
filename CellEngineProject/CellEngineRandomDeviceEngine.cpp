
#include <ctime>
#include <string>

#include "DestinationPlatform.h"

#include "CellEngineTypes.h"
#include "CellEngineConfigData.h"
#include "CellEngineRandomDeviceEngine.h"

using namespace std;
void CellEngineRandomDeviceEngine::RandomGeneratorSetSeedByRandomDevice()
{
    mt64R.seed(random_device{}());
}

void CellEngineRandomDeviceEngine::RandomGeneratorSetSeedByTime()
{
    mt64R.seed(time(nullptr));
}

void CellEngineRandomDeviceEngine::SavedGeneratedRandomValuesVectorToFile()
{
    try
    {
        LoggersManagerObject.Log(STREAM("SAVED GENERATED RANDOM VALUES VECTOR SIZE = " << SavedGeneratedRandomValuesVector.size()));

        LoggersManagerObject.Log(STREAM("SAVING OF SAVING GENERATED RANDOM VALUES VECTOR TO BINARY FILE"));

        string SavedGeneratedRandomValuesVectorFileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("binary") + OS_DIR_SEP + string("SavedGeneratedRandomValuesVectorFile.dat");
        ofstream SavedGeneratedRandomValuesVectorFile(SavedGeneratedRandomValuesVectorFileName, ios_base::out | ios_base::trunc | ios_base::binary);

        UnsignedInt SizeToSave = SavedGeneratedRandomValuesVector.size();
        SavedGeneratedRandomValuesVectorFile.write((char*)&SizeToSave, sizeof(SizeToSave));
        for (const auto& SavedGeneratedRandomValueObject : SavedGeneratedRandomValuesVector)
        {
            const SavedGeneratedRandomValue SavedGeneratedRandomValueLocalObject = SavedGeneratedRandomValueObject;
            SavedGeneratedRandomValuesVectorFile.write((char*)&SavedGeneratedRandomValueLocalObject, sizeof(SavedGeneratedRandomValue));
        }

        SavedGeneratedRandomValuesVectorFile.close();

        LoggersManagerObject.Log(STREAM("END OF SAVING GENERATED RANDOM VALUES VECTOR TO BINARY FILE"));
    }
    CATCH("saving generated random values vector to file");
}

template <template <class T> class U, class T>
T CellEngineRandomDeviceEngine::GetRandomValue(U<T>& UniformDistributionObject)
{
    T RandomValue = UniformDistributionObject(mt64R);

    if (CellEngineConfigDataObject.SaveGeneratedRandomValuesInVectorToFileBool == true)
    {
        if constexpr (is_same_v<T, UnsignedInt>)
            SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValue{ static_cast<UnsignedInt>(TypesOfSavedGeneratedRandomValue::UnsignedType), RandomValue, 0, 0 });
        if constexpr (is_same_v<T, SignedInt>)
            SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValue{ static_cast<UnsignedInt>(TypesOfSavedGeneratedRandomValue::SignedType), 0, RandomValue, 0 });
        if constexpr (is_same_v<T, float>)
            SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValue{ static_cast<UnsignedInt>(TypesOfSavedGeneratedRandomValue::FloatType), 0, 0, RandomValue });
    }

    return RandomValue;
}

template ElectricChargeType CellEngineRandomDeviceEngine::GetRandomValue<uniform_int_distribution, ElectricChargeType>(uniform_int_distribution<ElectricChargeType>& UniformDistributionObject);
template UnsignedInt CellEngineRandomDeviceEngine::GetRandomValue<uniform_int_distribution, UnsignedInt>(uniform_int_distribution<UnsignedInt>& UniformDistributionObject);
template SignedInt CellEngineRandomDeviceEngine::GetRandomValue<uniform_int_distribution, SignedInt>(uniform_int_distribution<SignedInt>& UniformDistributionObject);
template float CellEngineRandomDeviceEngine::GetRandomValue<uniform_real_distribution, float>(uniform_real_distribution<float>& UniformDistributionObject);

