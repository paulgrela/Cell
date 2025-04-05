
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

void CellEngineRandomDeviceEngine::WriteSavedGeneratedRandomValuesVectorToFile()
{
    try
    {
        LoggersManagerObject.Log(STREAM("SAVED GENERATED RANDOM VALUES VECTOR SIZE = " << SavedGeneratedRandomValuesVector.size()));

        LoggersManagerObject.Log(STREAM("WRITING SAVED GENERATED RANDOM VALUES VECTOR TO BINARY FILE"));

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
    CATCH("writing saved generated random values vector to file");
}

void CellEngineRandomDeviceEngine::ReadSavedGeneratedRandomValuesVectorFromFile()
{
    try
    {
        LoggersManagerObject.Log(STREAM("READ EARLIER SAVED GENERATED RANDOM VALUES VECTOR FROM BINARY FILE"));

        string SavedGeneratedRandomValuesVectorFileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("binary") + OS_DIR_SEP + string("SavedGeneratedRandomValuesVectorFile.dat");
        ifstream SavedGeneratedRandomValuesVectorFile(SavedGeneratedRandomValuesVectorFileName, ios_base::in | ios_base::binary);

        UnsignedInt SavedGeneratedRandomValuesVectorSize;
        SavedGeneratedRandomValuesVectorFile.read((char*)&SavedGeneratedRandomValuesVectorSize, sizeof(SavedGeneratedRandomValuesVectorSize));
        LoggersManagerObject.Log(STREAM("Number of saved generated random values vector size = " << SavedGeneratedRandomValuesVectorSize));

        for (UnsignedInt GeneObjectIndex = 1; GeneObjectIndex <= SavedGeneratedRandomValuesVectorSize; GeneObjectIndex++)
        {
            SavedGeneratedRandomValue SavedGeneratedRandomValueLocalObject;
            SavedGeneratedRandomValuesVectorFile.read((char*)&SavedGeneratedRandomValueLocalObject, sizeof(SavedGeneratedRandomValueLocalObject));
            SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValueLocalObject);
        }

        SavedGeneratedRandomValuesVectorFile.close();

        LoggersManagerObject.Log(STREAM("END OF READING SAVED GENERATED RANDOM VALUES VECTOR FROM BINARY FILE"));
    }
    CATCH("reading saved generated random values vector to file");
}

template <class T>
void CellEngineRandomDeviceEngine::SaveRandomValue(T RandomValue)
{
    if constexpr (is_same_v<T, UnsignedInt>)
        SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValue{ static_cast<uint16_t>(TypesOfSavedGeneratedRandomValue::UnsignedType), RandomValue, 0, 0 });
    if constexpr (is_same_v<T, SignedInt>)
        SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValue{ static_cast<uint16_t>(TypesOfSavedGeneratedRandomValue::SignedType), 0, RandomValue, 0 });
    if constexpr (is_same_v<T, float>)
        SavedGeneratedRandomValuesVector.emplace_back(SavedGeneratedRandomValue{ static_cast<uint16_t>(TypesOfSavedGeneratedRandomValue::FloatType), 0, 0, RandomValue });
}

template <template <class T> class U, class T>
T CellEngineRandomDeviceEngine::GetRandomValueInside1(U<T>& UniformDistributionObject)
{
    T RandomValue;

    const auto SavedGeneratedRandomValueLocalObject = SavedGeneratedRandomValuesVector[GetRandomValueIndex];
    switch (SavedGeneratedRandomValueLocalObject.Type)
    {
        case static_cast<uint16_t>(TypesOfSavedGeneratedRandomValue::UnsignedType) : RandomValue = SavedGeneratedRandomValueLocalObject.ValueUnsignedInt; break;
        case static_cast<uint16_t>(TypesOfSavedGeneratedRandomValue::SignedType) : RandomValue = SavedGeneratedRandomValueLocalObject.ValueSignedInt; break;
        case static_cast<uint16_t>(TypesOfSavedGeneratedRandomValue::FloatType) : RandomValue = SavedGeneratedRandomValueLocalObject.ValueFloat; break;
        default : break;
    }

    return RandomValue;
}

template <template <class T> class U, class T>
T CellEngineRandomDeviceEngine::GetRandomValueInside2(U<T>& UniformDistributionObject)
{
    T RandomValue;

    if constexpr (is_same_v<T, UnsignedInt>)
        RandomValue = SavedGeneratedRandomValuesVector[GetRandomValueIndex].ValueUnsignedInt;
    if constexpr (is_same_v<T, SignedInt>)
        RandomValue = SavedGeneratedRandomValuesVector[GetRandomValueIndex].ValueSignedInt;
    if constexpr (is_same_v<T, float>)
        RandomValue = SavedGeneratedRandomValuesVector[GetRandomValueIndex].ValueFloat;

    return RandomValue;
}

template <template <class T> class U, class T>
T CellEngineRandomDeviceEngine::GetRandomValue(U<T>& UniformDistributionObject)
{
    T RandomValue;

    if (GetRandomValuesFromSavedGeneratedRandomValuesInVectorBool == false)
    {
        RandomValue = UniformDistributionObject(mt64R);

        if (SaveGeneratedRandomValuesInVectorToFileBool == true)
            SaveRandomValue<T>(RandomValue);
    }
    else
    if (GetRandomValueIndex < SavedGeneratedRandomValuesVector.size())
    {
        RandomValue = GetRandomValueInside1<U, T>(UniformDistributionObject);

        GetRandomValueIndex++;
    }

    return RandomValue;
}

template ElectricChargeType CellEngineRandomDeviceEngine::GetRandomValue<uniform_int_distribution, ElectricChargeType>(uniform_int_distribution<ElectricChargeType>& UniformDistributionObject);
template UnsignedInt CellEngineRandomDeviceEngine::GetRandomValue<uniform_int_distribution, UnsignedInt>(uniform_int_distribution<UnsignedInt>& UniformDistributionObject);
template SignedInt CellEngineRandomDeviceEngine::GetRandomValue<uniform_int_distribution, SignedInt>(uniform_int_distribution<SignedInt>& UniformDistributionObject);
template float CellEngineRandomDeviceEngine::GetRandomValue<uniform_real_distribution, float>(uniform_real_distribution<float>& UniformDistributionObject);

