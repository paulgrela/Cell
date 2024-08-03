
#include <ctime>
#include "CellEngineRandomDeviceEngine.h"

void CellEngineRandomDeviceEngine::RandomGeneratorSetSeedByRandomDevice()
{
    mt64R.seed(std::random_device{}());
}

void CellEngineRandomDeviceEngine::RandomGeneratorSetSeedByTime()
{
    mt64R.seed(time(nullptr));
}
