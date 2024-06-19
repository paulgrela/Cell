
#include <cmath>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "ExceptionsMacro.h"
#include "../CellEngineChemicalReactionsManager.h"
#include "../CellEngineWellStirredChemicalReactionsSimulation.h"

using ::testing::Return;

void TestWellStirredChemicalReactionsSimulation(CellEngineWellStirredChemicalReactionsSimulation& CellEngineWellStirredChemicalReactionsSimulationObjectObject)
{
    ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", "H2O", 0, 0, "c", 100 });
    ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", "C6H12O6", 0, 0, "c", 50 });
    ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen", "02", "02", 0, 0, "c", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", "CO2", 0, 0, "c", 5 });
    ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", "CH2CH2", 0, 0, "c", 15 });
    ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", "CH3CH2(OH)", 0, 0, "c", 25 });
    ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", "CH3CHCH2", 0, 0, "c", 5 });
    ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", "HX", 0, 0, "c", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", "CH3CHXCH3", 0, 0, "c", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", "TEST", 0, 0, "c", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", "CH2CH2O", 0, 0, "c", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 11, "DNA", "CGATATTAAATAGGGCCT", "CGATATTAAATAGGGCCT", 0, 0, "c", 10 });

    ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1001, "", "C6H12O6 + O2 + ", {{1, 1, "", true }, {2, 6, "", true} }, {{3, 6, "", true }, {0, 6, "", true } }, nullptr));
    ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1002, "", "CH2CH2 + H2O + ", {{4, 1, "", true }, {0, 1, "", true} }, {{5, 1, "", true } }, nullptr));
    ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1003, "", "CH3CHCH2 + HX + ", {{6, 1, "", true }, {7, 1, "", true} }, {{8, 1, "", true } }, nullptr));
    ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1004, "", "CH2CH2 + O + ", {{4, 1, "", true }, {2, 1, "", true} }, {{10, 1, "", true } }, nullptr));
}

class MockCellEngineWellStirredChemicalReactionsSimulation : public CellEngineWellStirredChemicalReactionsSimulation
{
public:
    MOCK_METHOD(std::vector<UnsignedInt>, GetRandomParticles, (UnsignedInt));
};

TEST(CellEngineWellStirredChemicalReactionsSimulationTest, TryToDoRandomReactionTest1)
{
    MockCellEngineWellStirredChemicalReactionsSimulation MockCellEngineWellStirredChemicalReactionsSimulationObject;

    EXPECT_CALL(MockCellEngineWellStirredChemicalReactionsSimulationObject, GetRandomParticles(2)).Times(2).WillRepeatedly(Return(vector<UnsignedInt>{ 4, 0 }));

    TestWellStirredChemicalReactionsSimulation(MockCellEngineWellStirredChemicalReactionsSimulationObject);

    MockCellEngineWellStirredChemicalReactionsSimulationObject.Run(2);

    ASSERT_EQ(13, ParticlesKindsManagerObject.ParticlesKinds[4].Counter);
    ASSERT_EQ(98, ParticlesKindsManagerObject.ParticlesKinds[0].Counter);
    ASSERT_EQ(27, ParticlesKindsManagerObject.ParticlesKinds[5].Counter);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}