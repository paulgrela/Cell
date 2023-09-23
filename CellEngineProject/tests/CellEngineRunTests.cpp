
#include <cmath>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "ExceptionsMacro.h"
#include "../CellEngineWellStirredChemicalReactionsSimulation.h"

using ::testing::Return;

void TestWellStirredChemicalReactionsSimulation(CellEngineWellStirredChemicalReactionsSimulation& CellEngineWellStirredChemicalReactionsSimulationObjectObject)
{
    ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", 100 });
    ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", 50 });
    ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen", "02", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", 5 });
    ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", 15 });
    ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 25 });
    ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", 5 });
    ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 10 });
    ParticlesKindsManagerObject.AddParticleKind({ 11, "DNA", "CGATATTAAATAGGGCCT", 10 });

    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddChemicalReaction(Reaction(1001, "", "C6H12O6 + O2 + ", { { 1, 1, "", true }, { 2, 6, "", true} }, { { 3, 6, "", true }, { 0, 6, "", true } }));
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddChemicalReaction(Reaction(1002, "", "CH2CH2 + H2O + ", { { 4, 1, "", true }, { 0, 1, "", true} }, { { 5, 1, "", true } }));
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddChemicalReaction(Reaction(1003, "", "CH3CHCH2 + HX + ", { { 6, 1, "", true }, { 7, 1, "", true} }, { { 8, 1, "", true } }));
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddChemicalReaction(Reaction(1004, "", "CH2CH2 + O + ", { { 4, 1, "", true }, { 2, 1, "", true} }, { { 10, 1, "", true } }));
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