
#include <cmath>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../CellEngineWellStirredChemicalReactionsSimulation.h"

using ::testing::Return;

void TestWellStirredChemicalReactionsSimulation(CellEngineWellStirredChemicalReactionsSimulation& CellEngineWellStirredChemicalReactionsSimulationObjectObject)
{
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 0, "Water", "H2O", 100 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 1, "Glucose", "C6H12O6", 50 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 2, "Oxygen", "0", 10 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", 5 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 4, "Eten", "CH2CH2", 15 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 25 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", 5 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 7, "HX", "HX", 10 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 10 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 9, "Eten", "CH2CH2", 10 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 10 });
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddParticleKind({ 11, "DNA", "CGATATTAAATAGGGCCT", 10 });

    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddReaction(Reaction("C6H12O6 + O6 + ", { { 1, 1 }, { 2, 6 } }, { { 3, 6 }, { 2, 6 } }));
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddReaction(Reaction("CH2CH2 + H2O + ", { { 4, 1 }, { 0, 1 } }, { { 5, 1 } }));
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddReaction(Reaction("CH3CHCH2 + HX + ", { { 6, 1 }, { 7, 1 } }, { { 8, 1 } }));
    CellEngineWellStirredChemicalReactionsSimulationObjectObject.AddReaction(Reaction("CH2CH2 + O + ", { { 9, 1 }, { 2, 1 } }, { { 10, 1 } }));
}

class MockCellEngineWellStirredChemicalReactionsSimulation : public CellEngineWellStirredChemicalReactionsSimulation
{
public:
    MOCK_METHOD(std::vector<UnsignedInt>, GetRandomParticles, (UnsignedInt));
    //MOCK_METHOD1(GetRandomParticles, std::vector<UnsignedInt>(UnsignedInt));
    MOCK_METHOD0(start, bool());
};

TEST(CellEngineWellStirredChemicalReactionsSimulationTest, TryToDoRandomReactionTest1)
{
    MockCellEngineWellStirredChemicalReactionsSimulation MockCellEngineWellStirredChemicalReactionsSimulationObject;

    EXPECT_CALL(MockCellEngineWellStirredChemicalReactionsSimulationObject, GetRandomParticles(2)).Times(2).WillRepeatedly(Return(vector<UnsignedInt>{ 9, 0 }));

    TestWellStirredChemicalReactionsSimulation(MockCellEngineWellStirredChemicalReactionsSimulationObject);

    MockCellEngineWellStirredChemicalReactionsSimulationObject.Run(2);

    ASSERT_EQ(13, MockCellEngineWellStirredChemicalReactionsSimulationObject.Particles[4].Counter);
    ASSERT_EQ(98, MockCellEngineWellStirredChemicalReactionsSimulationObject.Particles[0].Counter);
    ASSERT_EQ(27, MockCellEngineWellStirredChemicalReactionsSimulationObject.Particles[5].Counter);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}