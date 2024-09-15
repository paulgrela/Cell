
#include <string>

#include "CellEngineParticle.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

#include "CellEngineCompiledDataCreator.h"

using namespace std;

void CellEngineCompiledDataCreator::AddSpecialParticlesKinds()
{
    try
    {
        EntityIdInt LocalRNAIdentifier = CellEngineConfigDataObject.RNAIdentifier;
        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::DNANucleotide, LocalRNAIdentifier, "RNANucleotide", "RNANucleotide", "rna", -1, 0, "c", 1);

        EntityIdInt LocalProteinInBuildingProcessIdentifier = CellEngineConfigDataObject.ProteinInBuildingProcessIdentifier;
        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::ProteinInBuildingProcess, LocalProteinInBuildingProcessIdentifier, "ProteinInBuildingProcess", "ProteinInBuildingProcess", "ProteinInBuildingProcess", -1, 0, "c", 1);

        LoggersManagerObject.Log(STREAM("ADDED SPECIAL PARTICLES KINDS"));
    }
    CATCH("adding special particles kinds")
};

void CellEngineCompiledDataCreator::AddParticlesKinds()
{
    try
    {
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", "H2O", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", "C6H12O6", 0, 0, "c", 0  });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen2", "02", "02", 0, 0, "c", 0  });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "CarbonDioxide", "CO2", "CO2", 0, 0, "c", 0   });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", "CH2CH2", 0, 0, "c", 0   });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", "CH3CH2(OH)", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", "CH3CHCH2", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", "HX", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2HalogenoPropan", "CH3CHXCH3", "CH3CHXCH3", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", "CH2CH2O", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", "0", 0, 0, "c", 0 });

        LoggersManagerObject.Log(STREAM("ADDED PARTICLES KINDS"));
    }
    CATCH("adding particles kinds")
};

void CellEngineCompiledDataCreator::AddChemicalReactions()
{
    try
    {
        const string StartSeq1 = "AAAWWTWTTTNNNAAANNNNNTTGACANNNNNNNNNNNNTGTGNTATAATNNNNNNANNNNNNNNNNN";
        const string StartSeq2 = "AAATATCTAACACCGTGCGTGTTGACTATTTTACCTCTGGCGGTGATAATGGTTGCATGTACTAAGGA";
        const string StartSeq3 = "TAAAATTTATCAAAAAGAGTATTGACTTAAAGTCTAACCTATAGGATACTTACAGCCATCGAGAGGGAC";
        const string StartSeq4 = "AAAATTATTTTAAATTTCCTCTTGACTAGGCCGGATAACTCCCTATAATGCGCCACCACTGACACGGAA";

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(151, "POLYMERASE RNA START SEQ", "DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, StartSeq1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links rna polymerase particle with DNA when found sequence and joins first nucleotide "));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(152, "POLYMERASE RNA START SEQ", "DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, StartSeq2, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links rna polymerase particle with DNA when found sequence and joins first nucleotide "));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(153, "POLYMERASE RNA START SEQ", "DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, StartSeq3, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links rna polymerase particle with DNA when found sequence and joins first nucleotide "));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(154, "POLYMERASE RNA START SEQ", "DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, StartSeq4, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links rna polymerase particle with DNA when found sequence and joins first nucleotide "));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(161, "POLYMERASE RNA CONTINUE", "ATP+DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("m_ppi_c")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(162, "POLYMERASE RNA CONTINUE", "CTP+DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("m_ppi_c")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(163, "POLYMERASE RNA CONTINUE", "GTP+DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("m_ppi_c")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(164, "POLYMERASE RNA CONTINUE", "TTP+DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("m_ppi_c")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(165, "POLYMERASE RNA CONTINUE", "UTP+DNA+rnapol+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RNAPolymerase")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("m_ppi_c")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));

        const string RNAStartSeq = "ACGTA";
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(171, "RIBOSOME RNA START SEQ", "GTP+R70S+RNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RIBOSOME")->EntityId, 1, "", false }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("M_gtp_c")->EntityId, 1, "", true }, { CellEngineConfigDataObject.RNAIdentifier, 1, RNAStartSeq, false } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("ProteinInBuildingProcess")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links ribosome particle with RNA when found sequence and joins first nucleotide "));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(181, "RIBOSOME RNA CONTINUE", "M_leutrna_c+R70S+RNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("particle_RIBOSOME")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.ProteinInBuildingProcessIdentifier } }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("M_leutrna_c")->EntityId, 1, "", true }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("M_gtp_c")->EntityId, 1, "", true } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("M_trnaleu_c")->EntityId, 1, "", false } }, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));

        ChemicalReactionsManagerObject.PreprocessChemicalReactions();

        LoggersManagerObject.Log(STREAM("ADDED CHEMICAL REACTIONS"));
    }
    CATCH("adding chemical reactions")
};

void CellEngineCompiledDataCreator::AddTestChemicalReactions()
{
    try
    {
        const string DNASequenceForTestFindingDNA = "CTCATTT";

        const string DNASequence1ForTestCutLink1 = "CTCATTTTTT";
        const string DNASequence2ForTestCutLink1 = "TTTCTA";
        const string DNASequence1ForTestCutLink1Any = "ANY";
        const string DNASequence1ForTestCutLink2 = "CTCATTT";

        const string DNASequenceForTestCutCrisper = "RNA";
        const string RNASequenceForTestCutCrisper = "ANY";

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1101, "STD ONLY WITH SEQ", "Ethanol+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", true }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequenceForTestFindingDNA, false } }, { { 10, 1, "", true } }, nullptr, "standard normal reaction example"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(10, "CUT 1 SEQ", "Ethanol+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequence1ForTestCutLink1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "only in presence of chosen sequence of DNA"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(20, "LINK 1 SEQ", "Ethanol+DNA+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequence1ForTestCutLink1, false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequence2ForTestCutLink1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInChosenPlaceSpecialReactionFunction, "links one strand of DNA with 2 endings with chosen sequences"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(30, "LINK 1 ANY", "Ethanol+DNA+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 2, DNASequence1ForTestCutLink1Any, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "link one strand of DNA with 2 endings with any sequence"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(40, "CUT 2 SEQ SHIFT 3 10", "Propen+DNA+", 2, 5, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=3 shift=10"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(41, "CUT 2 SEQ SHIFT 7 3", "Propen+DNA+", 3, 1, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=7 shift2=3"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(42, "CUT 2 SEQ SHIFT 10 3", "Propen+DNA+", 5, 2, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=10 shift2=3"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(60, "LINK 2 SEQ COMPLEMENT", "Propen+DNA+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "CTCATTT", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "TTTCT", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(61, "LINK 2 SEQ COMPLEMENT", "Propen+DNA+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "CTCATTTTTT", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "CTAAAT", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(70, "LINK 2 ANY COMPLEMENT", "Propen+DNA+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "ANY", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "ANY", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInAnyPlaceSpecialReactionFunction, "Links both DNA strands with any sequence if opposite joining strands when complementary"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(80, "LINK 2 ANY EQU SAME", "Propen+DNA+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "ANY", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "ANY", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "links both strands of DNA if they are cut equally in the same place"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(100, "CUT CRISPER 1", "Propen+RNA+DNA+", 3, 7, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequenceForTestCutCrisper, false }, { CellEngineConfigDataObject.DNAIdentifier, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of one strand of DNA with RNA as template"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(110, "CUT CRISPER 2", "Propen+RNA+DNA+", 3, 7, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, DNASequenceForTestCutCrisper, false }, { CellEngineConfigDataObject.DNAIdentifier, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of two strands of DNA with RNA as template"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(150, "POLYMERASE RNA START SEQ", "Ethanol+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", false }, { CellEngineConfigDataObject.DNAIdentifier, 1, "CTCATTT", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links particle with DNA when found sequence and joins first nucleotide "));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(160, "POLYMERASE RNA CONTINUE", "Ethanol+DNA+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1001, "STD", "Glucose+Oxygen2+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Glucose")->EntityId, 1, "", true }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("Oxygen2")->EntityId, 6, "", true } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("CarbonDioxide")->EntityId, 6, "", true }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("Water")->EntityId, 6, "", true } }, nullptr));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1002, "STD", "Eten+Water+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Eten")->EntityId, 1, "", true }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("Water")->EntityId, 1, "", true } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethanol")->EntityId, 1, "", true } }, nullptr));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1003, "STD", "Propen+HX+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Propen")->EntityId, 1, "", true }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("2HalogenPropan")->EntityId, 1, "", true } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("2HalogenPropan")->EntityId, 1, "", true } }, nullptr));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1004, "STD", "Eten+Oxygen+", { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Eten")->EntityId, 1, "", true }, { ParticlesKindsManagerObject.GetParticleKindFromStrId("Oxygen")->EntityId, 1, "", true } }, { { ParticlesKindsManagerObject.GetParticleKindFromStrId("Ethylene")->EntityId, 1, "", true } }, nullptr));

        ChemicalReactionsManagerObject.PreprocessChemicalReactions();

        LoggersManagerObject.Log(STREAM("ADDED TEST CHEMICAL REACTIONS"));
    }
    CATCH("adding chemical reactions")
};
