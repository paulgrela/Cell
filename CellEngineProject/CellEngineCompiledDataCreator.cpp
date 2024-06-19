
#include <string>

#include "CellEngineParticle.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

#include "CellEngineCompiledDataCreator.h"

using namespace std;

void CellEngineCompiledDataCreator::AddParticlesKinds()
{
    try
    {
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", "H2O", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", "C6H12O6", 0, 0, "c", 0  });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen2", "02", "02", 0, 0, "c", 0  });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", "CO2", 0, 0, "c", 0   });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", "CH2CH2", 0, 0, "c", 0   });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", "CH3CH2(OH)", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", "CH3CHCH2", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", "HX", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", "CH3CHXCH3", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", "TEST", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", "CH2CH2O", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", "0", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 12, "Polymerase", "POL", "POL", 0, 0, "c", 0 });


        //MOZE DEFINICJE PONIZSZE JUZ NIEPOTRZEBNE
        ParticlesKindsManagerObject.AddParticleKind({ CellEngineConfigDataObject.DNAIdentifier, "DNA", "DNA", "DNA", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10001, "DNA", "?", "?", 0, 0, "c", 0 });

        LoggersManagerObject.Log(STREAM("ADDED PARTICLES KINDS"));
    }
    CATCH("adding particles kinds and reactions")
};

void CellEngineCompiledDataCreator::AddChemicalReactions()
{
    try
    {
        const string DNASequenceForTestFindingDNA = "TACAAAAAAAGAGGTGTTAGC";

        const string DNASequence1ForTestCutLink1 = "TACAAAAAAAGAGGTGTT";
        const string DNASequence2ForTestCutLink1 = "AGCTCTTATTA";

        const string DNASequence1ForTestCutLink1Any = "ANY";

        const string DNASequence1ForTestCutLink2 = "TACAAAAAAAGAGGTGTT";

        const string DNASequenceForTestCutCrisper = "RNA";
        const string RNASequenceForTestCutCrisper = "ANY";

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1101, "STD ONLY WITH SEQ", "CH3CH2(OH) + DNA + ", {{5, 1, "", true }, {10001, 1, DNASequenceForTestFindingDNA, false } }, {{10, 1, "", true } }, nullptr, "standard normal reaction example"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(10, "CUT 1 SEQ", "CH3CH2(OH) + DNA + ", {{5, 1, "", false }, {10001, 1, DNASequence1ForTestCutLink1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "only in presence of chosen sequence of DNA"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(20, "LINK 1 SEQ", "CH3CH2(OH) + DNA + DNA + ", {{5, 1, "", false }, {10001, 1, DNASequence1ForTestCutLink1, false }, {10001, 1, DNASequence2ForTestCutLink1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInChosenPlaceSpecialReactionFunction, "links one strand of DNA with 2 endings with chosen sequences"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(30, "LINK 1 ANY", "CH3CH2(OH) + DNA + DNA + ", {{5, 1, "", false }, {10001, 2, DNASequence1ForTestCutLink1Any, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "link one strand of DNA with 2 endings with any sequence"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(40, "CUT 2 SEQ SHIFT 3 10", "CH3CHCH2 + DNA + ", 3, 10, {{5, 1, "", false }, {10001, 1, DNASequence1ForTestCutLink2, false } }, {{86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=3 shift=10"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(41, "CUT 2 SEQ SHIFT 7 3", "CH3CHCH2 + DNA + ", 7, 3, {{5, 1, "", false }, {10001, 1, DNASequence1ForTestCutLink2, false } }, {{86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=7 shift2=3"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(42, "CUT 2 SEQ SHIFT 10 3", "CH3CHCH2 + DNA + ", 10, 3, {{5, 1, "", false }, {10001, 1, DNASequence1ForTestCutLink2, false } }, {{86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=10 shift2=3"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(60, "LINK 2 SEQ COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", {{5, 1, "", false }, {10001, 1, "TACAAAAAAAGAGGTGTTAGC", false }, {10001, 1, "TCTTATT", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(61, "LINK 2 SEQ COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", {{5, 1, "", false }, {10001, 1, "TACAAAAAAAGAGGTGTTAGCTCTT", false }, {10001, 1, "ATTATGA", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(70, "LINK 2 ANY COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", {{5, 1, "", false }, {10001, 1, "ANY", false }, {10001, 1, "ANY", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInAnyPlaceSpecialReactionFunction, "Links both DNA strands with any sequence if opposite joining strands when complementary"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(80, "LINK 2 ANY EQU SAME", "CH3CHCH2 + DNA + DNA + ", {{5, 1, "", false }, {10001, 1, "ANY", false }, {10001, 1, "ANY", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "links both strands of DNA if they are cut equally in the same place"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(100, "CUT CRISPER 1", "CH3CHCH2 + RNA + DNA + ", 3, 7, {{5, 1, "", false }, {10001, 1, DNASequenceForTestCutCrisper, false }, {10001, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of one strand of DNA with RNA as template"));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(110, "CUT CRISPER 2", "CH3CHCH2 + RNA + DNA + ", 3, 7, {{5, 1, "", false }, {10001, 1, DNASequenceForTestCutCrisper, false }, {10001, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of two strands of DNA with RNA as template"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(150, "POLYMERASE DNA START SEQ", "CH3CHCH2 + DNA + ", {{5, 1, "", false }, {10001, 1, "TACAAAAAAAGAGGTGTTAGC", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAStartSpecialReactionFunction, "links particle with DNA when found sequence and joins first nucleotide "));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(160, "POLYMERASE DNA CONTINUE", "CH3CHCH2 + DNA + ", {{5, 1, "", false, {CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseRNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));

        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1001, "STD", "C6H12O6 + O2 + ", {{1, 1, "", true }, {2, 6, "", true } }, {{3, 6, "", true }, {0, 6, "", true } }, nullptr));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1002, "STD", "CH2CH2 + H2O + ", {{4, 1, "", true }, {0, 1, "", true } }, {{5, 1, "", true } }, nullptr));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1003, "STD", "CH3CHCH2 + HX + ", {{6, 1, "", true }, {7, 1, "", true } }, {{8, 1, "", true } }, nullptr));
        ChemicalReactionsManagerObject.AddChemicalReaction(ChemicalReaction(1004, "STD", "CH2CH2 + O + ", {{4, 1, "", true }, {11, 1, "", true } }, {{10, 1, "", true } }, nullptr));

        ChemicalReactionsManagerObject.PreprocessChemicalReactions();

        LoggersManagerObject.Log(STREAM("ADDED CHEMICAL REACTIONS"));
    }
    CATCH("adding particles kinds and reactions")
};
