
#ifndef CELL_ENGINE_AMINO_ACIDS_H
#define CELL_ENGINE_AMINO_ACIDS_H

#include <string>
#include <vector>
#include <unordered_map>

#include "CellEngineTypes.h"

class CellEngineAminoAcid
{
public:
    std::string OneLetterCode{};
    std::string Abbreviation{};
    std::string FullName{};
    std::string tRNAUncharged{};
    std::string tRNACharged{};
    std::string Name{};
    std::vector<std::string> CodonsDNA{};
    std::vector<std::string> CodonsRNA{};
    std::vector<GeneIdInt> GeneID{};
    std::vector<GeneIdInt> SynthaseGeneID{};
    EntityIdInt EntityId;
    EntityIdInt tRNAUnchargedNumId{};
    EntityIdInt tRNAChargedNumId{};
};

const CellEngineAminoAcid Alanine = { "A",	"Ala", 	"Alanine", 	 "M_trnaala_c","M_alatrna_c","M_ala__L_c",{ "GCA", "GCC", "GCG", "GCT" }, { "GCA", "GCC", "GCG", "GCU" }, { 719 }, { 163 }, 0, 0, 0 };

class Asparagine
{
    //B 	Asn	Asparagine 	AAC, AAT, GAC, GAT, M_trnaasp_c,M_asptrna_c,MMSYN1_0724,M_asp__L_c,JCVISYN3A_0287
};

class Cysteine
{
    //C 	Cys 	Cysteine 	TGC, TGT
};

class AsparticAcid
{
    //D 	Asp 	Aspartic acid 	GAC, GAT
};

// A 	Ala 	Alanine 	GCA, GCC, GCG, GCT
// B 	Asx 	Asparagine  AAC, AAT
// C 	Cys 	Cysteine 	TGC, TGT
// D 	Asp 	Aspartic acid 	GAC, GAT

// E 	Glu 	Glutamic acid 	GAA, GAG
// F 	Phe 	Phenylalanine 	TTC, TTT
// G 	Gly 	Glycine 	GGA, GGC, GGG, GGT
// H 	His 	Histidine 	CAC, CAT
// I 	Ile 	Isoleucine 	ATA, ATC, ATT
// K 	Lys 	Lysine 		AAA, AAG
// L 	Leu 	Leucine 	CTA, CTC, CTG, CTT, TTA, TTG
// M 	Met 	Methionine 	ATG
// N 	Asn 	Asparagine 	AAC, AAT
// P 	Pro 	Proline 	CCA, CCC, CCG, CCT
// Q 	Gln 	Glutamine 	CAA, CAG
// R 	Arg 	Arginine 	AGA, AGG, CGA, CGC, CGG, CGT
// S 	Ser 	Serine 	AGC, AGT, TCA, TCC, TCG, TCT
// T 	Thr 	Threonine 	ACA, ACC, ACG, ACT
// V 	Val 	Valine 	GTA, GTC, GTG, GTT
// W 	Trp 	Tryptophan 	TGG
// X 	X 	any codon 	NNN
// Y 	Tyr 	Tyrosine 	TAC, TAT
// Z 	Glx 	Glutamine or Glutamic acid 	CAA, CAG, GAA, GAG


inline std::unordered_map<std::string, const CellEngineAminoAcid*> CellEngineAminoAcids;

void a()
{
    CellEngineAminoAcids["GCA"] = &Alanine;
    CellEngineAminoAcids["GCC"] = &Alanine;
    CellEngineAminoAcids["GCG"] = &Alanine;
    CellEngineAminoAcids["GCT"] = &Alanine;
}

#endif
