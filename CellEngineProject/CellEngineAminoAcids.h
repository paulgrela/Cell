
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

const CellEngineAminoAcid Arginine = { "R", "Arg", "Arginine", "M_trnaarg_c", "M_argtrna_c",	"M_arg__L_c", { "AGA", "AGG", "CGA", "CGC", "CGG", "CGT" }, { "AGA", "AGG", "CGA", "CGC", "CGG", "CGU" }, {374, 717 }, { 535 }, 0, 0, 0 };

const CellEngineAminoAcid Asparagine = { "N", "Asn", "Asparagine", "M_trnaasn_c", "M_asntrna_c", "M_asn__L_c",{ "AAC", "AAT" }, { "AAC", "AAU" }, { 681 }, { 76 }, 0, 0, 0 };

const CellEngineAminoAcid AsparticAcid = { "D", "Asp", "AsparticAcid", "M_trnaasp_c" ,"M_asptrna_c", "M_asp__L_c", { "GAC", "GAT" }, { "GAC", "GAU" }, { 724 }, { 287 }, 0, 0, 0 };

const CellEngineAminoAcid Cysteine = { "C", 	"Cys", 	"Cysteine",  "M_trnacys_c", "M_cystrna_c","M_cys__L_c" ,{ "TGC", "TGT" }, { "UGC", "UGU" }, { 828 }, { 837 }, 0, 0, 0 };

const CellEngineAminoAcid Glutamine = { "Q", "Gln", "Glutamine", "M_trnagln_c", "M_glntrna_c", "M_gln__L_c", { "CAA", "CAG" }, { "CAA", "CAG" }, { 508 }, { 1 }, 0, 0, 0 };

const CellEngineAminoAcid GlutamicAcid = { "E", "Glu", "GlutamicAcid", "M_trnaglu_c", "M_glutrna_c", "M_glu__L_c" ,{ "GAA", "GAG" }, { "GAA", "GAG" }, { 680 }, {126 }, 0, 0, 0 };

const CellEngineAminoAcid Glycine = { "G", "Gly", "Glycine", "M_trnagly_c", "M_glytrna_c", "M_gly_c", { "GGA", "GGC", "GGG", "GGT" }, { "GGA", "GGC", "GGG", "GGU" }, { 295 }, { 405 }, 0, 0, 0 };

const CellEngineAminoAcid Histidine = { "H", "His", "Histidine", "M_trnahis_c", "M_histrna_c", "M_his__L_c", { "CAC", "CAT" }, { "CAC", "CAU" }, { 624 }, { 288 }, 0, 0, 0 };

const CellEngineAminoAcid Isoleucine = { "I", 	"Ile", "Isoleucine", "M_trnaile_c", "M_iletrna_c", "M_ile__L_c",	{ "ATA", "ATC", "ATT" }, { "AUA", "AUC", "AUU" }, { 635 }, { 519 }, 0, 0, 0 };

const CellEngineAminoAcid Leucine = { "L", "Leu", "Leucine", "M_trnaleu_c", "M_leutrna_c", "M_leu__L_c", { "CTA", "CTC", "CTG", "CTT", "TTA", "TTG" }, { "CUA", "CUC", "CUG", "CUU", "UUA", "UUG" }, { 70, 423, 506 }, { 634 }, 0, 0, 0 };

const CellEngineAminoAcid Lysine = { "K", "Lys", "Lysine", "M_trnalys_c", "M_lystrna_c", "M_lys__L_c", { "AAA", "AAG" }, { "AAA", "AAG" }, { 71, 507 }, { 64 }, 0, 0, 0 };

const CellEngineAminoAcid Methionine = { "M", "Met", "Methionine", "M_trnamet_c", "M_mettrna_c", "M_met__L_c", { "ATG" }, { "AUG" }, { 720, 721, 723 }, { 12 }, 0, 0, 0, };

const CellEngineAminoAcid Phenylalanine = { "F", "Phe", "Phenylalanine", "M_trnaphe_c", "M_phetrna_c", "M_phe__L_c", { "TTC", "TTT" }, { "UUC", "UUU" }, { 725 }, { 529 }, 0, 0, 0 };

const CellEngineAminoAcid Proline = { "P", "Pro", "Proline", "M_trnapro_c", "M_protrna_c", "M_pro__L_c", { "CCA", "CCC", "CCG", "CCT" }, { "CCA", "CCC", "CCG", "CCU" }, { 718 }, { 282 }, 0, 0, 0 };

const CellEngineAminoAcid Serine = { "S", "Ser", "Serine", "M_trnaser_c", "M_sertrna_c", "M_ser__L_c", { "AGC", "AGT", "TCA", "TCC", "TCG", "TCT" }, { "AGC", "AGU", "UCA", "UCC", "UCG", "UCU" }, { 280, 722 }, { 61 }, 0, 0, 0 };

const CellEngineAminoAcid Threonine = { "T", "Thr", "Threonine", "M_trnathr_c", "M_thrtrna_c", "M_thr__L_c", { "ACA", "ACC", "ACG", "ACT" }, { "ACA", "ACC", "ACG", "ACU" }, { 510, 678 }, { 222 }, 0, 0, 0 };

const CellEngineAminoAcid Tryptophan = { "W", "Trp", "Tryptophan", "M_trnatrp_c", "M_trptrna_c", "M_trp__L_c", { "TGG" }, { "UGG" }, { 618, 619 }, { 308 }, 0, 0, 0 };

const CellEngineAminoAcid Tyrosine = { "Y", "Tyr", "Tyrosine", "M_trnatyr_c", "M_tyrtrna_c", "M_tyr__L_c", { "TAC", "TAT" }, { "UAC", "UAU" }, { 509 }, { 613 }, 0, 0, 0 };

const CellEngineAminoAcid Valine = { "V", "Val", "Valine", "M_trnaval_c", "M_valtrna_c", "M_val__L_c", { "GTA", "GTC", "GTG", "GTT" }, { "GUA", "GUC", "GUG", "GUU" }, { 679 }, { 260 }, 0, 0, 0 };

inline std::unordered_map<std::string, const CellEngineAminoAcid*> CellEngineAminoAcids;

inline void MapAminoAcidsForProperCodons()
{
    CellEngineAminoAcids["GCA"] = &Alanine;
    CellEngineAminoAcids["GCC"] = &Alanine;
    CellEngineAminoAcids["GCG"] = &Alanine;
    CellEngineAminoAcids["GCT"] = &Alanine;
    CellEngineAminoAcids["GCU"] = &Alanine;

    CellEngineAminoAcids["AGA"] = &Arginine;
    CellEngineAminoAcids["AGG"] = &Arginine;
    CellEngineAminoAcids["CGA"] = &Arginine;
    CellEngineAminoAcids["CGC"] = &Arginine;
    CellEngineAminoAcids["CGG"] = &Arginine;
    CellEngineAminoAcids["CGT"] = &Arginine;
    CellEngineAminoAcids["CGU"] = &Arginine;

    CellEngineAminoAcids["AAC"] = &Asparagine;
    CellEngineAminoAcids["AAT"] = &Asparagine;
    CellEngineAminoAcids["AAU"] = &Asparagine;

    CellEngineAminoAcids["GAC"] = &AsparticAcid;
    CellEngineAminoAcids["GAT"] = &AsparticAcid;
    CellEngineAminoAcids["GAU"] = &AsparticAcid;

    CellEngineAminoAcids["TGC"] = &Cysteine;
    CellEngineAminoAcids["TGT"] = &Cysteine;
    CellEngineAminoAcids["UGC"] = &Cysteine;
    CellEngineAminoAcids["UGU"] = &Cysteine;

    CellEngineAminoAcids["CAA"] = &Glutamine;
    CellEngineAminoAcids["CAG"] = &Glutamine;

    CellEngineAminoAcids["GAA"] = &GlutamicAcid;
    CellEngineAminoAcids["GAG"] = &GlutamicAcid;

    CellEngineAminoAcids["GGA"] = &Glycine;
    CellEngineAminoAcids["GGC"] = &Glycine;
    CellEngineAminoAcids["GGG"] = &Glycine;
    CellEngineAminoAcids["GGT"] = &Glycine;
    CellEngineAminoAcids["GGU"] = &Glycine;

    CellEngineAminoAcids["CAC"] = &Histidine;
    CellEngineAminoAcids["CAT"] = &Histidine;
    CellEngineAminoAcids["CAU"] = &Histidine;

    CellEngineAminoAcids["ATA"] = &Isoleucine;
    CellEngineAminoAcids["ATC"] = &Isoleucine;
    CellEngineAminoAcids["ATT"] = &Isoleucine;
    CellEngineAminoAcids["AUA"] = &Isoleucine;
    CellEngineAminoAcids["AUC"] = &Isoleucine;
    CellEngineAminoAcids["AUU"] = &Isoleucine;

    CellEngineAminoAcids["CTA"] = &Leucine;
    CellEngineAminoAcids["CTC"] = &Leucine;
    CellEngineAminoAcids["CTG"] = &Leucine;
    CellEngineAminoAcids["CTT"] = &Leucine;
    CellEngineAminoAcids["TTA"] = &Leucine;
    CellEngineAminoAcids["TTG"] = &Leucine;
    CellEngineAminoAcids["CUA"] = &Leucine;
    CellEngineAminoAcids["CUC"] = &Leucine;
    CellEngineAminoAcids["CUG"] = &Leucine;
    CellEngineAminoAcids["CUU"] = &Leucine;
    CellEngineAminoAcids["UUA"] = &Leucine;
    CellEngineAminoAcids["UUG"] = &Leucine;

    CellEngineAminoAcids["AAA"] = &Lysine;
    CellEngineAminoAcids["AAG"] = &Lysine;

    CellEngineAminoAcids["ATG"] = &Methionine;
    CellEngineAminoAcids["AUG"] = &Methionine;

    CellEngineAminoAcids["TTC"] = &Phenylalanine;
    CellEngineAminoAcids["TTT"] = &Phenylalanine;
    CellEngineAminoAcids["UUC"] = &Phenylalanine;
    CellEngineAminoAcids["UUU"] = &Phenylalanine;

    CellEngineAminoAcids["CCA"] = &Proline;
    CellEngineAminoAcids["CCC"] = &Proline;
    CellEngineAminoAcids["CCG"] = &Proline;
    CellEngineAminoAcids["CCT"] = &Proline;
    CellEngineAminoAcids["CCU"] = &Proline;

    CellEngineAminoAcids["AGC"] = &Serine;
    CellEngineAminoAcids["AGT"] = &Serine;
    CellEngineAminoAcids["TCA"] = &Serine;
    CellEngineAminoAcids["TCC"] = &Serine;
    CellEngineAminoAcids["TCG"] = &Serine;
    CellEngineAminoAcids["TCT"] = &Serine;
    CellEngineAminoAcids["AGU"] = &Serine;
    CellEngineAminoAcids["UCA"] = &Serine;
    CellEngineAminoAcids["UCC"] = &Serine;
    CellEngineAminoAcids["UCG"] = &Serine;
    CellEngineAminoAcids["UCU"] = &Serine;

    CellEngineAminoAcids["ACA"] = &Threonine;
    CellEngineAminoAcids["ACC"] = &Threonine;
    CellEngineAminoAcids["ACG"] = &Threonine;
    CellEngineAminoAcids["ACT"] = &Threonine;
    CellEngineAminoAcids["ACU"] = &Threonine;

    CellEngineAminoAcids["TGG"] = &Tryptophan;
    CellEngineAminoAcids["UGG"] = &Tryptophan;

    CellEngineAminoAcids["TAC"] = &Tyrosine;
    CellEngineAminoAcids["TAT"] = &Tyrosine;
    CellEngineAminoAcids["UAC"] = &Tyrosine;
    CellEngineAminoAcids["UAU"] = &Tyrosine;

    CellEngineAminoAcids["GTA"] = &Valine;
    CellEngineAminoAcids["GTC"] = &Valine;
    CellEngineAminoAcids["GTG"] = &Valine;
    CellEngineAminoAcids["GTT"] = &Valine;
    CellEngineAminoAcids["GUA"] = &Valine;
    CellEngineAminoAcids["GUC"] = &Valine;
    CellEngineAminoAcids["GUG"] = &Valine;
    CellEngineAminoAcids["GUU"] = &Valine;
}

#endif
