
#ifndef CELL_ENGINE_AMINO_ACIDS_H
#define CELL_ENGINE_AMINO_ACIDS_H

#include <string>
#include <vector>
#include <unordered_map>

#include "CellEngineTypes.h"
#include "CellEngineParticlesKindsManager.h"

namespace CellEngineAminoAcids
{
    class CellEngineAminoAcid
    {
    public:
        std::string OneLetterCode{};
        std::string Abbreviation{};
        std::string FullName{};
        std::string tRNAUnchargedName{};
        std::string tRNAChargedName{};
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

    inline std::unordered_map<std::string, const CellEngineAminoAcid*> CellEngineAminoAcidsMap;

    inline void MapAminoAcidsForProperCodons()
    {
        CellEngineAminoAcidsMap["GCA"] = &Alanine;
        CellEngineAminoAcidsMap["GCC"] = &Alanine;
        CellEngineAminoAcidsMap["GCG"] = &Alanine;
        CellEngineAminoAcidsMap["GCT"] = &Alanine;
        CellEngineAminoAcidsMap["GCU"] = &Alanine;

        CellEngineAminoAcidsMap["AGA"] = &Arginine;
        CellEngineAminoAcidsMap["AGG"] = &Arginine;
        CellEngineAminoAcidsMap["CGA"] = &Arginine;
        CellEngineAminoAcidsMap["CGC"] = &Arginine;
        CellEngineAminoAcidsMap["CGG"] = &Arginine;
        CellEngineAminoAcidsMap["CGT"] = &Arginine;
        CellEngineAminoAcidsMap["CGU"] = &Arginine;

        CellEngineAminoAcidsMap["AAC"] = &Asparagine;
        CellEngineAminoAcidsMap["AAT"] = &Asparagine;
        CellEngineAminoAcidsMap["AAU"] = &Asparagine;

        CellEngineAminoAcidsMap["GAC"] = &AsparticAcid;
        CellEngineAminoAcidsMap["GAT"] = &AsparticAcid;
        CellEngineAminoAcidsMap["GAU"] = &AsparticAcid;

        CellEngineAminoAcidsMap["TGC"] = &Cysteine;
        CellEngineAminoAcidsMap["TGT"] = &Cysteine;
        CellEngineAminoAcidsMap["UGC"] = &Cysteine;
        CellEngineAminoAcidsMap["UGU"] = &Cysteine;

        CellEngineAminoAcidsMap["CAA"] = &Glutamine;
        CellEngineAminoAcidsMap["CAG"] = &Glutamine;

        CellEngineAminoAcidsMap["GAA"] = &GlutamicAcid;
        CellEngineAminoAcidsMap["GAG"] = &GlutamicAcid;

        CellEngineAminoAcidsMap["GGA"] = &Glycine;
        CellEngineAminoAcidsMap["GGC"] = &Glycine;
        CellEngineAminoAcidsMap["GGG"] = &Glycine;
        CellEngineAminoAcidsMap["GGT"] = &Glycine;
        CellEngineAminoAcidsMap["GGU"] = &Glycine;

        CellEngineAminoAcidsMap["CAC"] = &Histidine;
        CellEngineAminoAcidsMap["CAT"] = &Histidine;
        CellEngineAminoAcidsMap["CAU"] = &Histidine;

        CellEngineAminoAcidsMap["ATA"] = &Isoleucine;
        CellEngineAminoAcidsMap["ATC"] = &Isoleucine;
        CellEngineAminoAcidsMap["ATT"] = &Isoleucine;
        CellEngineAminoAcidsMap["AUA"] = &Isoleucine;
        CellEngineAminoAcidsMap["AUC"] = &Isoleucine;
        CellEngineAminoAcidsMap["AUU"] = &Isoleucine;

        CellEngineAminoAcidsMap["CTA"] = &Leucine;
        CellEngineAminoAcidsMap["CTC"] = &Leucine;
        CellEngineAminoAcidsMap["CTG"] = &Leucine;
        CellEngineAminoAcidsMap["CTT"] = &Leucine;
        CellEngineAminoAcidsMap["TTA"] = &Leucine;
        CellEngineAminoAcidsMap["TTG"] = &Leucine;
        CellEngineAminoAcidsMap["CUA"] = &Leucine;
        CellEngineAminoAcidsMap["CUC"] = &Leucine;
        CellEngineAminoAcidsMap["CUG"] = &Leucine;
        CellEngineAminoAcidsMap["CUU"] = &Leucine;
        CellEngineAminoAcidsMap["UUA"] = &Leucine;
        CellEngineAminoAcidsMap["UUG"] = &Leucine;

        CellEngineAminoAcidsMap["AAA"] = &Lysine;
        CellEngineAminoAcidsMap["AAG"] = &Lysine;

        CellEngineAminoAcidsMap["ATG"] = &Methionine;
        CellEngineAminoAcidsMap["AUG"] = &Methionine;

        CellEngineAminoAcidsMap["TTC"] = &Phenylalanine;
        CellEngineAminoAcidsMap["TTT"] = &Phenylalanine;
        CellEngineAminoAcidsMap["UUC"] = &Phenylalanine;
        CellEngineAminoAcidsMap["UUU"] = &Phenylalanine;

        CellEngineAminoAcidsMap["CCA"] = &Proline;
        CellEngineAminoAcidsMap["CCC"] = &Proline;
        CellEngineAminoAcidsMap["CCG"] = &Proline;
        CellEngineAminoAcidsMap["CCT"] = &Proline;
        CellEngineAminoAcidsMap["CCU"] = &Proline;

        CellEngineAminoAcidsMap["AGC"] = &Serine;
        CellEngineAminoAcidsMap["AGT"] = &Serine;
        CellEngineAminoAcidsMap["TCA"] = &Serine;
        CellEngineAminoAcidsMap["TCC"] = &Serine;
        CellEngineAminoAcidsMap["TCG"] = &Serine;
        CellEngineAminoAcidsMap["TCT"] = &Serine;
        CellEngineAminoAcidsMap["AGU"] = &Serine;
        CellEngineAminoAcidsMap["UCA"] = &Serine;
        CellEngineAminoAcidsMap["UCC"] = &Serine;
        CellEngineAminoAcidsMap["UCG"] = &Serine;
        CellEngineAminoAcidsMap["UCU"] = &Serine;

        CellEngineAminoAcidsMap["ACA"] = &Threonine;
        CellEngineAminoAcidsMap["ACC"] = &Threonine;
        CellEngineAminoAcidsMap["ACG"] = &Threonine;
        CellEngineAminoAcidsMap["ACT"] = &Threonine;
        CellEngineAminoAcidsMap["ACU"] = &Threonine;

        CellEngineAminoAcidsMap["TGG"] = &Tryptophan;
        CellEngineAminoAcidsMap["UGG"] = &Tryptophan;

        CellEngineAminoAcidsMap["TAC"] = &Tyrosine;
        CellEngineAminoAcidsMap["TAT"] = &Tyrosine;
        CellEngineAminoAcidsMap["UAC"] = &Tyrosine;
        CellEngineAminoAcidsMap["UAU"] = &Tyrosine;

        CellEngineAminoAcidsMap["GTA"] = &Valine;
        CellEngineAminoAcidsMap["GTC"] = &Valine;
        CellEngineAminoAcidsMap["GTG"] = &Valine;
        CellEngineAminoAcidsMap["GTT"] = &Valine;
        CellEngineAminoAcidsMap["GUA"] = &Valine;
        CellEngineAminoAcidsMap["GUC"] = &Valine;
        CellEngineAminoAcidsMap["GUG"] = &Valine;
        CellEngineAminoAcidsMap["GUU"] = &Valine;
    }

    inline void SetNumIdFromName(const CellEngineAminoAcid& CellEngineAminoAcidObject)
    {
        const_cast<CellEngineAminoAcid&>(CellEngineAminoAcidObject).EntityId = ParticlesKindsManagerObject.GetParticleKindFromStrId(CellEngineAminoAcidObject.Name)->EntityId;
        const_cast<CellEngineAminoAcid&>(CellEngineAminoAcidObject).tRNAUnchargedNumId = ParticlesKindsManagerObject.GetParticleKindFromStrId(CellEngineAminoAcidObject.tRNAUnchargedName)->EntityId;
        const_cast<CellEngineAminoAcid&>(CellEngineAminoAcidObject).tRNAChargedNumId = ParticlesKindsManagerObject.GetParticleKindFromStrId(CellEngineAminoAcidObject.tRNAChargedName)->EntityId;
    }

    inline void SetProperParticleKindEntityIdForAminoAcids()
    {
        SetNumIdFromName(Alanine);
        SetNumIdFromName(Arginine);
        SetNumIdFromName(Asparagine);
        SetNumIdFromName(AsparticAcid);
        SetNumIdFromName(Cysteine);
        SetNumIdFromName(Glutamine);
        SetNumIdFromName(GlutamicAcid);
        SetNumIdFromName(Glycine);
        SetNumIdFromName(Histidine);
        SetNumIdFromName(Isoleucine);
        SetNumIdFromName(Leucine);
        SetNumIdFromName(Lysine);
        SetNumIdFromName(Methionine);
        SetNumIdFromName(Phenylalanine);
        SetNumIdFromName(Proline);
        SetNumIdFromName(Serine);
        SetNumIdFromName(Threonine);
        SetNumIdFromName(Tryptophan);
        SetNumIdFromName(Tyrosine);
        SetNumIdFromName(Valine);
    }

    inline bool IsAminoAcid(const EntityIdInt EntityId)
    {

        return true;
    }

    inline bool IstRNAWithAminoAcidForCodon(const EntityIdInt ParticleEntityId, const std::string& Codon)
    {
        const auto AminoAcidIter = CellEngineAminoAcidsMap.find(Codon);
        if (AminoAcidIter != CellEngineAminoAcidsMap.end())
            return (AminoAcidIter->second->EntityId == ParticleEntityId);

        return false;
    }
};

#endif
