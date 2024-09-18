
#ifndef CELL_ENGINE_AMINO_ACIDS_H
#define CELL_ENGINE_AMINO_ACIDS_H

#include <string>
#include <vector>
#include <unordered_map>

#include "CellEngineTypes.h"
#include "CellEngineParticlesKindsManager.h"

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

class CellEngineAminoAcidsManager
{
    const CellEngineAminoAcid Alanine = { "A", "Ala", 	"Alanine", 	 "M_trnaala_c","M_alatrna_c","M_ala__L_c",{ "GCA", "GCC", "GCG", "GCT" }, { "GCA", "GCC", "GCG", "GCU" }, { 719 }, { 163 }, 0, 0, 0 };
    const CellEngineAminoAcid Arginine = { "R", "Arg", "Arginine", "M_trnaarg_c", "M_argtrna_c",	"M_arg__L_c", { "AGA", "AGG", "CGA", "CGC", "CGG", "CGT" }, { "AGA", "AGG", "CGA", "CGC", "CGG", "CGU" }, {374, 717 }, { 535 }, 0, 0, 0 };
    const CellEngineAminoAcid Asparagine = { "N", "Asn", "Asparagine", "M_trnaasn_c", "M_asntrna_c", "M_asn__L_c",{ "AAC", "AAT" }, { "AAC", "AAU" }, { 681 }, { 76 }, 0, 0, 0 };
    const CellEngineAminoAcid AsparticAcid = { "D", "Asp", "AsparticAcid", "M_trnaasp_c" ,"M_asptrna_c", "M_asp__L_c", { "GAC", "GAT" }, { "GAC", "GAU" }, { 724 }, { 287 }, 0, 0, 0 };
    const CellEngineAminoAcid Cysteine = { "C", "Cys", 	"Cysteine",  "M_trnacys_c", "M_cystrna_c","M_cys__L_c" ,{ "TGC", "TGT" }, { "UGC", "UGU" }, { 828 }, { 837 }, 0, 0, 0 };
    const CellEngineAminoAcid Glutamine = { "Q", "Gln", "Glutamine", "M_trnagln_c", "M_glntrna_c", "M_gln__L_c", { "CAA", "CAG" }, { "CAA", "CAG" }, { 508 }, { 1 }, 0, 0, 0 };
    const CellEngineAminoAcid GlutamicAcid = { "E", "Glu", "GlutamicAcid", "M_trnaglu_c", "M_glutrna_c", "M_glu__L_c" ,{ "GAA", "GAG" }, { "GAA", "GAG" }, { 680 }, {126 }, 0, 0, 0 };
    const CellEngineAminoAcid Glycine = { "G", "Gly", "Glycine", "M_trnagly_c", "M_glytrna_c", "M_gly_c", { "GGA", "GGC", "GGG", "GGT" }, { "GGA", "GGC", "GGG", "GGU" }, { 295 }, { 405 }, 0, 0, 0 };
    const CellEngineAminoAcid Histidine = { "H", "His", "Histidine", "M_trnahis_c", "M_histrna_c", "M_his__L_c", { "CAC", "CAT" }, { "CAC", "CAU" }, { 624 }, { 288 }, 0, 0, 0 };
    const CellEngineAminoAcid Isoleucine = { "I", "Ile", "Isoleucine", "M_trnaile_c", "M_iletrna_c", "M_ile__L_c",	{ "ATA", "ATC", "ATT" }, { "AUA", "AUC", "AUU" }, { 635 }, { 519 }, 0, 0, 0 };
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

    const std::unordered_map<std::string, const CellEngineAminoAcid*> CellEngineNamesForAminoAcidsMap =
    {
        { "Alanine", &Alanine },
        { "Arginine", &Arginine },
        { "Asparagine", &Asparagine },
        { "AsparticAcid", &AsparticAcid },
        { "Cysteine", &Cysteine },
        { "Glutamine", &Glutamine },
        { "GlutamicAcid", &GlutamicAcid },
        { "Glycine", &Glycine },
        { "Histidine", &Histidine },
        { "Isoleucine", &Isoleucine },
        { "Leucine", &Leucine },
        { "Lysine", &Lysine },
        { "Methionine", &Methionine },
        { "Phenylalanine", &Phenylalanine },
        { "Proline", &Proline },
        { "Serine", &Serine },
        { "Threonine", &Threonine },
        { "Tryptophan", &Tryptophan },
        { "Tyrosine", &Tyrosine },
        { "Valine", &Valine}
    };

    std::unordered_map<std::string, const CellEngineAminoAcid*> CellEngineCodonsForAminoAcidsMap;

    std::unordered_map<EntityIdInt, const CellEngineAminoAcid*> CellEngineAminoAcidsIdMap;
    std::unordered_map<EntityIdInt, const CellEngineAminoAcid*> CellEnginetRNAChargedIdMap;
    std::unordered_map<EntityIdInt, const CellEngineAminoAcid*> CellEnginetRNAUnchargedIdMap;

    static void SetNumIdFromName(const CellEngineAminoAcid& CellEngineAminoAcidObject)
    {
        const_cast<CellEngineAminoAcid&>(CellEngineAminoAcidObject).EntityId = ParticlesKindsManagerObject.GetParticleKindFromStrId(CellEngineAminoAcidObject.Name)->EntityId;
        const_cast<CellEngineAminoAcid&>(CellEngineAminoAcidObject).tRNAUnchargedNumId = ParticlesKindsManagerObject.GetParticleKindFromStrId(CellEngineAminoAcidObject.tRNAUnchargedName)->EntityId;
        const_cast<CellEngineAminoAcid&>(CellEngineAminoAcidObject).tRNAChargedNumId = ParticlesKindsManagerObject.GetParticleKindFromStrId(CellEngineAminoAcidObject.tRNAChargedName)->EntityId;
    }
public:
    void MapAminoAcidsIdToAminoAcidsObject()
    {
        for (auto& AminoAcidObject : CellEngineNamesForAminoAcidsMap)
        {
            SetNumIdFromName(*AminoAcidObject.second);

            CellEngineAminoAcidsIdMap[AminoAcidObject.second->EntityId] = AminoAcidObject.second;
            CellEnginetRNAChargedIdMap[AminoAcidObject.second->tRNAChargedNumId] = AminoAcidObject.second;
            CellEnginetRNAChargedIdMap[AminoAcidObject.second->tRNAUnchargedNumId] = AminoAcidObject.second;
        }
    }

    void MapAminoAcidsForProperCodons()
    {
        CellEngineCodonsForAminoAcidsMap["GCA"] = &Alanine;
        CellEngineCodonsForAminoAcidsMap["GCC"] = &Alanine;
        CellEngineCodonsForAminoAcidsMap["GCG"] = &Alanine;
        CellEngineCodonsForAminoAcidsMap["GCT"] = &Alanine;
        CellEngineCodonsForAminoAcidsMap["GCU"] = &Alanine;

        CellEngineCodonsForAminoAcidsMap["AGA"] = &Arginine;
        CellEngineCodonsForAminoAcidsMap["AGG"] = &Arginine;
        CellEngineCodonsForAminoAcidsMap["CGA"] = &Arginine;
        CellEngineCodonsForAminoAcidsMap["CGC"] = &Arginine;
        CellEngineCodonsForAminoAcidsMap["CGG"] = &Arginine;
        CellEngineCodonsForAminoAcidsMap["CGT"] = &Arginine;
        CellEngineCodonsForAminoAcidsMap["CGU"] = &Arginine;

        CellEngineCodonsForAminoAcidsMap["AAC"] = &Asparagine;
        CellEngineCodonsForAminoAcidsMap["AAT"] = &Asparagine;
        CellEngineCodonsForAminoAcidsMap["AAU"] = &Asparagine;

        CellEngineCodonsForAminoAcidsMap["GAC"] = &AsparticAcid;
        CellEngineCodonsForAminoAcidsMap["GAT"] = &AsparticAcid;
        CellEngineCodonsForAminoAcidsMap["GAU"] = &AsparticAcid;

        CellEngineCodonsForAminoAcidsMap["TGC"] = &Cysteine;
        CellEngineCodonsForAminoAcidsMap["TGT"] = &Cysteine;
        CellEngineCodonsForAminoAcidsMap["UGC"] = &Cysteine;
        CellEngineCodonsForAminoAcidsMap["UGU"] = &Cysteine;

        CellEngineCodonsForAminoAcidsMap["CAA"] = &Glutamine;
        CellEngineCodonsForAminoAcidsMap["CAG"] = &Glutamine;

        CellEngineCodonsForAminoAcidsMap["GAA"] = &GlutamicAcid;
        CellEngineCodonsForAminoAcidsMap["GAG"] = &GlutamicAcid;

        CellEngineCodonsForAminoAcidsMap["GGA"] = &Glycine;
        CellEngineCodonsForAminoAcidsMap["GGC"] = &Glycine;
        CellEngineCodonsForAminoAcidsMap["GGG"] = &Glycine;
        CellEngineCodonsForAminoAcidsMap["GGT"] = &Glycine;
        CellEngineCodonsForAminoAcidsMap["GGU"] = &Glycine;

        CellEngineCodonsForAminoAcidsMap["CAC"] = &Histidine;
        CellEngineCodonsForAminoAcidsMap["CAT"] = &Histidine;
        CellEngineCodonsForAminoAcidsMap["CAU"] = &Histidine;

        CellEngineCodonsForAminoAcidsMap["ATA"] = &Isoleucine;
        CellEngineCodonsForAminoAcidsMap["ATC"] = &Isoleucine;
        CellEngineCodonsForAminoAcidsMap["ATT"] = &Isoleucine;
        CellEngineCodonsForAminoAcidsMap["AUA"] = &Isoleucine;
        CellEngineCodonsForAminoAcidsMap["AUC"] = &Isoleucine;
        CellEngineCodonsForAminoAcidsMap["AUU"] = &Isoleucine;

        CellEngineCodonsForAminoAcidsMap["CTA"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CTC"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CTG"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CTT"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["TTA"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["TTG"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CUA"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CUC"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CUG"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["CUU"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["UUA"] = &Leucine;
        CellEngineCodonsForAminoAcidsMap["UUG"] = &Leucine;

        CellEngineCodonsForAminoAcidsMap["AAA"] = &Lysine;
        CellEngineCodonsForAminoAcidsMap["AAG"] = &Lysine;

        CellEngineCodonsForAminoAcidsMap["ATG"] = &Methionine;
        CellEngineCodonsForAminoAcidsMap["AUG"] = &Methionine;

        CellEngineCodonsForAminoAcidsMap["TTC"] = &Phenylalanine;
        CellEngineCodonsForAminoAcidsMap["TTT"] = &Phenylalanine;
        CellEngineCodonsForAminoAcidsMap["UUC"] = &Phenylalanine;
        CellEngineCodonsForAminoAcidsMap["UUU"] = &Phenylalanine;

        CellEngineCodonsForAminoAcidsMap["CCA"] = &Proline;
        CellEngineCodonsForAminoAcidsMap["CCC"] = &Proline;
        CellEngineCodonsForAminoAcidsMap["CCG"] = &Proline;
        CellEngineCodonsForAminoAcidsMap["CCT"] = &Proline;
        CellEngineCodonsForAminoAcidsMap["CCU"] = &Proline;

        CellEngineCodonsForAminoAcidsMap["AGC"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["AGT"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["TCA"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["TCC"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["TCG"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["TCT"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["AGU"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["UCA"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["UCC"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["UCG"] = &Serine;
        CellEngineCodonsForAminoAcidsMap["UCU"] = &Serine;

        CellEngineCodonsForAminoAcidsMap["ACA"] = &Threonine;
        CellEngineCodonsForAminoAcidsMap["ACC"] = &Threonine;
        CellEngineCodonsForAminoAcidsMap["ACG"] = &Threonine;
        CellEngineCodonsForAminoAcidsMap["ACT"] = &Threonine;
        CellEngineCodonsForAminoAcidsMap["ACU"] = &Threonine;

        CellEngineCodonsForAminoAcidsMap["TGG"] = &Tryptophan;
        CellEngineCodonsForAminoAcidsMap["UGG"] = &Tryptophan;

        CellEngineCodonsForAminoAcidsMap["TAC"] = &Tyrosine;
        CellEngineCodonsForAminoAcidsMap["TAT"] = &Tyrosine;
        CellEngineCodonsForAminoAcidsMap["UAC"] = &Tyrosine;
        CellEngineCodonsForAminoAcidsMap["UAU"] = &Tyrosine;

        CellEngineCodonsForAminoAcidsMap["GTA"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GTC"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GTG"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GTT"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GUA"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GUC"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GUG"] = &Valine;
        CellEngineCodonsForAminoAcidsMap["GUU"] = &Valine;
    }

    bool IsAminoAcid(const EntityIdInt EntityId)
    {
        return CellEngineAminoAcidsIdMap.find(EntityId) != CellEngineAminoAcidsIdMap.end();
    }

    bool IstRNAUncharged(const EntityIdInt EntityId)
    {
        return CellEnginetRNAUnchargedIdMap.find(EntityId) != CellEnginetRNAUnchargedIdMap.end();
    }

    bool IstRNACharged(const EntityIdInt EntityId)
    {
        return CellEnginetRNAChargedIdMap.find(EntityId) != CellEnginetRNAChargedIdMap.end();
    }

    bool IstRNAWithAminoAcidForCodon(const EntityIdInt ParticleEntityId, const std::string& Codon)
    {
        const auto AminoAcidIter = CellEngineCodonsForAminoAcidsMap.find(Codon);
        if (AminoAcidIter != CellEngineCodonsForAminoAcidsMap.end())
            return (AminoAcidIter->second->EntityId == ParticleEntityId);

        return false;
    }
};

inline CellEngineAminoAcidsManager CellEngineAminoAcidsManagerObject;

#endif
