
#ifndef CELL_ENGINE_NUCLEIC_ACIDS_BASIC_OPERATIONS_H
#define CELL_ENGINE_NUCLEIC_ACIDS_BASIC_OPERATIONS_H

#include "CellEngineParticle.h"

using namespace std;

class CellEngineNucleicAcidsBasicOperations
{
protected:
    static inline string GetPairedSequenceStr(string SequenceStr)
    {
        try
        {
            for (auto& SequenceNucleotideChar : SequenceStr)
                SequenceNucleotideChar = CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(SequenceNucleotideChar)));
        }
        CATCH("getting paired sequence")

        return SequenceStr;
    }
protected:
    static inline Particle* GoSomeNucleotides(Particle* Particle::*Direction, const UnsignedInt LengthOfSequence, Particle& ParticleObjectForReaction)
    {
        Particle* ParticlePtr = &ParticleObjectForReaction;

        try
        {
            UnsignedInt NucleotidesCounter = 1;

            while (NucleotidesCounter < LengthOfSequence + 1 && ParticleObjectForReaction.*Direction != nullptr && ParticlePtr != nullptr)
            {
                ParticlePtr = ParticlePtr->*Direction;
                NucleotidesCounter++;
            }
        }
        CATCH("going some nucleotides")

        return ParticlePtr;
    }
protected:
    static inline Particle* GoToGenomeIndex(Particle* Particle::*Direction, const UnsignedInt TargetGenomeIndex, const UnsignedInt MaximalLengthOfSequence, Particle& ParticleObjectForReaction)
    {
        Particle* ParticlePtr = &ParticleObjectForReaction;

        try
        {
            UnsignedInt NucleotidesCounter = 1;

            while (NucleotidesCounter < MaximalLengthOfSequence + 1 && ParticleObjectForReaction.*Direction != nullptr && ParticlePtr != nullptr && ParticlePtr->GenomeIndex != TargetGenomeIndex)
            {
                ParticlePtr = ParticlePtr->*Direction;
                NucleotidesCounter++;
            }
        }
        CATCH("going some nucleotides")

        return ParticlePtr;
    }
protected:
    static inline tuple<Particle*, Particle*, UnsignedInt, string, vector<ChainIdInt>> GetNucleotidesSequence(Particle* Particle::*Direction, const UnsignedInt LengthOfSequence, Particle& ParticleObjectForReaction, const bool ToString, const bool ToVector, bool (*Predicate)(const Particle*))
    {
        string NucleotidesSequenceToCompareString;
        vector<ChainIdInt> NucleotidesSequenceToCompareVector;

        Particle* ParticlePtr = &ParticleObjectForReaction;
        Particle* ParticlePtrPrev = &ParticleObjectForReaction;
        UnsignedInt NucleotidesCounter = 1;

        try
        {
            while (NucleotidesCounter < LengthOfSequence + 1 && ParticleObjectForReaction.*Direction != nullptr && ParticlePtr != nullptr && Predicate(ParticlePtr) == true)
            {
                if (ToString == true)
                    NucleotidesSequenceToCompareString += CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticlePtr->ChainId);
                if (ToVector == true)
                    NucleotidesSequenceToCompareVector.emplace_back(ParticlePtr->ChainId);
                ParticlePtrPrev = ParticlePtr;
                ParticlePtr = ParticlePtr->*Direction;
                NucleotidesCounter++;
            }
        }
        CATCH("getting nucleotides sequence forward")

        return { ParticlePtr, ParticlePtrPrev, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector };
    }
protected:
    static inline void CutDNAPrev(Particle* NucleotideObjectForReactionPtr)
    {
        try
        {
            if (NucleotideObjectForReactionPtr->Prev != nullptr && NucleotideObjectForReactionPtr->Prev->Next != nullptr)
            {
                NucleotideObjectForReactionPtr->Prev->Next = nullptr;
                NucleotideObjectForReactionPtr->Prev = nullptr;
            }
        }
        CATCH("cutting DNA")
    }
protected:
    static inline void CutDNANext(Particle* NucleotideObjectForReactionPtr)
    {
        try
        {
            if (NucleotideObjectForReactionPtr->Next != nullptr && NucleotideObjectForReactionPtr->Next->Prev != nullptr)
            {
                NucleotideObjectForReactionPtr->Next->Prev = nullptr;
                NucleotideObjectForReactionPtr->Next = nullptr;
            }
        }
        CATCH("cutting DNA")
    }
protected:
    static inline void LinkDNA(Particle* NucleotideObjectForReactionPtr1, Particle* NucleotideObjectForReactionPtr2)
    {
        try
        {
            if (NucleotideObjectForReactionPtr1 != nullptr && NucleotideObjectForReactionPtr2 != nullptr)
            {
                NucleotideObjectForReactionPtr1->Prev = NucleotideObjectForReactionPtr2;
                NucleotideObjectForReactionPtr2->Next = NucleotideObjectForReactionPtr1;
            }
        }
        CATCH("linking DNA")
    }
protected:
    static inline void SeparateTwoPairedDNANucleotides(Particle* ParticleNucleotide)
    {
        try
        {
            if (ParticleNucleotide->PairedNucleotidePtr != nullptr && ParticleNucleotide->PairedNucleotidePtr->PairedNucleotidePtr != nullptr)
            {
                ParticleNucleotide->PairedNucleotidePtr->PairedNucleotidePtr = nullptr;
                ParticleNucleotide->PairedNucleotidePtr = nullptr;
            }
        }
        CATCH("separating dna strands")
    }
protected:
    static inline void SeparateDNAStrands(Particle* Particle::*Direction, Particle* NucleotidePtr, const UnsignedInt LengthOfStrand)
    {
        try
        {
            for (UnsignedInt DiffPos = 0; DiffPos < LengthOfStrand; DiffPos++)
            {
                SeparateTwoPairedDNANucleotides(NucleotidePtr);
                NucleotidePtr = NucleotidePtr->*Direction;
            }
        }
        CATCH("separating dna strands")
    }
protected:
    static inline void JoinDNAStrands(Particle* Particle::*Direction, Particle* Strand1, Particle* Strand2)
    {
        try
        {
            while (Strand1->PairedNucleotidePtr == nullptr)
            {
                Strand1->PairedNucleotidePtr = Strand2;
                Strand2->PairedNucleotidePtr = Strand1;
                Strand1 = Strand1->*Direction;
                Strand2 = Strand2->*Direction;
            }
        }
        CATCH("joining dna strands")
    }
protected:
    static inline void DeleteLinkedParticlesPointersList(Particle& ParticleObject)
    {
        try
        {
            for (auto &LinkedParticlesPointersListObject: ParticleObject.LinkedParticlesPointersList)
                if (LinkedParticlesPointersListObject != nullptr)
                {
                    for (auto &LinkedParticlesPointersListObjectInternal: LinkedParticlesPointersListObject->LinkedParticlesPointersList)
                        if (LinkedParticlesPointersListObjectInternal->Index == ParticleObject.Index)
                            LinkedParticlesPointersListObjectInternal = nullptr;
                    erase_if(LinkedParticlesPointersListObject->LinkedParticlesPointersList, [](const Particle *item){ return (item == nullptr); });

                    LinkedParticlesPointersListObject = nullptr;
                }
            ParticleObject.LinkedParticlesPointersList.clear();
        }
        CATCH("delete linked particles pointers list");
    }
};

#endif
