#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <span>
#include <cmath>
#include <array>
#include <string>
#include <vector>

#include "Logger.h"

using SignedInt = int64_t;
using UnsignedInt = uint64_t;

class Combinations
{
private:
    const UnsignedInt MaxSizeOfArray = 8;

    const UnsignedInt MaxSizeOfCombinationsArray = 32;

    std::vector<UnsignedInt> CombinationsArray;
public:
    Combinations() : CombinationsArray(MaxSizeOfCombinationsArray)
    {
    }
private:
    static void ShowArray(std::span<UnsignedInt> Array, UnsignedInt StartPos, UnsignedInt Counter, std::string_view AdditionalText);
public:
    void AddTwoBitArrays(const std::vector<UnsignedInt>& MainArray, const std::vector<UnsignedInt>& ArrayToAdd, std::vector<UnsignedInt>& ResultArray) const;
    void Add1ToMainArray(const std::vector<UnsignedInt>& MainArray, const std::vector<UnsignedInt>& ArrayToAdd, std::vector<UnsignedInt>& ResultArray) const;
public:
    void Add1ToMainArrayByAddTwoBitsArraysTestDemo(UnsignedInt NumberOfBits) const;
    void Add1ToMainArrayTestDemo(UnsignedInt NumberOfBits) const;
    void AddTwoDifferentBitsArraysTestDemo() const;
private:
    void RotateWithPrint(UnsignedInt R, UnsignedInt N, UnsignedInt K);
public:
    void GenerateAllKElementsCombinationsFromNElementsFirstWayWithPrintTestDemo(UnsignedInt N, UnsignedInt K);
private:
    void Rotate(UnsignedInt R);
public:
    void GenerateAllKElementsCombinationsFromNElementsFirstWayTestDemo(UnsignedInt N, UnsignedInt K);
    void GenerateAllKElementsCombinationsFromNElementsFirstWayAndHalfTestDemo(UnsignedInt N, UnsignedInt K);
public:
    static UnsignedInt SetKBitsInNumber(UnsignedInt N, UnsignedInt K);
    static UnsignedInt SetAllBitsInNumber(UnsignedInt N, UnsignedInt K);
    static UnsignedInt NumberOfCombinations(UnsignedInt NP, UnsignedInt KP);
    static std::string CreateBoolStringFromInt64BitState(UnsignedInt Number);
    static UnsignedInt NextNumberWithTheSameNumberOf1Bits(UnsignedInt Number);
public:
    static void GenerateAllKElementsCombinationsFromNElementsSecondWayTestDemo(UnsignedInt N, UnsignedInt K);
};

class Permutations
{
public:
    /* Glowny algorytm polega na rotowaniu */
    /*
        1234 -> 1 234 ->
                rotacja 4 na pocz¹tek 34 czyli dostajê 43
                        1 2 34
                        1 2 43
               rotacja 4 na pocz¹tek potrójnego
               1 423 ->
               rotacja 3 na pocz¹tek 23 czyli dostajê 32
                        1 4 23
                        1 4 32
               rotacja 3 na pocz¹tek potrójnego
               1 342 ->
               rotacja 2 na pocz¹tek 42 czyli dostajê 24
                        1 3 42
                        1 3 24
        Potem rotacja 1234 na 4123 i od poczatku podzial na string krótszy gdzie znów ta rotacja
    */
public:
    static void RotateForPermutation(UnsignedInt Number, std::string& PermutationStringLocal);
public:
    void Permute(UnsignedInt Number, std::string& PermutationStringLocal, UnsignedInt& Line);
public:
    void PermutationsTestDemo(std::string& PermutationString);
};

class Variations
{
    /* To sa wariacje z powtorzeniami - kazdy k wyrazowy ciag ze zbioru n elementowego czyli n do k-tej */
    /* czyli tu 4 do 3 czyli k = 3 a n = 4 - mam 4 cyfry - to zbior n = 4 i mam wyrazy 3 elementowe czyli k = 3 */
private:
    static std::string GetResultAsString(UnsignedInt KWords, std::vector<UnsignedInt>& Repeat);
public:
    static void VariationsFrom1ToEndIn3LoopsTestDemo();
    static void VariationsFrom1ToEndIn1LoopTestDemo();
    static void VariationsFrom1ToMIn3LoopsTestDemo(const std::vector<UnsignedInt>& Lengths);
    static void VariationsFrom1ToMIn1LoopTestDemo(const std::vector<UnsignedInt>& Lengths);
    static void VariationsFromNtoMIn3LoopsTestDemo(const std::vector<UnsignedInt>& Starts, const std::vector<UnsignedInt>& Lengths);
    static void VariationsFromNtoMIn1LoopTestDemo(const std::vector<UnsignedInt>& Starts, const std::vector<UnsignedInt>& Lengths);
};

void ShowDemoTest();

#endif
