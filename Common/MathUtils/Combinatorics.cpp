
#include <span>
#include <cmath>
#include <array>
#include <string>
#include <vector>

#include "Logger.h"
#include "Combinatorics.h"

using namespace std;

using SignedInt = int64_t;
using UnsignedInt = uint64_t;

//const UnsignedInt MaxSizeOfArray = 8;
//
//const UnsignedInt MaxSizeOfCombinationArray = 16;
//
//UnsignedInt CombinationsArray[MaxSizeOfCombinationArray + 1];

class Combinations
{
//    static constexpr UnsignedInt MaxSizeOfArray = 8;
//
//    static constexpr UnsignedInt MaxSizeOfCombinationArray = 32;
//
//    static UnsignedInt CombinationsArray[MaxSizeOfCombinationArray + 1];
    const UnsignedInt MaxSizeOfArray = 8;

    const UnsignedInt MaxSizeOfCombinationArray = 32;

    //vector<UnsignedInt> CombinationsArray[MaxSizeOfCombinationArray + 1];
    vector<UnsignedInt> CombinationsArray;
public:
    Combinations() : CombinationsArray(MaxSizeOfCombinationArray)
    {
    }
public:

    static void ShowArray(span<UnsignedInt> Array, UnsignedInt StartPos, UnsignedInt Counter, string_view AdditionalText)
    {
        string ArrayStr;

        for (UnsignedInt Pos = StartPos; Pos < Counter; Pos++)
            ArrayStr += to_string(Array[Pos]);

        LoggersManagerObject.Log(STREAM(ArrayStr << AdditionalText));
    }

    //void AddTwoBitArrays(const UnsignedInt *MainArray, const UnsignedInt *ArrayToAdd, UnsignedInt *ResultArray) const
    void AddTwoBitArrays(const vector<UnsignedInt>& MainArray, const vector<UnsignedInt>& ArrayToAdd, vector<UnsignedInt>& ResultArray) const
    {
        UnsignedInt Remember = 0;

        for (SignedInt Pos = (SignedInt)MaxSizeOfArray - 1; Pos >= 0; Pos--)
        {
            if (MainArray[Pos] == 0 && ArrayToAdd[Pos] == 0)
            {
                ResultArray[Pos] = Remember;
                Remember = 0;
            }
            else
            if (MainArray[Pos] == 0 && ArrayToAdd[Pos] == 1 && Remember == 0)
            {
                ResultArray[Pos] = 1;
                Remember = 0;
            }
            else
            if (MainArray[Pos] == 0 && ArrayToAdd[Pos] == 1 && Remember == 1)
            {
                ResultArray[Pos] = 0;
                Remember = 1;
            }
            else
            if (MainArray[Pos] == 1 && ArrayToAdd[Pos] == 0 && Remember == 0)
            {
                ResultArray[Pos] = 1;
                Remember = 0;
            }
            else
            if (MainArray[Pos] == 1 && ArrayToAdd[Pos] == 0 && Remember == 1)
            {
                ResultArray[Pos] = 0;
                Remember = 1;
            }
            else
            if (MainArray[Pos] == 1 && ArrayToAdd[Pos] == 1 && Remember == 0)
            {
                ResultArray[Pos] = 0;
                Remember = 1;
            }
            else
            if (MainArray[Pos] == 1 && ArrayToAdd[Pos] == 1 && Remember == 1)
            {
                ResultArray[Pos] = 1;
                Remember = 1;
            }
        }
    }

    //void Add1ToMainArray(const UnsignedInt *MainArray, const UnsignedInt *ArrayToAdd, UnsignedInt *ResultArray) const
    void Add1ToMainArray(const vector<UnsignedInt>& MainArray, const vector<UnsignedInt>& ArrayToAdd, vector<UnsignedInt>& ResultArray) const
    {
        UnsignedInt Remember = 0;
        if (MainArray[MaxSizeOfArray - 1] == 0 && ArrayToAdd[MaxSizeOfArray - 1] == 0)
            ResultArray[MaxSizeOfArray - 1] = 0;
        else
        if (MainArray[MaxSizeOfArray - 1] == 0 && ArrayToAdd[MaxSizeOfArray - 1] == 1)
            ResultArray[MaxSizeOfArray - 1] = 1;
        else
        if (MainArray[MaxSizeOfArray - 1] == 1 && ArrayToAdd[MaxSizeOfArray - 1] == 0)
            ResultArray[MaxSizeOfArray - 1] = 1;
        else
        if (MainArray[MaxSizeOfArray - 1] == 1 && ArrayToAdd[MaxSizeOfArray - 1] == 1)
        {
            ResultArray[MaxSizeOfArray - 1] = 0;
            Remember = 1;
        }

        for (SignedInt Pos = (SignedInt)MaxSizeOfArray - 2; Pos >= 0; Pos--)
        {
            if (MainArray[Pos] == 0)
            {
                ResultArray[Pos] = Remember;
                Remember = 0;
            }
            else
            if (MainArray[Pos] == 1 && Remember == 0)
            {
                ResultArray[Pos] = 1;
                Remember = 0;
            }
            else
            if (MainArray[Pos] == 1 && Remember == 1)
            {
                ResultArray[Pos] = 0;
                Remember = 1;
            }
        }
    }

    void Add1ToMainArrayByAddTwoBitsArraysTestDemo(const UnsignedInt NumberOfBits) const
    {
        //UnsignedInt ArrayToAddCopy[MaxSizeOfArray];
        vector<UnsignedInt> ArrayToAddCopy(MaxSizeOfArray);

        //array<UnsignedInt> MainArray[MaxSizeOfArray] = {0, 0, 0, 0, 0, 0, 0, 0};
        vector<UnsignedInt> MainArray = { 0, 0, 0, 0, 0, 0, 0, 0 };
        vector<UnsignedInt> ArrayToAdd = { 0, 0, 0, 0, 0, 0, 0, 1 };

        //UnsignedInt ResultArray[MaxSizeOfArray];
        vector<UnsignedInt> ResultArray(MaxSizeOfArray);

        for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
            ArrayToAddCopy[Pos] = ArrayToAdd[Pos];

        for (UnsignedInt BitPos = 1; BitPos <= NumberOfBits; BitPos++)
        {
            for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
                ArrayToAdd[Pos] = ArrayToAddCopy[Pos];

            AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);

            ShowArray(ResultArray, 0, MaxSizeOfArray, "");

            for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
                MainArray[Pos] = ResultArray[Pos];
        }
    }

    void Add1ToMainArrayTestDemo(const UnsignedInt NumberOfBits) const
    {
        vector<UnsignedInt> ArrayToAddCopy(MaxSizeOfArray);

        vector<UnsignedInt> MainArray = { 0, 0, 0, 0, 0, 0, 0, 0 };
        vector<UnsignedInt> ArrayToAdd = { 0, 0, 0, 0, 0, 0, 0, 1 };

        vector<UnsignedInt> ResultArray(MaxSizeOfArray);

        for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
            ArrayToAddCopy[Pos] = ArrayToAdd[Pos];

        for (UnsignedInt BitPos = 1; BitPos <= NumberOfBits; BitPos++)
        {
            for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
                ArrayToAdd[Pos] = ArrayToAddCopy[Pos];

            Add1ToMainArray(MainArray, ArrayToAdd, ResultArray);

            ShowArray(ResultArray, 0, MaxSizeOfArray, "");

            for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
                MainArray[Pos] = ResultArray[Pos];
        }
    }

    void AddTwoDifferentBitsArraysTestDemo() const
    {
        vector<UnsignedInt> MainArray = { 0, 0, 0, 0, 1, 1, 0, 1 };
        vector<UnsignedInt> ArrayToAdd = { 0, 0, 0, 0, 0, 1, 1,  1 };

        vector<UnsignedInt> ResultArray(MaxSizeOfArray);

        AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);

        ShowArray(MainArray, 0, MaxSizeOfArray, "");
        ShowArray(ArrayToAdd, 0, MaxSizeOfArray, "");
        ShowArray(ResultArray, 0, MaxSizeOfArray, "");
    }

    /* COMBINATIONS */
    /*
    011101 6
    011101 PREV, j + 1 = 3
    011101 IN TRAIN PREV, R = 3 i = 3
    011101 IN TRAIN AFTER, R = 3 i = 3
    011101 IN TRAIN PREV, R = 3 i = 2
    001101 IN TRAIN AFTER, R = 3 i = 2
    101101 AFTER TOTAL, R = 3, WSTAWIENIE PAM NA MIEJSCE
    101101 AFTER
    */
    /* PONIZSZE DWIE FUNKCJE TO TO SAMO CO DWIE JESZCZE NASTEPNE TYLKO Z WYPISYWANIEM STANU ALGORYTMU przez ShowArray */

    void RotateWithPrint(UnsignedInt R, const UnsignedInt N, const UnsignedInt K)
    {
        UnsignedInt Pam = CombinationsArray[R];

        for (UnsignedInt i = R; i > 1; i--)
        {
            ShowArray(CombinationsArray, 1, N + 1, " IN TRAIN PREV, R = " + to_string(R) + " i = " + to_string(i));
            CombinationsArray[i] = CombinationsArray[i - 1];
            ShowArray(CombinationsArray, 1, N + 1, " IN TRAIN AFTER, R = " + to_string(R) + " i = " + to_string(i));
        }

        CombinationsArray[1] = Pam;

        ShowArray(CombinationsArray, 1, N + 1, " AFTER TOTAL, R = " + to_string(R));
    }

    void AllKElementsCombinationsFromNElementsFirstWayWithPrintTestDemo(const UnsignedInt N, const UnsignedInt K)
    {
        UnsignedInt Line = 1;

        for (UnsignedInt i = 1; i <= N; i++)
            CombinationsArray[i] = (i <= K);

        UnsignedInt j = N - 1;

        bool LoopCondition = true;

        while (LoopCondition)
        {
            ShowArray(CombinationsArray, 1, N + 1, " " + to_string(Line++));

            if (j >= N || N == K)
            {
                LoopCondition = false;
                goto EndLabel;
            }
            ShowArray(CombinationsArray, 1, N + 1, " PREV, j + 1 = " + to_string(j + 1));

            RotateWithPrint(j + 1, N, K);

            ShowArray(CombinationsArray, 1, N + 1, " AFTER");

            for (UnsignedInt i = 1; i < N; i++)
                if (CombinationsArray[i] == 0 && CombinationsArray[i + 1] == 1)
                {
                    j = i + 1;
                    break;
                }
            EndLabel:;
        }
    }

    void Rotate(UnsignedInt R)
    {
        UnsignedInt Pam = CombinationsArray[R];

        for (UnsignedInt i = R; i > 1; i--)
            CombinationsArray[i] = CombinationsArray[i - 1];

        CombinationsArray[1] = Pam;
    }

    void AllKElementsCombinationsFromNElementsFirstWayTestDemo(const UnsignedInt N, const UnsignedInt K)
    {
        UnsignedInt Line = 1;

        for (UnsignedInt i = 1; i <= N; i++)
            CombinationsArray[i] = i <= K;

        UnsignedInt j = N - 1;

        bool LoopCondition = true;

        while (LoopCondition == true)
        {
            ShowArray(CombinationsArray, 1, N + 1, " " + to_string(Line++));

            if (j >= N || N == K)
            {
                LoopCondition = false;
                goto EndLabel;
            }

            Rotate(j + 1);

            for (UnsignedInt i = 1; i < N; i++)
                if (CombinationsArray[i] == 0 && CombinationsArray[i + 1] == 1)
                {
                    j = i + 1;
                    break;
                }
            EndLabel:;
        }
    }

    void AllKElementsCombinationsFromNElementsFirstWayAndHalfTestDemo(const UnsignedInt N, const UnsignedInt K)
    {
        UnsignedInt Line = 1;

        for (UnsignedInt i = 1; i <= N; i++)
            CombinationsArray[i] = (i <= K);

        UnsignedInt j = N - 1;

        while (true)
        {
            ShowArray(CombinationsArray, 1, N + 1, " " + to_string(Line++));

            if (j >= N || N == K)
                break;

            Rotate(j + 1);

            for (UnsignedInt i = 1; i < N; i++)
                if (CombinationsArray[i] == 0 && CombinationsArray[i + 1] == 1)
                {
                    j = i + 1;
                    break;
                }
        }
    }

    static UnsignedInt NumberOfCombinations(UnsignedInt NP, UnsignedInt KP)
    {
        UnsignedInt NFactorial = 1;
        for (UnsignedInt Number = 1; Number <= NP; Number++)
            NFactorial *= Number;

        UnsignedInt KFactorial = 1;
        for (UnsignedInt Number = 1; Number <= KP; Number++)
            KFactorial *= Number;

        UnsignedInt NMinusKFactorial = 1;
        for (UnsignedInt Number = 1; Number <= (NP - KP); Number++)
            NMinusKFactorial *= Number;

        return NFactorial / (KFactorial * NMinusKFactorial);
    }

    static string CreateBoolStringFromInt64BitState(UnsignedInt Number)
    {
        string NumberString;

        UnsignedInt LenOfInt = sizeof(UnsignedInt) * 8;
        for (UnsignedInt i = LenOfInt - 1; i >= 1; i--)
            (Number >> (LenOfInt - i - 1)) << (LenOfInt - 1) ? NumberString += "1" : NumberString += "0";

        return NumberString;
    }

    #pragma warn -ngu
    static UnsignedInt NextNumberWithTheSameNumberOf1Bits(UnsignedInt Number)
    {
        UnsignedInt Smallest, Ripple, Ones;

        Smallest = Number & -Number;
        Ripple = Number + Smallest;
        Ones = Number ^ Ripple;
        Ones = (Ones >> 2) / Smallest;

        return Ripple | Ones;
    }

    void AllKElementsCombinationsFromNElementsSecondWayTestDemo(const UnsignedInt N, const UnsignedInt K)
    {
        UnsignedInt Line = 1;
        for (UnsignedInt i = 1; i <= N; i++)
            CombinationsArray[i] = (i <= K);

        UnsignedInt AimNumber = 0;
        for (UnsignedInt i = 1; i <= N; i++)
            AimNumber |= (UnsignedInt) pow(2, i - 1);

        ShowArray(CombinationsArray, 1, N + 1, "");

        LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(AimNumber)));
        LoggersManagerObject.Log(STREAM(""));

        UnsignedInt Number = 0;
        for (UnsignedInt i = 1; i <= N; i++)
            if (CombinationsArray[i] == 1)
                Number |= (UnsignedInt) pow(2, i - 1);

        LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) << "    " << to_string(Line++)));

        while (Number < AimNumber)
        {
            Number = NextNumberWithTheSameNumberOf1Bits(Number);
            if (Number < AimNumber)
                LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) << "    " << to_string(Line++)));
        }

        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM(to_string(NumberOfCombinations(N, K))));
    }
};







/* PERMUTATIONS */
/* Algorytm polega na rotowaniu */
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
    Potem rotacja 1234 na 4123 i od pocz¹tku podzia³ na string krótsdzy gdzie znów ta rotacja
*/

void RotateForPermutation(UnsignedInt Number, string& PermutationStringLocal)
{
    UnsignedInt PermutationStringLength = PermutationStringLocal.length();

    char TempChar = PermutationStringLocal[PermutationStringLength];

    for (UnsignedInt i = PermutationStringLength - 1; i >= Number; i--)
        PermutationStringLocal[i + 1] = PermutationStringLocal[i];

    PermutationStringLocal[Number] = TempChar;
}

void PermutationFunction(UnsignedInt Number, string& PermutationStringLocal, UnsignedInt& Line)
{
    if (Number == 1)
        LoggersManagerObject.Log(STREAM(PermutationStringLocal + " LINE = " + to_string(Line++)));
    else
    {
        UnsignedInt PermutationStringLength = PermutationStringLocal.length();

        for (UnsignedInt i = 1; i < Number; i++)
        {
            PermutationFunction(Number - 1, PermutationStringLocal, Line);
            RotateForPermutation(PermutationStringLength - Number + 1 + 1, PermutationStringLocal);
        }
    }
}

void Permutations1TestDemo()
{
    LoggersManagerObject.Log(STREAM("PERMUTATIONS:"));

    UnsignedInt Line = 1;
    string PermutationString = "ABCDE";
    UnsignedInt PermutationStringLength = PermutationString.length();
    PermutationFunction(PermutationStringLength + 1, PermutationString, Line);
}



















//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//VARIATIONS

void VariationsTestDemo()
{
    /* To s wariacje z powtórzeniami - ka¿dy k wyrazowy ci¹g ze zbioru n elementowego czyli n do k-tej */
    /* czyli tu 4 do 3 czyli k = 3 a n = 4 - mam 4 cyfry - to zbiór n = 4 i mam wyrazy 3 elementowe czyli k = 3 */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("WARIACJE Z POWTÓRZENIAMI W KILKU PÊTLACH:"));
    LoggersManagerObject.Log(STREAM(""));

    UnsignedInt Line = 1;
    for(UnsignedInt i1 = 1; i1 <= 4; i1++)
        for(UnsignedInt i2 = 1; i2 <= 4; i2++)
            for(UnsignedInt i3 = 1; i3 <= 4; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) << to_string(i2) << to_string(i3) << " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest n^k razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("WARIACJE Z POWTÓRZENIAMI W JEDNEJ PÊTLI ZAMIAST KILKU PÊTLI:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    UnsignedInt SizeOfSetN = 4;
    UnsignedInt KWords = 3;
    UnsignedInt TotalNumberOfVariations = pow(SizeOfSetN, KWords);
    vector<UnsignedInt> Repeat(KWords + 1);
    for(UnsignedInt ki = 1; ki <= KWords; ki++)
        Repeat[ki] = 1;
    for(UnsignedInt Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        /* Poni¿sze linijki to jedynie wyœwietlenie rezultatu */

        string SKWord = "";
        for(UnsignedInt ki = 1; ki <= KWords; ki++)
            SKWord += to_string(Repeat[ki]);
        LoggersManagerObject.Log(STREAM(SKWord << " LINE = " << to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzednia cyfra jest przewijana */

        Repeat[KWords]++;
        for(UnsignedInt ki = KWords; ki >= 1; ki--)
            if(Repeat[ki] % (SizeOfSetN + 1) == 0)
            {
                Repeat[ki] = 1;
                Repeat[ki - 1]++;
            }
    }
}

void LoopsFrom1ToMInOneLoopTestDemo()
{
    /* Zamiana kilku pêtli ró¿nej d³ugoœci zaczynanych od 1 o kroku 1 na jedn¹ pêtlê */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("W KILKU PÊTLACH RO¯NEJ D£UGOŒCI OD 1:"));
    LoggersManagerObject.Log(STREAM(""));

    UnsignedInt Line = 1;
    UnsignedInt Length[4] = {0, 4, 5, 7};

    for(UnsignedInt i1 = 1; i1 <= Length[1]; i1++)
        for(UnsignedInt i2 = 1; i2 <= Length[2]; i2++)
            for(UnsignedInt i3 = 1; i3 <= Length[3]; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) << to_string(i2) << to_string(i3) << " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest (Length[1] * Length[2] * Length[3]) razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("POPRZEDNIE PÊTLE W JEDNEJ PÊTLI ZAMIAST KILKU:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    UnsignedInt KWords = 3;
    UnsignedInt TotalNumberOfVariations = Length[1] * Length[2] * Length[3];
    vector<UnsignedInt> Repeat(KWords + 1);
    for(UnsignedInt ki = 1; ki <= KWords; ki++)
        Repeat[ki] = 1;
    for(UnsignedInt Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        /* Poni¿sze linijki to jedynie wyœwietlenie rezultatu */

        string SKWord = "";
        for(UnsignedInt ki = 1; ki <= KWords; ki++)
            SKWord += to_string(Repeat[ki]);
        LoggersManagerObject.Log(STREAM(SKWord + " LINE = " + to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzednia cyfra jest przewijana */

        Repeat[KWords]++;
        for(UnsignedInt ki = KWords; ki >= 1; --ki)
            if(Repeat[ki] % (Length[ki] + 1) == 0)
            {
                Repeat[ki] = 1;
                Repeat[ki - 1]++;
            }
    }
}

void LoopsFromNtoMInOneLoopTestDemo()
{
    /* Zamiana kilku pêtli ró¿nej d³ugoœci zaczynanych od Start[i] o kroku 1 na jedn¹ pêtlê */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("W KILKU PÊTLACH RO¯NEJ D£UGOŒCI OD 1:"));
    LoggersManagerObject.Log(STREAM(""));

    UnsignedInt Line = 1;
    UnsignedInt Length[4] = {0, 4, 5, 7};
    UnsignedInt Start[4] = {0, 2, 3, 5};

    for(UnsignedInt i1 = Start[1]; i1 <= Length[1]; i1++)
        for(UnsignedInt i2 = Start[2]; i2 <= Length[2]; i2++)
            for(UnsignedInt i3 = Start[3]; i3 <= Length[3]; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) << to_string(i2) << to_string(i3) << " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest (Length[1] * Length[2] * Length[3]) razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("POPRZEDNIE PÊTLE W JEDNEJ PÊTLI ZAMIAST KILKU:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    UnsignedInt KWords = 3;

    /* UnsignedInt TotalNumberOfVariations = (Length[1] - Start[1] + 1) * (Length[2] - Start[2] + 1) * (Length[3] - Start[3] + 1); */
    UnsignedInt TotalNumberOfVariations = 1;
    for(UnsignedInt ki = 1; ki <= KWords; ki++)
        TotalNumberOfVariations *= (Length[ki] - Start[ki] + 1);

    vector<UnsignedInt> Repeat(KWords + 1);
    for(UnsignedInt ki = 1; ki <= KWords; ki++)
        Repeat[ki] = Start[ki];
    for(UnsignedInt Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        /* Poni¿sze linijki to jedynie wyœwietlenie rezultatu */

        string SKWord;
        for(UnsignedInt ki = 1; ki <= KWords; ki++)
            SKWord += to_string(Repeat[ki]);
        LoggersManagerObject.Log(STREAM(SKWord << " LINE = " << to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzednia cyfra jest przewijana */

        Repeat[KWords]++;
        for(UnsignedInt ki = KWords; ki >= 1; --ki)
            if(Repeat[ki] % (Length[ki] + 1) == 0)
            {
                Repeat[ki] = Start[ki];
                Repeat[ki - 1]++;
            }
    }
}











void ShowDemoTest()
{
    Combinations CombinationsObject;

    LoggersManagerObject.Log(STREAM("T1"));
    CombinationsObject.Add1ToMainArrayByAddTwoBitsArraysTestDemo(16);
    LoggersManagerObject.Log(STREAM("T2"));
    CombinationsObject.Add1ToMainArrayTestDemo(16);
    LoggersManagerObject.Log(STREAM("T3"));
    CombinationsObject.AddTwoDifferentBitsArraysTestDemo();

    LoggersManagerObject.Log(STREAM("T4"));
    CombinationsObject.AllKElementsCombinationsFromNElementsFirstWayTestDemo(6, 4);
    LoggersManagerObject.Log(STREAM("T5"));
    CombinationsObject.AllKElementsCombinationsFromNElementsFirstWayAndHalfTestDemo(6, 4);
    LoggersManagerObject.Log(STREAM("T6"));
    CombinationsObject.AllKElementsCombinationsFromNElementsSecondWayTestDemo(6, 4);

    LoggersManagerObject.Log(STREAM("T7"));
    Permutations1TestDemo();

    LoggersManagerObject.Log(STREAM("T8"));
    VariationsTestDemo();
    LoggersManagerObject.Log(STREAM("T9"));
    LoopsFrom1ToMInOneLoopTestDemo();
    LoggersManagerObject.Log(STREAM("T10"));
    LoopsFromNtoMInOneLoopTestDemo();
}