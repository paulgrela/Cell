
#include <cmath>
#include <string>
#include <vector>

#include "Logger.h"
#include "Combinatorics.h"

using SignedInt = int64_t;
using UnsignedInt = uint64_t;

using namespace std;

const UnsignedInt MaxArr = 8;

void ShowArray(UnsignedInt* Array, UnsignedInt Start, UnsignedInt Counter, string_view LineStr)
{
    string s;
    for (UnsignedInt p = Start; p < Counter; p++)
        s += to_string(Array[p]);
    LoggersManagerObject.Log(STREAM(s << LineStr));
}

void AddTwoBitArrays(const UnsignedInt* MainArray, const UnsignedInt* ArrayToAdd, UnsignedInt* ResultArray)
{
    UnsignedInt Remember = 0;
    for(SignedInt Pos = MaxArr - 1; Pos >= 0; Pos--)
    {
        if(MainArray[Pos] == 0 && ArrayToAdd[Pos] == 0)
        {
            ResultArray[Pos] = Remember;
            Remember = 0;
        }
        else
        if(MainArray[Pos] == 0 && ArrayToAdd[Pos] == 1 && Remember == 0)
        {
            ResultArray[Pos] = 1;
            Remember = 0;
        }
        else
        if(MainArray[Pos] == 0 && ArrayToAdd[Pos] == 1 && Remember == 1)
        {
            ResultArray[Pos] = 0;
            Remember = 1;
        }
        else
        if(MainArray[Pos] == 1 && ArrayToAdd[Pos] == 0 && Remember == 0)
        {
            ResultArray[Pos] = 1;
            Remember = 0;
        }
        else
        if(MainArray[Pos] == 1 && ArrayToAdd[Pos] == 0 && Remember == 1)
        {
            ResultArray[Pos] = 0;
            Remember = 1;
        }
        else
        if(MainArray[Pos] == 1 && ArrayToAdd[Pos] == 1 && Remember == 0)
        {
            ResultArray[Pos] = 0;
            Remember = 1;
        }
        else
        if(MainArray[Pos] == 1 && ArrayToAdd[Pos] == 1 && Remember == 1)
        {
            ResultArray[Pos] = 1;
            Remember = 1;
        }
    }
}

void Add1ToMainArray(const UnsignedInt* MainArray, const UnsignedInt* ArrayToAdd, UnsignedInt* ResultArray)
{
    UnsignedInt Remember = 0;
    if(MainArray[MaxArr - 1] == 0 && ArrayToAdd[MaxArr - 1] == 0)
        ResultArray[MaxArr - 1] = 0;
    else
    if(MainArray[MaxArr - 1] == 0 && ArrayToAdd[MaxArr - 1] == 1)
        ResultArray[MaxArr - 1] = 1;
    else
    if(MainArray[MaxArr - 1] == 1 && ArrayToAdd[MaxArr - 1] == 0)
        ResultArray[MaxArr - 1] = 1;
    else
    if(MainArray[MaxArr - 1] == 1 && ArrayToAdd[MaxArr - 1] == 1)
    {
        ResultArray[MaxArr - 1] = 0;
        Remember = 1;
    }

    for (SignedInt Pos = MaxArr - 2; Pos >= 0; Pos--)
    {
        if(MainArray[Pos] == 0)
        {
            ResultArray[Pos] = Remember;
            Remember = 0;
        }
        else
        if(MainArray[Pos] == 1 && Remember == 0)
        {
            ResultArray[Pos] = 1;
            Remember = 0;
        }
        else
        if(MainArray[Pos] == 1 && Remember == 1)
        {
            ResultArray[Pos] = 0;
            Remember = 1;
        }
    }
}

void Add1ToMainArrayByAddTwoBitsArraysButtonClick()
{
    UnsignedInt MainArrayCopy[MaxArr];
    UnsignedInt ArrayToAddCopy[MaxArr];

    UnsignedInt MainArray[MaxArr] =  { 0,0,0,0, 0, 0, 0, 0 };
    UnsignedInt ArrayToAdd[MaxArr] = { 0,0,0,0, 0, 0, 0, 1 };

    UnsignedInt ResultArray[MaxArr];

    for (UnsignedInt p = 0; p < MaxArr; p++)
    {
        MainArrayCopy[p] = MainArray[p];
        ArrayToAddCopy[p] = ArrayToAdd[p];
    }
    for (UnsignedInt m = 1; m <= 16; m++)
    {
        for(UnsignedInt p = 0; p < MaxArr; p++)
            ArrayToAdd[p] = ArrayToAddCopy[p];
        AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);
        ShowArray(ResultArray, 0, MaxArr, "");
        for(UnsignedInt p = 0; p < MaxArr; p++)
            MainArray[p] = ResultArray[p];
    }
}

void Add1ToMainArrayButtonClick()
{
    UnsignedInt MainArrayCopy[MaxArr];
    UnsignedInt ArrayToAddCopy[MaxArr];

    UnsignedInt MainArray[MaxArr] =  { 0,0,0,0, 0, 0, 0, 0 };
    UnsignedInt ArrayToAdd[MaxArr] = { 0,0,0,0, 0, 0, 0, 1 };

    UnsignedInt ResultArray[MaxArr];

    for (UnsignedInt p = 0; p < MaxArr; p++)
    {
        MainArrayCopy[p] = MainArray[p];
        ArrayToAddCopy[p] = ArrayToAdd[p];
    }
    for(UnsignedInt m = 1; m <= 16; m++)
    {
        for (UnsignedInt p = 0; p < MaxArr; p++)
            ArrayToAdd[p] = ArrayToAddCopy[p];
        Add1ToMainArray(MainArray, ArrayToAdd, ResultArray);
        ShowArray(ResultArray, 0, MaxArr, "");
        for (UnsignedInt p = 0; p < MaxArr; p++)
            MainArray[p] = ResultArray[p];
    }
}

void AddTwoDifferentBitsArraysButtonClick()
{
    UnsignedInt MainArrayCopy[MaxArr];
    UnsignedInt ArrayToAddCopy[MaxArr];

    UnsignedInt MainArray[MaxArr] =  {0,0,0,0, 1, 1, 0, 1};
    UnsignedInt ArrayToAdd[MaxArr] = {0,0,0,0, 0, 1, 1, 1};

    UnsignedInt ResultArray[MaxArr];

    for(UnsignedInt p = 0; p < MaxArr; p++)
    {
        MainArrayCopy[p] = MainArray[p];
        ArrayToAddCopy[p] = ArrayToAdd[p];
    }
    AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);

    ShowArray(MainArray, 0, MaxArr, "");
    ShowArray(ArrayToAdd, 0, MaxArr, "");
    ShowArray(ResultArray, 0, MaxArr, "");
}

/*
011101 6
011101 PREV, j + 1 = 3
011101 IN TRAIN PREV, R = 3 i = 3
011101 IN TRAIN AFTER, R = 3 i = 3
011101 IN TRAIN PREV, R = 3 i = 2
001101 IN TRAIN AFTER, R = 3 i = 2
101101 AFTER TOTAL, R = 3, WSTAWIENIE PAM NA MIEJSCE
101101 AFTER

//PONI�SZE DWIE FUNKCJE TO TO SAMO CO DWIE JESZCZE NAST�PNE TYLKO Z WYPISYWANIEM STANU ALGORYTMU przez ShowArray

void  Rotate(UnsignedInt R)
{
    UnsignedInt Pam = ArraySet[R];
    for(UnsignedInt i = R; i > 1; i--)
    {
        ShowArray(ArraySet, 1, N + 1, " IN TRAIN PREV, R = " + UnsignedIntToStr(R) + " i = " + UnsignedIntToStr(i));
        ArraySet[i] = ArraySet[i - 1];
        ShowArray(ArraySet, 1, N + 1, " IN TRAIN AFTER, R = " + UnsignedIntToStr(R) + " i = " + UnsignedIntToStr(i));
    }
    ArraySet[1] = Pam;
    ShowArray(ArraySet, 1, N + 1, " AFTER TOTAL, R = " + UnsignedIntToStr(R));
}

void AllKElementsCombinationsFromNElementsFirstWayButtonClick()
{
    UnsignedInt Line = 1;
    CombinationsMemo->Clear();
    for(UnsignedInt i = 1; i <= N; i++)
        ArraySet[i] = i <= K;
    UnsignedInt j = N - 1;
    bool LoopCondition = true;
    while(LoopCondition)
    {
        ShowArray(ArraySet, 1, N + 1, " " + to_string(Line++));
        if(j >= N || N == K)
        {
            LoopCondition = false;
            goto EndLabel;
        }
        ShowArray(ArraySet, 1, N + 1, " PREV, j + 1 = " + to_string(j + 1));
        Rotate(j + 1);
        ShowArray(ArraySet, 1, N + 1, " AFTER");
        UnsignedInt i;
        for(i = 1; i < N; i++)
            if(ArraySet[i] == 0 && ArraySet[i + 1] == 1)
            {
                j = i + 1;
                break;
            }
        EndLabel:
    }
}
*/

const UnsignedInt N = 6;
const UnsignedInt K = 4;

UnsignedInt ArraySet[N + 1];

void  Rotate(UnsignedInt R)
{
    UnsignedInt Pam = ArraySet[R];
    for(UnsignedInt i = R; i > 1; i--)
        ArraySet[i] = ArraySet[i - 1];
    ArraySet[1] = Pam;
}

void AllKElementsCombinationsFromNElementsFirstWayButtonClick()
{
    UnsignedInt Line = 1;
    for(UnsignedInt i = 1; i <= N; i++)
        ArraySet[i] = i <= K;
    UnsignedInt j = N - 1;
    bool LoopCondition = true;
    while(LoopCondition)
    {
        ShowArray(ArraySet, 1, N + 1, " " + to_string(Line++));
        if(j >= N || N == K)
        {
            LoopCondition = false;
            goto EndLabel;
        }
        Rotate(j + 1);
        UnsignedInt i;
        for(i = 1; i < N; i++)
            if(ArraySet[i] == 0 && ArraySet[i + 1] == 1)
            {
                j = i + 1;
                break;
            }
        EndLabel:;
    }
}

void AllKElementsCombinationsFromNElementsFirstWayAndHalfButtonClick()
{
    UnsignedInt Line = 1;
    for(UnsignedInt i = 1; i <= N; i++)
        ArraySet[i] = i <= K;
    UnsignedInt j = N - 1;
    while(true)
    {
        ShowArray(ArraySet, 1, N + 1, " " + to_string(Line++));
        if(j >= N || N == K)
            break;
        Rotate(j + 1);
        UnsignedInt i;
        for(i = 1; i < N; i++)
            if(ArraySet[i] == 0 && ArraySet[i + 1] == 1)
            {
                j = i + 1;
                break;
            }
    }
}



UnsignedInt NumberOfCombinations(UnsignedInt NP, UnsignedInt KP)
{
    UnsignedInt NSilnia = 1;
    for(UnsignedInt Number = 1; Number <= NP; Number++)
        NSilnia *= Number;
    UnsignedInt KSilnia = 1;
    for(UnsignedInt Number = 1; Number <= KP; Number++)
        KSilnia *= Number;
    UnsignedInt NMinusKSilnia = 1;
    for(UnsignedInt Number = 1; Number <= (NP - KP); Number++)
        NMinusKSilnia *= Number;
    return NSilnia / (KSilnia * NMinusKSilnia);
}

string CreateBoolStringFromInt64BitState(UnsignedInt Number)
{
    string NumberString;
    UnsignedInt LenOfInt = sizeof(UnsignedInt) * 8;
    for(UnsignedInt i = LenOfInt - 1; i >= 1; i--)
        (Number >> (LenOfInt - i - 1)) << (LenOfInt - 1) ? NumberString += "1" : NumberString += "0";
    return NumberString;
}

#pragma warn -ngu
UnsignedInt NextNumberWithTheSameNumberOf1Bits(UnsignedInt Number)
{
    UnsignedInt Smallest, Ripple, Ones;
    Smallest = Number & -Number;
    Ripple = Number + Smallest;
    Ones = Number ^ Ripple;
    Ones = (Ones >> 2) / Smallest;
    return Ripple | Ones;
}

void AllKElementsCombinationsFromNElementsSecondWayButtonClick()
{
    UnsignedInt Line = 1;
    for(UnsignedInt i = 1; i <= N; i++)
        ArraySet[i] = (i <= K);

    UnsignedInt AimNumber = 0;
    for(UnsignedInt i = 1; i <= N; i++)
        AimNumber |= (UnsignedInt)pow(2, i - 1);
    ShowArray(ArraySet, 1, N + 1, "");
    LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(AimNumber)));
    LoggersManagerObject.Log(STREAM(""));

    UnsignedInt Number = 0;
    for(UnsignedInt i = 1; i <= N; i++)
        if(ArraySet[i] == 1)
            Number |= (UnsignedInt)pow(2, i - 1);
    LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) << "    " << to_string(Line++)));
    while(Number < AimNumber)
    {
        Number = NextNumberWithTheSameNumberOf1Bits(Number);
        if(Number < AimNumber)
            LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) << "    " << to_string(Line++)));
    }

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM(to_string(NumberOfCombinations(N,K))));
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//PERMUTATIONS - VARIATIONS

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

void Rotate1(UnsignedInt Number, string& PermutationStringLocal)
{
    UnsignedInt PermutationStringLength = PermutationStringLocal.length();

    char TempChar = PermutationStringLocal[PermutationStringLength];
    for(UnsignedInt i = PermutationStringLength - 1; i >= Number; i--)
        PermutationStringLocal[i + 1] = PermutationStringLocal[i];
    PermutationStringLocal[Number] = TempChar;
}

void PermutationFunction(UnsignedInt Number, string& PermutationStringLocal, UnsignedInt& Line)
{
    if(Number == 1)
        LoggersManagerObject.Log(STREAM(PermutationStringLocal + " LINE = " + to_string(Line++)));
    else
    {
        UnsignedInt PermutationStringLength = PermutationStringLocal.length();

        for(UnsignedInt i = 1; i < Number; i++)
        {
            PermutationFunction(Number - 1, PermutationStringLocal, Line);
            Rotate1(PermutationStringLength - Number + 1 + 1, PermutationStringLocal);
        }
    }
}

void Permutations1ButtonClick()
{
    LoggersManagerObject.Log(STREAM("PERMUTATIONS:"));

    UnsignedInt Line = 1;
    string PermutationString = "ABCDE";
    UnsignedInt PermutationStringLength = PermutationString.length();
    PermutationFunction(PermutationStringLength + 1, PermutationString, Line);
}



















//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//VARIATIONS

void VariationsButtonClick1()
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

void LoopsFrom1ToMInOneLoopButtonClick()
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

void LoopsFromNtoMInOneLoopButtonClick()
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
    LoggersManagerObject.Log(STREAM("T1"));
    Add1ToMainArrayByAddTwoBitsArraysButtonClick();
    LoggersManagerObject.Log(STREAM("T2"));
    Add1ToMainArrayButtonClick();
    LoggersManagerObject.Log(STREAM("T3"));
    AddTwoDifferentBitsArraysButtonClick();

    LoggersManagerObject.Log(STREAM("T4"));
    AllKElementsCombinationsFromNElementsFirstWayButtonClick();
    LoggersManagerObject.Log(STREAM("T5"));
    AllKElementsCombinationsFromNElementsFirstWayAndHalfButtonClick();
    LoggersManagerObject.Log(STREAM("T6"));
    AllKElementsCombinationsFromNElementsSecondWayButtonClick();

    LoggersManagerObject.Log(STREAM("T7"));
    Permutations1ButtonClick();

    LoggersManagerObject.Log(STREAM("T8"));
    VariationsButtonClick1();
    LoggersManagerObject.Log(STREAM("T9"));
    LoopsFrom1ToMInOneLoopButtonClick();
    LoggersManagerObject.Log(STREAM("T10"));
    LoopsFromNtoMInOneLoopButtonClick();
}