
#include <cmath>
#include <string>

#include "../CellEngineTypes.h"

#include "Logger.h"
#include "Combinatorics.h"


using namespace std;

void ShowArray(int* Array, int Start, int Counter, string LineStr)
{
    string s = "";
    for(int p = Start; p < Counter; p++)
        s += Array[p];
    LoggersManagerObject.Log(STREAM(s + LineStr));
}

const int MaxArr = 8;

/*
void AddTwoBitArrays(int MainArray[max], int ArrayToAdd[max], int ResultArray[max])
{
    int Remember = 0;
    for(int Pos = max - 1; Pos >= 0; Pos--)
    {
        if(Remember == 1 && MainArray[Pos] == 0)
        {
            MainArray[Pos] = 1;
            Remember = 0;
        }
        else
        if(Remember == 1 && MainArray[Pos] == 1 && ArrayToAdd[Pos] == 1)
        {
            ResultArray[Pos] = 1;
            continue;
        }

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
    }
}
*/
void AddTwoBitArrays(int MainArray[MaxArr], int ArrayToAdd[MaxArr], int ResultArray[MaxArr])
{
    int Remember = 0;
    for(int Pos = MaxArr - 1; Pos >= 0; Pos--)
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

void Add1ToMainArray(int MainArray[MaxArr], int ArrayToAdd[MaxArr], int ResultArray[MaxArr])
{
    int Remember = 0;
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

    for(int Pos = MaxArr - 2; Pos >= 0; Pos--)
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
    //CombinationsMemo->Font->Color = clBlue;

    int MainArrayCopy[MaxArr];
    int ArrayToAddCopy[MaxArr];

    int MainArray[MaxArr] =  {0,0,0,0, 0, 0, 0, 0};
    int ArrayToAdd[MaxArr] = {0,0,0,0, 0, 0, 0, 1};

    int ResultArray[MaxArr];

    //CombinationsMemo->Clear();
    for(int p = 0; p < MaxArr; p++)
    {
        MainArrayCopy[p] = MainArray[p];
        ArrayToAddCopy[p] = ArrayToAdd[p];
    }
    for(int m = 1; m <= 16; m++)
    {
        for(int p = 0; p < MaxArr; p++)
            ArrayToAdd[p] = ArrayToAddCopy[p];
        AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);
        ShowArray(ResultArray, 0, MaxArr, "");
        for(int p = 0; p < MaxArr; p++)
            MainArray[p] = ResultArray[p];
    }
}

void Add1ToMainArrayButtonClick()
{
    //CombinationsMemo->Font->Color = clPurple;

    int MainArrayCopy[MaxArr];
    int ArrayToAddCopy[MaxArr];

    int MainArray[MaxArr] =  {0,0,0,0, 0, 0, 0, 0};
    int ArrayToAdd[MaxArr] = {0,0,0,0, 0, 0, 0, 1};

    int ResultArray[MaxArr];

    //CombinationsMemo->Clear();
    for(int p = 0; p < MaxArr; p++)
    {
        MainArrayCopy[p] = MainArray[p];
        ArrayToAddCopy[p] = ArrayToAdd[p];
    }
    for(int m = 1; m <= 16; m++)
    {
        for(int p = 0; p < MaxArr; p++)
            ArrayToAdd[p] = ArrayToAddCopy[p];
        Add1ToMainArray(MainArray, ArrayToAdd, ResultArray);
        ShowArray(ResultArray, 0, MaxArr, "");
        for(int p = 0; p < MaxArr; p++)
            MainArray[p] = ResultArray[p];
    }
}

void AddTwoDifferentBitsArraysButtonClick()
{
    //CombinationsMemo->Font->Color = clBlack;

    int MainArrayCopy[MaxArr];
    int ArrayToAddCopy[MaxArr];

    int MainArray[MaxArr] =  {0,0,0,0, 1, 1, 0, 1};
    int ArrayToAdd[MaxArr] = {0,0,0,0, 0, 1, 1, 1};

    int ResultArray[MaxArr];

    //CombinationsMemo->Clear();
    for(int p = 0; p < MaxArr; p++)
    {
        MainArrayCopy[p] = MainArray[p];
        ArrayToAddCopy[p] = ArrayToAdd[p];
    }
    AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);
    for(int i = 1; i <= 3; i++)
    {
        string s = "";
        for(int p = 0; p < MaxArr; p++)
        {
            if(i == 1)
                s += MainArrayCopy[p];
            else
            if(i == 2)
                s += ArrayToAddCopy[p];
            else
            if(i == 3)
                s += ResultArray[p];
        }
        LoggersManagerObject.Log(STREAM(s));
    }
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

void  Rotate(int R)
{
    int Pam = ArraySet[R];
    for(int i = R; i > 1; i--)
    {
        ShowArray(ArraySet, 1, N + 1, " IN TRAIN PREV, R = " + IntToStr(R) + " i = " + IntToStr(i));
        ArraySet[i] = ArraySet[i - 1];
        ShowArray(ArraySet, 1, N + 1, " IN TRAIN AFTER, R = " + IntToStr(R) + " i = " + IntToStr(i));
    }
    ArraySet[1] = Pam;
    ShowArray(ArraySet, 1, N + 1, " AFTER TOTAL, R = " + IntToStr(R));
}

void AllKElementsCombinationsFromNElementsFirstWayButtonClick()
{
    int Line = 1;
    CombinationsMemo->Clear();
    for(int i = 1; i <= N; i++)
        ArraySet[i] = i <= K;
    int j = N - 1;
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
        int i;
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

const int N = 6;
const int K = 4;

int ArraySet[N + 1];

void  Rotate(int R)
{
    int Pam = ArraySet[R];
    for(int i = R; i > 1; i--)
        ArraySet[i] = ArraySet[i - 1];
    ArraySet[1] = Pam;
}

void AllKElementsCombinationsFromNElementsFirstWayButtonClick()
{
    int Line = 1;
    //CombinationsMemo->Clear();
    for(int i = 1; i <= N; i++)
        ArraySet[i] = i <= K;
    int j = N - 1;
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
        int i;
        for(i = 1; i < N; i++)
            if(ArraySet[i] == 0 && ArraySet[i + 1] == 1)
            {
                j = i + 1;
                break;
            }
        EndLabel:;
    }
    //EndLabel:;
}

void AllKElementsCombinationsFromNElementsFirstWayAndHalfButtonClick()
{
    int Line = 1;
    //CombinationsMemo->Clear();
    for(int i = 1; i <= N; i++)
        ArraySet[i] = i <= K;
    int j = N - 1;
    while(true)
    {
        ShowArray(ArraySet, 1, N + 1, " " + to_string(Line++));
        if(j >= N || N == K)
            break;
        Rotate(j + 1);
        int i;
        for(i = 1; i < N; i++)
            if(ArraySet[i] == 0 && ArraySet[i + 1] == 1)
            {
                j = i + 1;
                break;
            }
    }
}

UnsignedInt NumberOfCombinations(UnsignedInt N, UnsignedInt K)
{
    UnsignedInt NSilnia = 1;
    for(int Number = 1; Number <= N; Number++)
        NSilnia *= Number;
    UnsignedInt KSilnia = 1;
    for(int Number = 1; Number <= K; Number++)
        KSilnia *= Number;
    UnsignedInt NMinusKSilnia = 1;
    for(int Number = 1; Number <= (N - K); Number++)
        NMinusKSilnia *= Number;
    return NSilnia / (KSilnia * NMinusKSilnia);
}

string CreateBoolStringFromInt64BitState(UnsignedInt Number)
{
    string NumberString = "";
    unsigned int LenOfInt = sizeof(UnsignedInt) * 8;
    for(unsigned int i = LenOfInt - 1; i >= 1; i--)
        (Number >> (LenOfInt - i - 1)) << (LenOfInt - 1) ? NumberString += "1" : NumberString += "0";
    return NumberString;
}

#pragma warn -ngu
unsigned int NextNumberWithTheSameNumberOf1Bits(unsigned int Number)
{
    unsigned int Smallest, Ripple, Ones;
    Smallest = Number & -Number;
    Ripple = Number + Smallest;
    Ones = Number ^ Ripple;
    Ones = (Ones >> 2) / Smallest;
    return Ripple | Ones;
}

void AllKElementsCombinationsFromNElementsSecondWayButtonClick()
{
    int Line = 1;
    //CombinationsMemo->Clear();
    for(int i = 1; i <= N; i++)
        ArraySet[i] = i <= K;

    unsigned int AimNumber = 0;
    for(int i = 1; i <= N; i++)
        AimNumber |= (UnsignedInt)pow(2, i - 1);
    ShowArray(ArraySet, 1, N + 1, "");
    LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(AimNumber)));
    LoggersManagerObject.Log(STREAM(""));

    unsigned int Number = 0;
    for(int i = 1; i <= N; i++)
        if(ArraySet[i] == 1)
            Number |= (UnsignedInt)pow(2, i - 1);
    LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) + "    " + to_string(Line++)));
    while(Number < AimNumber)
    {
        Number = NextNumberWithTheSameNumberOf1Bits(Number);
        if(Number < AimNumber)
            LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) + "    " + to_string(Line++)));
    }

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM(to_string(NumberOfCombinations(N,K))));
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//PERMUTATIONS






















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

void Rotate(int Number, string& PermutationStringLocal)
{
    int PermutationStringLength = PermutationStringLocal.length();

    char TempChar = PermutationStringLocal[PermutationStringLength];
    for(int i = PermutationStringLength - 1; i >= Number; i--)
        PermutationStringLocal[i + 1] = PermutationStringLocal[i];
    PermutationStringLocal[Number] = TempChar;
}

void PermutationFunction(int Number, string& PermutationStringLocal, int& Line)
{
    if(Number == 1)
        LoggersManagerObject.Log(STREAM(PermutationStringLocal + " LINE = " + to_string(Line++)));
    else
    {
        int PermutationStringLength = PermutationStringLocal.length();

        for(int i = 1; i < Number; i++)
        {
            PermutationFunction(Number - 1, PermutationStringLocal, Line);
            Rotate(PermutationStringLength - Number + 1 + 1, PermutationStringLocal);
        }
    }
}

void Permut1ButtonClick()
{
    //PermutationsTForm->PermutationsMemo->Lines->Clear();
    int Line = 1;
    string PermutationString = "ABCDE";
    int PermutationStringLength = PermutationString.length();
    PermutationFunction(PermutationStringLength + 1, PermutationString, Line);
}


void VariationsButtonClick()
{
    /* To s¹ wariacje z powtórzeniami - ka¿dy k wyrazowy ci¹g ze zbioru n elementowego czyli n do k-tej */
    /* czyli tu 4 do 3 czyli k = 3 a n = 4 - mam 4 cyfry - to zbiór n = 4 i mam wyrazy 3 elementowe czyli k = 3 */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("W KILKU PÊTLACH:"));
    LoggersManagerObject.Log(STREAM(""));

    int Line = 1;
    for(int i1 = 1; i1 <= 4; i1++)
        for(int i2 = 1; i2 <= 4; i2++)
            for(int i3 = 1; i3 <= 4; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) + to_string(i2) + to_string(i3) + " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest n^k razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("W JEDNEJ PÊTLI:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    int SizeOfSetN = 4;
    int KWords = 3;
    int TotalNumberOfVariations = pow(SizeOfSetN, KWords);
    std::vector<int> Repeat(KWords + 1);
    for(int ki = 1; ki <= KWords; ki++)
        Repeat[ki] = 1;
    for(int Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        string SKWord = "";
        for(int ki = 1; ki <= KWords; ki++)
            SKWord += Repeat[ki];
        LoggersManagerObject.Log(STREAM(SKWord + " LINE = " + to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzedni jest przewijany */

        Repeat[KWords]++;
        for(int ki = KWords; ki >= 1; ki--)
            if(Repeat[ki] % (SizeOfSetN + 1) == 0)
            {
                Repeat[ki] = 1;
                Repeat[ki - 1]++;
            }
    }
}















//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//PERMUTATIONS - VARIATIONS

#include <cmath>
#include <string>
#include <vector>
//#include "PermutationsVariationsLoopsForm.h"

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

void Rotate1(int Number, string& PermutationStringLocal)
{
    int PermutationStringLength = PermutationStringLocal.length();

    char TempChar = PermutationStringLocal[PermutationStringLength];
    for(int i = PermutationStringLength - 1; i >= Number; i--)
        PermutationStringLocal[i + 1] = PermutationStringLocal[i];
    PermutationStringLocal[Number] = TempChar;
}

void PermutationFunction1(int Number, string& PermutationStringLocal, int& Line)
{
    if(Number == 1)
        LoggersManagerObject.Log(STREAM(PermutationStringLocal + " LINE = " + to_string(Line++)));
    else
    {
        int PermutationStringLength = PermutationStringLocal.length();

        for(int i = 1; i < Number; i++)
        {
            PermutationFunction1(Number - 1, PermutationStringLocal, Line);
            Rotate1(PermutationStringLength - Number + 1 + 1, PermutationStringLocal);
        }
    }
}

void Permutations1ButtonClick()
{
    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("PERMUTACJE:"));
    LoggersManagerObject.Log(STREAM(""));

    int Line = 1;
    string PermutationString = "ABCDE";
    int PermutationStringLength = PermutationString.length();
    PermutationFunction1(PermutationStringLength + 1, PermutationString, Line);
}


void VariationsButtonClick1()
{
    /* To s¹ wariacje z powtórzeniami - ka¿dy k wyrazowy ci¹g ze zbioru n elementowego czyli n do k-tej */
    /* czyli tu 4 do 3 czyli k = 3 a n = 4 - mam 4 cyfry - to zbiór n = 4 i mam wyrazy 3 elementowe czyli k = 3 */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("WARIACJE Z POWTÓRZENIAMI W KILKU PÊTLACH:"));
    LoggersManagerObject.Log(STREAM(""));

    int Line = 1;
    for(int i1 = 1; i1 <= 4; i1++)
        for(int i2 = 1; i2 <= 4; i2++)
            for(int i3 = 1; i3 <= 4; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) + to_string(i2) + to_string(i3) + " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest n^k razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("WARIACJE Z POWTÓRZENIAMI W JEDNEJ PÊTLI ZAMIAST KILKU PÊTLI:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    int SizeOfSetN = 4;
    int KWords = 3;
    int TotalNumberOfVariations = pow(SizeOfSetN, KWords);
    std::vector<int> Repeat(KWords + 1);
    for(int ki = 1; ki <= KWords; ki++)
        Repeat[ki] = 1;
    for(int Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        /* Poni¿sze linijki to jedynie wyœwietlenie rezultatu */

        string SKWord = "";
        for(int ki = 1; ki <= KWords; ki++)
            SKWord += Repeat[ki];
        LoggersManagerObject.Log(STREAM(SKWord + " LINE = " + to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzednia cyfra jest przewijana */

        Repeat[KWords]++;
        for(int ki = KWords; ki >= 1; ki--)
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

    int Line = 1;
    int Length[4] = {0, 4, 5, 7};

    for(int i1 = 1; i1 <= Length[1]; i1++)
        for(int i2 = 1; i2 <= Length[2]; i2++)
            for(int i3 = 1; i3 <= Length[3]; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) + to_string(i2) + to_string(i3) + " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest (Length[1] * Length[2] * Length[3]) razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("POPRZEDNIE PÊTLE W JEDNEJ PÊTLI ZAMIAST KILKU:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    int KWords = 3;
    int TotalNumberOfVariations = Length[1] * Length[2] * Length[3];
    std::vector<int> Repeat(KWords + 1);
    for(int ki = 1; ki <= KWords; ki++)
        Repeat[ki] = 1;
    for(int Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        /* Poni¿sze linijki to jedynie wyœwietlenie rezultatu */

        string SKWord = "";
        for(int ki = 1; ki <= KWords; ki++)
            SKWord += Repeat[ki];
        LoggersManagerObject.Log(STREAM(SKWord + " LINE = " + to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzednia cyfra jest przewijana */

        Repeat[KWords]++;
        for(int ki = KWords; ki >= 1; --ki)
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

    int Line = 1;
    int Length[4] = {0, 4, 5, 7};
    int Start[4] = {0, 2, 3, 5};

    for(int i1 = Start[1]; i1 <= Length[1]; i1++)
        for(int i2 = Start[2]; i2 <= Length[2]; i2++)
            for(int i3 = Start[3]; i3 <= Length[3]; i3++)
                LoggersManagerObject.Log(STREAM(to_string(i1) + to_string(i2) + to_string(i3) + " LINE = " +to_string(Line++)));

    /* To samo poni¿ej - jedna pêtla zastêpuje 3 gdy jest (Length[1] * Length[2] * Length[3]) razy d³u¿sza */

    LoggersManagerObject.Log(STREAM(""));
    LoggersManagerObject.Log(STREAM("POPRZEDNIE PÊTLE W JEDNEJ PÊTLI ZAMIAST KILKU:"));
    LoggersManagerObject.Log(STREAM(""));

    Line = 1;
    int KWords = 3;

    /* int TotalNumberOfVariations = (Length[1] - Start[1] + 1) * (Length[2] - Start[2] + 1) * (Length[3] - Start[3] + 1); */
    int TotalNumberOfVariations = 1;
    for(int ki = 1; ki <= KWords; ki++)
        TotalNumberOfVariations *= (Length[ki] - Start[ki] + 1);

    std::vector<int> Repeat(KWords + 1);
    for(int ki = 1; ki <= KWords; ki++)
        Repeat[ki] = Start[ki];
    for(int Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
    {
        /* Poni¿sze linijki to jedynie wyœwietlenie rezultatu */

        string SKWord = "";
        for(int ki = 1; ki <= KWords; ki++)
            SKWord += Repeat[ki];
        LoggersManagerObject.Log(STREAM(SKWord + " LINE = " + to_string(Line++)));

        /* to jest zwiêkszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
        /* pozosta³e zwiêkszane jedynie gdy poprzednia cyfra jest przewijana */

        Repeat[KWords]++;
        for(int ki = KWords; ki >= 1; --ki)
            if(Repeat[ki] % (Length[ki] + 1) == 0)
            {
                Repeat[ki] = Start[ki];
                Repeat[ki - 1]++;
            }
    }
}

