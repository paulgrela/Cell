
#include "ExceptionsMacro.h"

#include "Combinatorics.h"

using namespace std;

#define SHOW_EVERY_STEP_OF_GENERATION_COMBINATIONS_

void Combinations::ShowArray(span<UnsignedInt> Array, UnsignedInt StartPos, UnsignedInt Counter, string_view AdditionalText)
{
    try
    {
        string ArrayStr;

        for (UnsignedInt Pos = StartPos; Pos < Counter; Pos++)
            ArrayStr += to_string(Array[Pos]);

        LoggersManagerObject.Log(STREAM(ArrayStr << AdditionalText));
    }
    CATCH("showing arrays")
}

void Combinations::AddTwoBitArrays(const vector<UnsignedInt>& MainArray, const vector<UnsignedInt>& ArrayToAdd, vector<UnsignedInt>& ResultArray) const
{
    try
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
    CATCH("adding two bits arrays")
}

void Combinations::Add1ToMainArray(const vector<UnsignedInt>& MainArray, const vector<UnsignedInt>& ArrayToAdd, vector<UnsignedInt>& ResultArray) const
{
    try
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
    CATCH("adding 1 to main array")
}

void Combinations::Add1ToMainArrayByAddTwoBitsArraysTestDemo(const UnsignedInt NumberOfBits) const
{
    try
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

            AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);

            ShowArray(ResultArray, 0, MaxSizeOfArray, "");

            for (UnsignedInt Pos = 0; Pos < MaxSizeOfArray; Pos++)
                MainArray[Pos] = ResultArray[Pos];
        }
    }
    CATCH("adding 1 to main array by adding two bits arrays test demo")
}

void Combinations::Add1ToMainArrayTestDemo(const UnsignedInt NumberOfBits) const
{
    try
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
    CATCH("adding 1 to main array test demo")
}

void Combinations::AddTwoDifferentBitsArraysTestDemo() const
{
    try
    {
        vector<UnsignedInt> MainArray = { 0, 0, 0, 0, 1, 1, 0, 1 };
        vector<UnsignedInt> ArrayToAdd = { 0, 0, 0, 0, 0, 1, 1, 1 };

        vector<UnsignedInt> ResultArray(MaxSizeOfArray);

        AddTwoBitArrays(MainArray, ArrayToAdd, ResultArray);

        ShowArray(MainArray, 0, MaxSizeOfArray, "");
        ShowArray(ArrayToAdd, 0, MaxSizeOfArray, "");
        ShowArray(ResultArray, 0, MaxSizeOfArray, "");
    }
    CATCH("adding two different bits array test demo")
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
*/
/* PONIZSZE DWIE FUNKCJE TO TO SAMO CO DWIE JESZCZE NASTEPNE TYLKO Z WYPISYWANIEM STANU ALGORYTMU PRZEZ ShowArray */

void Combinations::RotateWithPrint(UnsignedInt R, const UnsignedInt N, const UnsignedInt K)
{
    try
    {
        UnsignedInt Pam = CombinationsArray[R];

        for (UnsignedInt i = R; i > 1; i--)
        {
            #ifdef SHOW_EVERY_STEP_OF_GENERATION_COMBINATIONS
            ShowArray(CombinationsArray, 1, N + 1, " IN TRAIN PREV, R = " + to_string(R) + " i = " + to_string(i));
            #endif

            CombinationsArray[i] = CombinationsArray[i - 1];

            #ifdef SHOW_EVERY_STEP_OF_GENERATION_COMBINATIONS
            ShowArray(CombinationsArray, 1, N + 1, " IN TRAIN AFTER, R = " + to_string(R) + " i = " + to_string(i));
            #endif
        }

        CombinationsArray[1] = Pam;

        #ifdef SHOW_EVERY_STEP_OF_GENERATION_COMBINATIONS
        ShowArray(CombinationsArray, 1, N + 1, " AFTER TOTAL, R = " + to_string(R));
        #endif
    }
    CATCH("rotating for combinations with printing")
}

void Combinations::GenerateAllKElementsCombinationsFromNElementsFirstWayWithPrintTestDemo(UnsignedInt N, UnsignedInt K)
{
    try
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

            #ifdef SHOW_EVERY_STEP_OF_GENERATION_COMBINATIONS
            ShowArray(CombinationsArray, 1, N + 1, " PREV, j + 1 = " + to_string(j + 1));
            #endif

            RotateWithPrint(j + 1, N, K);

            #ifdef SHOW_EVERY_STEP_OF_GENERATION_COMBINATIONS
            ShowArray(CombinationsArray, 1, N + 1, " AFTER");
            #endif

            for (UnsignedInt i = 1; i < N; i++)
                if (CombinationsArray[i] == 0 && CombinationsArray[i + 1] == 1)
                {
                    j = i + 1;
                    break;
                }
            EndLabel:;
        }
    }
    CATCH("generating all k elements combinations from n elements first way with printing")
}

void Combinations::Rotate(UnsignedInt R)
{
    try
    {
        UnsignedInt Pam = CombinationsArray[R];

        for (UnsignedInt i = R; i > 1; i--)
            CombinationsArray[i] = CombinationsArray[i - 1];

        CombinationsArray[1] = Pam;
    }
    CATCH("rotating for combinations")
}

void Combinations::GenerateAllKElementsCombinationsFromNElementsFirstWayTestDemo(UnsignedInt N, UnsignedInt K)
{
    try
    {
        UnsignedInt Line = 1;

        for (UnsignedInt i = 1; i <= N; i++)
            CombinationsArray[i] = (i <= K);

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
    CATCH("generating all k elements combinations from n elements first way")
}

void Combinations::GenerateAllKElementsCombinationsFromNElementsFirstWayAndHalfTestDemo(UnsignedInt N, UnsignedInt K)
{
    try
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
    CATCH("generating all k elements combinations from n elements first way and half")
}

UnsignedInt Combinations::NumberOfCombinations(UnsignedInt NP, UnsignedInt KP)
{
    UnsignedInt NFactorial = 1;
    UnsignedInt KFactorial = 1;
    UnsignedInt NMinusKFactorial = 1;

    try
    {
        for (UnsignedInt Number = 1; Number <= NP; Number++)
            NFactorial *= Number;

        for (UnsignedInt Number = 1; Number <= KP; Number++)
            KFactorial *= Number;

        for (UnsignedInt Number = 1; Number <= (NP - KP); Number++)
            NMinusKFactorial *= Number;
    }
    CATCH("counting number of combinations")

    return NFactorial / (KFactorial * NMinusKFactorial);
}

string Combinations::CreateBoolStringFromInt64BitState(UnsignedInt Number)
{
    string NumberString;

    try
    {
        UnsignedInt LenOfInt = sizeof(UnsignedInt) * 8;
        for (UnsignedInt i = LenOfInt - 1; i >= 1; i--)
            (Number >> (LenOfInt - i - 1)) << (LenOfInt - 1) ? NumberString += "1" : NumberString += "0";
    }
    CATCH("creating bool string from int64 bit state")

    return NumberString;
}

#pragma warn -ngu
UnsignedInt Combinations::NextNumberWithTheSameNumberOf1Bits(UnsignedInt Number)
{
    UnsignedInt Smallest, Ripple, Ones;

    try
    {
        Smallest = Number & -Number;
        Ripple = Number + Smallest;
        Ones = Number ^ Ripple;
        Ones = (Ones >> 2) / Smallest;
    }
    CATCH("next number with the same number of bits")

    return Ripple | Ones;
}

UnsignedInt Combinations::SetAllBitsInNumber(UnsignedInt N, UnsignedInt K)
{
    UnsignedInt Number = 0;

    try
    {
        for (UnsignedInt i = 1; i <= N; i++)
            Number |= (UnsignedInt) pow(2, i - 1);
    }
    CATCH("setting k bits in number")

    return Number;
}

UnsignedInt Combinations::SetKBitsInNumber(UnsignedInt N, UnsignedInt K)
{
    UnsignedInt Number = 0;

    try
    {
        for (UnsignedInt i = 1; i <= N; i++)
            if (i <= K)
                Number |= (UnsignedInt) pow(2, i - 1);
    }
    CATCH("setting k bits in number")

    return Number;
}

void Combinations::GenerateAllKElementsCombinationsFromNElementsSecondWayTestDemo(UnsignedInt N, UnsignedInt K)
{
    try
    {
        UnsignedInt Line = 1;
        UnsignedInt AimNumber = SetAllBitsInNumber(N, K);

        LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(AimNumber)));
        LoggersManagerObject.Log(STREAM(""));

        UnsignedInt Number = SetKBitsInNumber(N, K);

        LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) << "    " << to_string(Line++)));

        while (Number < AimNumber)
        {
            Number = NextNumberWithTheSameNumberOf1Bits(Number);
            if (Number < AimNumber)
                LoggersManagerObject.Log(STREAM(CreateBoolStringFromInt64BitState(Number) << "    " << to_string(Line++)));
        }

        LoggersManagerObject.Log(STREAM(to_string(NumberOfCombinations(N, K))));
    }
    CATCH("generating all k elements combinations from n elements second way")
}

void Permutations::RotateForPermutation(UnsignedInt Number, string& PermutationStringLocal)
{
    try
    {
        UnsignedInt PermutationStringLength = PermutationStringLocal.length();

        char TempChar = PermutationStringLocal[PermutationStringLength];

        for (UnsignedInt i = PermutationStringLength - 1; i >= Number; i--)
            PermutationStringLocal[i + 1] = PermutationStringLocal[i];

        PermutationStringLocal[Number] = TempChar;
    }
    CATCH("rotating for permutation")
}

void Permutations::Permute(UnsignedInt Number, string& PermutationStringLocal, UnsignedInt& Line)
{
    try
    {
        if (Number == 1)
            LoggersManagerObject.Log(STREAM(PermutationStringLocal + " LINE = " + to_string(Line++)));
        else
            for (UnsignedInt i = 1; i < Number; i++)
            {
                Permute(Number - 1, PermutationStringLocal, Line);
                RotateForPermutation(PermutationStringLocal.length() - Number + 1 + 1, PermutationStringLocal);
            }
    }
    CATCH("permuting")
}

void Permutations::PermutationsTestDemo(string& PermutationString)
{
    try
    {
        UnsignedInt Line = 1;
        Permute(PermutationString.length() + 1, PermutationString, Line);
    }
    CATCH("execution permutations test demo")
}

string Variations::GetResultAsString(const UnsignedInt KWords, std::vector<UnsignedInt>& Repeat)
{
    string SKWord;

    try
    {
        for (UnsignedInt ki = 1; ki <= KWords; ki++)
            SKWord += to_string(Repeat[ki]);
    }
    CATCH("getting result as string for variations")

    return SKWord;
}

void Variations::VariationsFrom1ToEndIn3LoopsTestDemo()
{
    try
    {
        LoggersManagerObject.Log(STREAM("VARIATIONS WITH REPETITIONS IN SEVERAL LOOPS WITH STEP = 1:"));

        UnsignedInt Line = 1;
        for (UnsignedInt i1 = 1; i1 <= 4; i1++)
            for (UnsignedInt i2 = 1; i2 <= 4; i2++)
                for (UnsignedInt i3 = 1; i3 <= 4; i3++)
                    LoggersManagerObject.Log(STREAM(to_string(i1) << to_string(i2) << to_string(i3) << " LINE = " + to_string(Line++)));
    }
    CATCH("generating variations from 1 to end in 3 loops");
}

void Variations::VariationsFrom1ToEndIn1LoopTestDemo()
{
    try
    {
        LoggersManagerObject.Log(STREAM("Variations with repetitions - one loop replaces N loops when this 1 loop is n^k times longer:"));
        LoggersManagerObject.Log(STREAM("VARIATIONS WITH REPETITIONS IN SEVERAL LOOP REPLACING SEVERAL LOOPS:"));

        UnsignedInt Line = 1;
        UnsignedInt SizeOfSetN = 4;
        UnsignedInt KWords = 3;
        UnsignedInt TotalNumberOfVariations = pow(SizeOfSetN, KWords);
        vector<UnsignedInt> Repeat(KWords + 1);
        for (UnsignedInt ki = 1; ki <= KWords; ki++)
            Repeat[ki] = 1;
        for (UnsignedInt Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
        {
            LoggersManagerObject.Log(STREAM(GetResultAsString(KWords, Repeat)  << " LINE = " << to_string(Line++)));

            /* to jest zwiekszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
            /* pozostale zwiekszane jedynie gdy poprzednia cyfra jest przewijana */

            Repeat[KWords]++;
            for (UnsignedInt ki = KWords; ki >= 1; ki--)
                if (Repeat[ki] % (SizeOfSetN + 1) == 0)
                {
                    Repeat[ki] = 1;
                    Repeat[ki - 1]++;
                }
        }
    }
    CATCH("generating variations from 1 to end in one loop")
}

void Variations::VariationsFrom1ToMIn3LoopsTestDemo(const vector<UnsignedInt>& Lengths)
{
    try
    {
        LoggersManagerObject.Log(STREAM("VARIATIONS WITH REPETITIONS IN SEVERAL DIFFERENT LOOPS FROM 1 TO LENGTHS [i] WITH STEP = 1"));

        UnsignedInt Line = 1;
        for (UnsignedInt i1 = 1; i1 <= Lengths[1]; i1++)
            for (UnsignedInt i2 = 1; i2 <= Lengths[2]; i2++)
                for (UnsignedInt i3 = 1; i3 <= Lengths[3]; i3++)
                    LoggersManagerObject.Log(STREAM(to_string(i1) << to_string(i2) << to_string(i3) << " LINE = " + to_string(Line++)));
    }
    CATCH("generating variations from 1 to m in 3 loops")
}

void Variations::VariationsFrom1ToMIn1LoopTestDemo(const vector<UnsignedInt>& Lengths)
{
    try
    {
        LoggersManagerObject.Log(STREAM("Variations with repetitions - one loop replaces N loops when this 1 loop is (Length[1] * Length[2] * Length[3] * ... * Length[N]) times longer"));
        LoggersManagerObject.Log(STREAM("VARIATIONS WITH REPETITIONS IN ONE LOOP REPLACING SEVERAL LOOPS:"));

        UnsignedInt Line = 1;
        UnsignedInt KWords = 3;

        UnsignedInt TotalNumberOfVariations = Lengths[1] * Lengths[2] * Lengths[3];
        vector<UnsignedInt> Repeat(KWords + 1);
        for (UnsignedInt ki = 1; ki <= KWords; ki++)
            Repeat[ki] = 1;

        for (UnsignedInt Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
        {
            LoggersManagerObject.Log(STREAM(GetResultAsString(KWords, Repeat)  << " LINE = " << to_string(Line++)));

            /* to jest zwiekszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
            /* pozostale zwiekszane jedynie gdy poprzednia cyfra jest przewijana */

            Repeat[KWords]++;
            for (UnsignedInt ki = KWords; ki >= 1; --ki)
                if (Repeat[ki] % (Lengths[ki] + 1) == 0)
                {
                    Repeat[ki] = 1;
                    Repeat[ki - 1]++;
                }
        }
    }
    CATCH("generating variations from 1 to m in one loop")
}

void Variations::VariationsFromNtoMIn3LoopsTestDemo(const vector<UnsignedInt>& Starts, const vector<UnsignedInt>& Lengths)
{
    try
    {
        LoggersManagerObject.Log(STREAM("VARIATIONS WITH REPETITIONS IN SEVERAL DIFFERENT LOOPS FROM START[i] TO LENGTHS [i] WITH STEP = 1"));

        UnsignedInt Line = 1;
        for (UnsignedInt i1 = Starts[1]; i1 <= Lengths[1]; i1++)
            for (UnsignedInt i2 = Starts[2]; i2 <= Lengths[2]; i2++)
                for (UnsignedInt i3 = Starts[3]; i3 <= Lengths[3]; i3++)
                    LoggersManagerObject.Log(STREAM(to_string(i1) << to_string(i2) << to_string(i3) << " LINE = " + to_string(Line++)));
    }
    CATCH("generating variations from n to m in 3 loops")
}

void Variations::VariationsFromNtoMIn1LoopTestDemo(const vector<UnsignedInt>& Starts, const vector<UnsignedInt>& Lengths)
{
    try
    {
        LoggersManagerObject.Log(STREAM("Variations with repetitions - one loop replaces N loops when this 1 loop is ((Length[1] - Start[1]) * (Length[2] - Start[2]) * (Length[3] - Start[3]) * ... * (Length[N] - Start[N])) times longer:"));
        LoggersManagerObject.Log(STREAM("VARIATIONS WITH REPETITIONS IN SEVERAL LOOP REPLACING SEVERAL LOOPS:"));

        UnsignedInt Line = 1;
        UnsignedInt KWords = 3;

        /* UnsignedInt TotalNumberOfVariations = (Length[1] - Start[1] + 1) * (Length[2] - Start[2] + 1) * (Length[3] - Start[3] + 1); */
        UnsignedInt TotalNumberOfVariations = 1;
        for (UnsignedInt ki = 1; ki <= KWords; ki++)
            TotalNumberOfVariations *= (Lengths[ki] - Starts[ki] + 1);

        vector<UnsignedInt> Repeat(KWords + 1);
        for (UnsignedInt ki = 1; ki <= KWords; ki++)
            Repeat[ki] = Starts[ki];

        for (UnsignedInt Variation = 1; Variation <= TotalNumberOfVariations; Variation++)
        {
            LoggersManagerObject.Log(STREAM(GetResultAsString(KWords, Repeat)  << " LINE = " << to_string(Line++)));

            /* to jest zwiekszane zawsze bo ostatni jest zawsze, ale co pewien czas zerowane */
            /* pozostale zwiekszane jedynie gdy poprzednia cyfra jest przewijana */

            Repeat[KWords]++;
            for (UnsignedInt ki = KWords; ki >= 1; --ki)
                if (Repeat[ki] % (Lengths[ki] + 1) == 0)
                {
                    Repeat[ki] = Starts[ki];
                    Repeat[ki - 1]++;
                }
        }
    }
    CATCH("generating variations from n to m in one loop")
}

void ShowCombinationsDemoTest()
{
    try
    {
        Combinations CombinationsObject;

        LoggersManagerObject.Log(STREAM("T1"));
        CombinationsObject.Add1ToMainArrayByAddTwoBitsArraysTestDemo(16);
        LoggersManagerObject.Log(STREAM("T2"));
        CombinationsObject.Add1ToMainArrayTestDemo(16);
        LoggersManagerObject.Log(STREAM("T3"));
        CombinationsObject.AddTwoDifferentBitsArraysTestDemo();

        LoggersManagerObject.Log(STREAM("T4"));
        CombinationsObject.GenerateAllKElementsCombinationsFromNElementsFirstWayTestDemo(6, 4);
        LoggersManagerObject.Log(STREAM("T5"));
        CombinationsObject.GenerateAllKElementsCombinationsFromNElementsFirstWayAndHalfTestDemo(6, 4);
        LoggersManagerObject.Log(STREAM("T6"));
        CombinationsObject.GenerateAllKElementsCombinationsFromNElementsSecondWayTestDemo(6, 4);
    }
    CATCH("showing combinations demo test")
}

void ShowPermutationsDemoTest()
{
    try
    {
        Permutations PermutationsObject;
        LoggersManagerObject.Log(STREAM("T7"));
        string PermutationString = "ABCDE";
        PermutationsObject.PermutationsTestDemo(PermutationString);
    }
    CATCH("showing permutations demo test")
}

void ShowVariationsDemoTest()
{
    try
    {
        LoggersManagerObject.Log(STREAM("T8"));
        Variations::VariationsFrom1ToEndIn3LoopsTestDemo();
        Variations::VariationsFrom1ToEndIn1LoopTestDemo();
        LoggersManagerObject.Log(STREAM("T9"));
        Variations::VariationsFrom1ToMIn3LoopsTestDemo({ 0, 4, 5, 7 });
        Variations::VariationsFrom1ToMIn1LoopTestDemo({ 0, 4, 5, 7 });
        LoggersManagerObject.Log(STREAM("T10"));
        Variations::VariationsFromNtoMIn3LoopsTestDemo({ 0, 2, 3, 5 }, { 0, 4, 5, 7 });
        Variations::VariationsFromNtoMIn1LoopTestDemo({ 0, 2, 3, 5 },{ 0, 4, 5, 7 });
    }
    CATCH("showing variations demo test")
}

void ShowDemoTest()
{
    try
    {
        ShowCombinationsDemoTest();
        ShowPermutationsDemoTest();
        ShowVariationsDemoTest();
    }
    CATCH("showing combinatorics demo test")
}