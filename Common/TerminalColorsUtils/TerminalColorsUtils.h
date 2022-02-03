#pragma once

#ifndef _TERMINAL_COLORS_UTILS_H_
#define _TERMINAL_COLORS_UTILS_H_

#include <iostream>

#ifdef WINDOWS_PLATFORM

#include <windows.h>

namespace terminal_colors_utils
{
    inline std::ostream& red(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_INTENSITY);
        return s;
    }

    inline std::ostream& blue(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
        return s;
    }

    inline std::ostream& green(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN | FOREGROUND_INTENSITY);
        return s;
    }

    inline std::ostream& yellow(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY);
        return s;
    }

    inline std::ostream& white(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
        return s;
    }

    inline std::ostream& black(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, 0 | FOREGROUND_INTENSITY);
        return s;
    }

    inline std::ostream& background_red(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, BACKGROUND_RED);
        return s;
    }

    inline std::ostream& background_blue(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, BACKGROUND_BLUE);
        return s;
    }

    inline std::ostream& background_green(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, BACKGROUND_GREEN);
        return s;
    }

    inline std::ostream& background_yellow(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, BACKGROUND_GREEN | BACKGROUND_RED);
        return s;
    }

    inline std::ostream& background_white(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE);
        return s;
    }

    inline std::ostream& background_black(std::ostream& s)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, 0);
        return s;
    }

    struct TerminalColor
    {
        TerminalColor(WORD attribute) :m_color(attribute) {};
        WORD m_color;
    };

    template <class _Elem, class _Traits>
    std::basic_ostream<_Elem, _Traits>& operator<<(std::basic_ostream<_Elem, _Traits>& i, TerminalColor& c)
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(hStdout, c.m_color);
        return i;
    }
};

#endif

#ifdef UNIX_PLATFORM

namespace terminal_colors_utils
{
    inline std::ostream& white(std::ostream& s)
    {
        return s << "\033[0m";
    }

    inline std::ostream& magneta(std::ostream& s)
    {
        return s << "\033[0;35m";
    }

    inline std::ostream& cyan(std::ostream& s)
    {
        return s << "\033[0;36m";
    }

    inline std::ostream& blue(std::ostream& s)
    {
        return s << "\033[0;34m";
    }

    inline std::ostream& green(std::ostream& s)
    {
        return s << "\033[1;32m";
    }

    inline std::ostream& red(std::ostream& s)
    {
        return s << "\033[0;31m";
    }

    inline std::ostream& yellow(std::ostream& s)
    {
        return s << "\033[1;33m";
    }

    inline std::ostream& background_white(std::ostream& s)
    {
        return s << "\033[0m";
    }

    inline std::ostream& background_black(std::ostream& s)
    {
        return s << "\033[0m";
    }
};

#endif

#endif
