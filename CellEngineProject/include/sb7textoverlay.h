
#ifndef __SB7TEXTOVERLAY_H__
#define __SB7TEXTOVERLAY_H__

#include <sb7.h>

namespace sb7
{
    class TextOverlay
    {
    public:
        TextOverlay() : CursorX(0), CursorY(0)
        {
        }

        void Init(int width, int height, const char* font = nullptr);
        void TearDown();
        void Draw();

        void DrawText(const char* StringParameter, int X, int Y);
        void Print(const char* StringParameter);
        void Scroll(int Lines);
        void MoveCursor(int x, int y);
        void Clear();

    private:
        GLuint TextBuffer;
        GLuint FontTexture;
        GLuint VAO;

        GLuint TextProgram;
        char* ScreenBuffer;
        int BufferWidth;
        int BufferHeight;
        bool Dirty;
        int CursorX;
        int CursorY;
    };
}

#endif
