
#ifndef __SB7TEXTOVERLAY_H__
#define __SB7TEXTOVERLAY_H__

#include <sb7.h>

namespace sb7
{
    class text_overlay
    {
    public:
        text_overlay() : cursor_x(0), cursor_y(0)
        {
        }

        void init(int width, int height, const char* font = nullptr);
        void teardown();
        void draw();

        void drawText(const char* str, int x, int y);
        void print(const char* str);
        void scroll(int lines);
        void moveCursor(int x, int y);
        void clear();

    private:
        GLuint text_buffer;
        GLuint font_texture;
        GLuint vao;

        GLuint text_program;
        char *screen_buffer;
        int buffer_width;
        int buffer_height;
        bool dirty;
        int cursor_x;
        int cursor_y;
    };
}

#endif
