
#ifndef SB6KTX_H
#define SB6KTX_H

namespace sb7::ktx::File
{
    struct Header
    {
        unsigned char Identifier[12];
        unsigned int Endianness;
        unsigned int glType;
        unsigned int glTypeSize;
        unsigned int glFormat;
        unsigned int glInternalFormat;
        unsigned int glBaseInternalFormat;
        unsigned int PixelWidth;
        unsigned int PixelHeight;
        unsigned int PixelDepth;
        unsigned int ArrayElements;
        unsigned int Faces;
        unsigned int mipLevels;
        unsigned int KeyPairBytes;
    };

    union KeyValuePair
    {
        unsigned int Size;
        unsigned char RawBytes[4];
    };

    unsigned int Load(const char * FileName, unsigned int Tex = 0);
    bool Save(const char * FileName, unsigned int Target, unsigned int Tex);
}

#endif
