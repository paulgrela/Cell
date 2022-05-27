
#ifndef __SB6MFILE_H__
#define __SB6MFILE_H__

#define SB6M_FOURCC(a,b,c,d) ( ((unsigned int)(a) << 0) | ((unsigned int)(b) << 8) | ((unsigned int)(c) << 16) | ((unsigned int)(d) << 24) )
#define SB6M_MAGIC SB6M_FOURCC('S','B','6','M')

#ifdef _MSC_VER
#pragma pack (push, 1)
#endif

typedef enum SB6M_CHUNK_TYPE_t
{
    SB6M_CHUNK_TYPE_INDEX_DATA = SB6M_FOURCC('I','N','D','X'),
    SB6M_CHUNK_TYPE_VERTEX_DATA = SB6M_FOURCC('V','R','T','X'),
    SB6M_CHUNK_TYPE_VERTEX_ATTRIBS = SB6M_FOURCC('A','T','R','B'),
    SB6M_CHUNK_TYPE_SUB_OBJECT_LIST = SB6M_FOURCC('O','L','S','T'),
    SB6M_CHUNK_TYPE_COMMENT = SB6M_FOURCC('C','M','N','T'),
    SB6M_CHUNK_TYPE_DATA = SB6M_FOURCC('D','A','T','A')
} 
SB6M_CHUNK_TYPE;

typedef struct SB6M_HEADER_t
{
    union
    {
        unsigned int Magic;
        char MagicName[4];
    };
    unsigned int Size;
    unsigned int NumberOfChunks;
    unsigned int Flags;
} 
SB6M_HEADER;

typedef struct SB6M_CHUNK_HEADER_t
{
    union
    {
        unsigned int ChunkType;
        char ChunkName[4];
    };
    unsigned int Size;
} 
SB6M_CHUNK_HEADER;

typedef struct SB6M_CHUNK_INDEX_DATA_t
{
    SB6M_CHUNK_HEADER Header;
    unsigned int IndexType;
    unsigned int IndexCount;
    unsigned int IndexDataOffset;
} 
SB6M_CHUNK_INDEX_DATA;

typedef struct SB6M_CHUNK_VERTEX_DATA_t
{
    SB6M_CHUNK_HEADER Header;
    unsigned int DataSize;
    unsigned int DataOffset;
    unsigned int TotalVertices;
} 
SB6M_CHUNK_VERTEX_DATA;

typedef struct SB6M_VERTEX_ATTRIB_DECL_t
{
    char Name[64];
    unsigned int Size;
    unsigned int Type;
    unsigned int Stride;
    unsigned int Flags;
    unsigned int DataOffset;
} 
SB6M_VERTEX_ATTRIB_DECL;

#define SB6M_VERTEX_ATTRIB_FLAG_NORMALIZED      0x00000001
#define SB6M_VERTEX_ATTRIB_FLAG_INTEGER         0x00000002

typedef struct SB6M_VERTEX_ATTRIB_CHUNK_t
{
    SB6M_CHUNK_HEADER Header;
    unsigned int AttribCount;
    SB6M_VERTEX_ATTRIB_DECL AttribData[1];
} 
SB6M_VERTEX_ATTRIB_CHUNK;

typedef enum SB6M_DATA_ENCODING_t
{
    SB6M_DATA_ENCODING_RAW = 0
} 
SB6M_DATA_ENCODING;

typedef struct SB6M_DATA_CHUNK_t
{
    SB6M_CHUNK_HEADER Header;
    unsigned int Encoding;
    unsigned int DataOffset;
    unsigned int DataLength;
} 
SB6M_DATA_CHUNK;

typedef struct SB6M_SUB_OBJECT_DECL_t
{
    unsigned int First;
    unsigned int Count;
} 
SB6M_SUB_OBJECT_DECL;

typedef struct SB6M_CHUNK_SUB_OBJECT_LIST_t
{
    SB6M_CHUNK_HEADER Header;
    unsigned int Count;
    SB6M_SUB_OBJECT_DECL SubObject[1];
} 
SB6M_CHUNK_SUB_OBJECT_LIST;

typedef struct SB6M_CHUNK_COMMENT_t
{
    SB6M_CHUNK_HEADER Header;
    char Comment[1];
} 
SB6M_CHUNK_COMMENT;

#ifdef _MSC_VER
#pragma pack (pop)
#endif

#endif
