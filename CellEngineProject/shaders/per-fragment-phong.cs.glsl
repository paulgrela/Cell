#version 450

//INPUT DATA

layout(local_size_x = 256) in;

struct GPUAtom
{
    float X;
    float Y;
    float Z;
    //float ColorR;
    //float ColorG;
    //float ColorB;

    //short ColorR;
    //short ColorG;
    //short ColorB;
    uint ColorR;
    uint ColorG;
    uint ColorB;
    uint _padding[2];

    //float _padding1[2];
    //vec2 _padding1;
};

struct GPUParticle
{
    uint EntityId;
    uint ChainId;
    uint Index;
    uint AtomOffset;
    uint AtomCount;
    uint _padding[3];
};

layout(std430, binding = 0) readonly buffer ParticleBufferIn
{
    GPUParticle ParticlesIn[];
};

layout(std430, binding = 1) readonly buffer AtomBufferIn
{
    GPUAtom AtomsIn[];
};

uniform vec3 Center;
uniform mat4 ViewMatrix;

//OUT DATA

struct ParticleOut
{
    mat4 MoveMatrix;
    vec3 Color;
    float padding;
};

layout(std430, binding = 2) buffer ParticlesBufferOut
{
    ParticleOut ParticlesOut[];
};

void main()
{
    uint ParticleId = gl_GlobalInvocationID.x;
    if (ParticleId >= ParticlesIn.length())
        return;

    GPUParticle ParticleIn = ParticlesIn[ParticleId];

    for (uint AtomIndex = 0; AtomIndex < ParticleIn.AtomCount; AtomIndex++)
    {
        uint AtomOffsetIndexOut = ParticleIn.AtomOffset + AtomIndex;
        GPUAtom AtomIn = AtomsIn[AtomOffsetIndexOut];

        vec3 AtomPosition = vec3(AtomIn.X, AtomIn.Y, AtomIn.Z);

        mat4 ModelMatrix = mat4(1.0);
        ModelMatrix[3][0] = AtomPosition.x - Center.x;
        ModelMatrix[3][1] = AtomPosition.y - Center.y;
        ModelMatrix[3][2] = AtomPosition.z - Center.z;

        ParticlesOut[AtomOffsetIndexOut].MoveMatrix = ViewMatrix * ModelMatrix;
        //ParticlesOut[AtomOffsetIndexOut].Color = vec4(AtomIn.ColorR / 255.0, AtomIn.ColorG / 255.0, AtomIn.ColorB, 1.0 / 255.0);
        ParticlesOut[AtomOffsetIndexOut].Color = vec3(AtomIn.ColorR, AtomIn.ColorG, AtomIn.ColorB);
    }
}
