#version 450

layout(local_size_x = 256) in;

struct GPUAtom
{
    float X;
    float Y;
    float Z;
    float ColorR;
    float ColorG;
    float ColorB;
    float _padding1[2];
};

struct GPUParticle
{
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

uniform float AtomSize;
uniform float CameraXPosition;
uniform float CameraYPosition;
uniform float CameraZPosition;

uniform vec3 Center;
uniform mat4 ViewMatrix;

struct ParticleOut
{
    mat4 MoveMatrix;
    vec3 Color;
    uint ParticleIdL;
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

        ModelMatrix[0][0] = AtomSize;
        ModelMatrix[1][1] = AtomSize;
        ModelMatrix[2][2] = AtomSize;

        ModelMatrix[3][0] = AtomPosition.x - CameraXPosition - Center.x;
        ModelMatrix[3][1] = AtomPosition.y - CameraYPosition - Center.y;
        ModelMatrix[3][2] = AtomPosition.z - CameraZPosition - Center.z;

        ParticlesOut[AtomOffsetIndexOut].MoveMatrix = ViewMatrix * ModelMatrix;

        ParticlesOut[AtomOffsetIndexOut].Color = vec3(AtomIn.ColorR, AtomIn.ColorG, AtomIn.ColorB);

        ParticlesOut[AtomOffsetIndexOut].ParticleIdL = AtomOffsetIndexOut;
    }
}
