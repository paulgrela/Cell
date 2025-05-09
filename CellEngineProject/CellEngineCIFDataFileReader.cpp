
#include <regex>
#include <fstream>

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include "ExceptionsMacro.h"

#include "StringUtils.h"

#include "CellEngineUseful.h"

#include "CellEngineConstants.h"
#include "CellEngineConfigData.h"
#include "CellEngineParticlesKindsManager.h"

#include "CellEngineCIFDataFileReader.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

constexpr bool WriteErrorsInCIFFileToScreen = false;

using namespace std;
using namespace string_utils;

CellEngineAtom CellEngineCIFDataFileReader::ParseRecord(const char* LocalCIFRecord)
{
    CellEngineAtom CellEngineAtomObject{};

    try
    {
        string RecordStr = LocalCIFRecord;

        vector<string> AtomFields = split(RecordStr, " ");

        strncpy(CellEngineAtomObject.Name, AtomFields[3].c_str(), AtomFields[3].length() + 1);
        strncpy(CellEngineAtomObject.ResName, AtomFields[5].c_str(), AtomFields[5].length() + 1);
        strncpy(CellEngineAtomObject.Chain, AtomFields[6].c_str(), AtomFields[6].length() + 1);
        CellEngineAtomObject.EntityId = stoi(AtomFields[7]);

        CellEngineAtomObject.SetAtomPositionsData(stof(AtomFields[10]), stof(AtomFields[11]), stof(AtomFields[12]));

        auto AtomKindObjectIterator = ParticlesKindsManagerObject.GetGraphicAtomKindDataFromAtomName(AtomFields[3][0]);
        CellEngineAtomObject.AtomColor = AtomKindObjectIterator->Color;
        #ifdef EXTENDED_RAM_MEMORY
        CellEngineAtomObject.SizeXAtom = AtomKindObjectIterator->SizeX;
        CellEngineAtomObject.SizeYAtom = AtomKindObjectIterator->SizeY;
        CellEngineAtomObject.SizeZAtom = AtomKindObjectIterator->SizeZ;
        #endif

        #ifdef EXTENDED_RAM_MEMORY
        auto ParticleKindObject = ParticlesKindsManagerObject.GetGraphicParticleKind(CellEngineAtomObject.EntityId);
        CellEngineAtomObject.SizeXParticle = ParticleKindObject.SizeX;
        CellEngineAtomObject.SizeYParticle = ParticleKindObject.SizeY;
        CellEngineAtomObject.SizeZParticle = ParticleKindObject.SizeZ;
        #endif
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

void CopyFieldsToMatrix(TransformationMatrix3x4& TransformationMatrix3x4Object, const vector<string>& MatrixFields)
{
    try
    {
        constexpr UnsignedInt Shift = 6;
        TransformationMatrix3x4Object.Matrix[0][0] = stof(MatrixFields[Shift + 0]);
        TransformationMatrix3x4Object.Matrix[0][1] = stof(MatrixFields[Shift + 1]);
        TransformationMatrix3x4Object.Matrix[0][2] = stof(MatrixFields[Shift + 2]);
        TransformationMatrix3x4Object.Matrix[0][3] = stof(MatrixFields[Shift + 3]);
        TransformationMatrix3x4Object.Matrix[1][0] = stof(MatrixFields[Shift + 4]);
        TransformationMatrix3x4Object.Matrix[1][1] = stof(MatrixFields[Shift + 5]);
        TransformationMatrix3x4Object.Matrix[1][2] = stof(MatrixFields[Shift + 6]);
        TransformationMatrix3x4Object.Matrix[1][3] = stof(MatrixFields[Shift + 7]);
        TransformationMatrix3x4Object.Matrix[2][0] = stof(MatrixFields[Shift + 8]);
        TransformationMatrix3x4Object.Matrix[2][1] = stof(MatrixFields[Shift + 9]);
        TransformationMatrix3x4Object.Matrix[2][2] = stof(MatrixFields[Shift + 10]);
        TransformationMatrix3x4Object.Matrix[2][3] = stof(MatrixFields[Shift + 11]);
    }
    CATCH("copying fields to matrix")
}

glm::vec3 CountResultPositionsFromTransformationMatrix(const std::unordered_map<UnsignedInt, TransformationMatrix3x4>::iterator& TransformationMatrixIterator, const CellEngineAtom& AppliedAtom)
{
    glm::vec3 Result{};

    try
    {
        auto TransformationMatrix = glm::mat3();

        TransformationMatrix[0][0] = TransformationMatrixIterator->second.Matrix[0][0];
        TransformationMatrix[0][1] = TransformationMatrixIterator->second.Matrix[1][0];
        TransformationMatrix[0][2] = TransformationMatrixIterator->second.Matrix[2][0];

        TransformationMatrix[1][0] = TransformationMatrixIterator->second.Matrix[0][1];
        TransformationMatrix[1][1] = TransformationMatrixIterator->second.Matrix[1][1];
        TransformationMatrix[1][2] = TransformationMatrixIterator->second.Matrix[2][1];

        TransformationMatrix[2][0] = TransformationMatrixIterator->second.Matrix[0][2];
        TransformationMatrix[2][1] = TransformationMatrixIterator->second.Matrix[1][2];
        TransformationMatrix[2][2] = TransformationMatrixIterator->second.Matrix[2][2];

        Result = TransformationMatrix * glm::vec3(AppliedAtom.X, AppliedAtom.Y, AppliedAtom.Z) + glm::vec3(TransformationMatrixIterator->second.Matrix[0][3], TransformationMatrixIterator->second.Matrix[1][3], TransformationMatrixIterator->second.Matrix[2][3]);
    }
    CATCH("copying from matrix to matrix")

    return Result;
}

void AssociateAutinNameWithIllinoisName(unordered_map<string, string>& AutinIllinoisNamesMap)
{
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.DNA'"] = "DNANucleotide";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.RNA'"] = "RNANucleotide";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.DNA_POL'"] = "particle_DNAPolymerase";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.RNA_POL'"] = "particle_RNAPolymerase";
    AutinIllinoisNamesMap["'root.syn3A.membrane.inner_membrane'"] = "LipidInnerMembrane";
    AutinIllinoisNamesMap["'root.syn3A.membrane.outer_membrane'"] = "LipidOuterMembrane";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.RIBOSOME_70S'"] = "particle_RIBOSOME";

    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.GLY_tRNASYNTH'"] = "JCVISYN3A_0405";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.THR_tRNASYNTH'"] = "JCVISYN3A_0222";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.SER_tRNASYNTH'"] = "JCVISYN3A_0061";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.HIS_tRNASYNTH'"] = "JCVISYN3A_0288";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.ALA_tRNASYNTH'"] = "JCVISYN3A_0163";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.TRP_tRNASYNTH'"] = "JCVISYN3A_0308";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.MET_tRNASYNTH'"] = "JCVISYN3A_0012";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.LEU_tRNASYNTH'"] = "JCVISYN3A_0634";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.CYS_tRNASYNTH'"] = "JCVISYN3A_0637";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.GLU_tRNASYNTH'"] = "JCVISYN3A_0126";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.TYR_tRNASYNTH'"] = "JCVISYN3A_0613";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.VAL_tRNASYNTH'"] = "JCVISYN3A_0260";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.PRO_tRNASYNTH'"] = "JCVISYN3A_0282";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.ILE_tRNASYNTH'"] = "JCVISYN3A_0519";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.ARG_tRNASYNTH'"] = "JCVISYN3A_0535";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.PHE_tRNASYNTH'"] = "JCVISYN3A_0529";
    AutinIllinoisNamesMap["'root.syn3A.interior.proteins.ASP_tRNASYNTH'"] = "JCVISYN3A_0287";
}

EntityIdInt GetParticleKindIdFromGeneIdOrName(const string& ParticleKindName, const EntityIdInt ParticleKindIdParam, const std::unordered_map<EntityIdInt, UniqueIdInt>& ProteinIdFromGeneIdTranslator, const unordered_map<string, string>& AutinIllinoisNameMap)
{
    try
    {
        if (auto GeneId = ProteinIdFromGeneIdTranslator.find(ParticleKindIdParam); GeneId != ProteinIdFromGeneIdTranslator.end())
             if (auto ParticleKindId = ParticlesKindsManagerObject.GetParticleKindFromGeneId(GeneId->second); ParticleKindId.has_value() == true)
                 return ParticleKindId->EntityId;

        if (auto IllinoisNameIter = AutinIllinoisNameMap.find(ParticleKindName); IllinoisNameIter != AutinIllinoisNameMap.end())
            return ParticlesKindsManagerObject.GetParticleKindFromStrId(IllinoisNameIter->second)->EntityId;
    }
    CATCH("translating entity id")

    return ParticlesKindsManagerObject.GetParticleKindFromStrId("M_coa_c")->EntityId;
}

void AddParticleKindGraphicDataFromConfigXMLData(const EntityIdInt EntityId)
{
    try
    {
        if (auto ParticleKindObjectIterator = ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.find(EntityId); ParticleKindObjectIterator == ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.end())
        {
            auto OthersParticleKindObjectIterator = ParticlesKindsManagerObject.GraphicParticlesKindsFromConfigXML.find(10000);
            ParticlesKindsManagerObject.ParticlesKinds[EntityId].GraphicData = ParticleKindGraphicData{ EntityId, OthersParticleKindObjectIterator->second.Visible, false, OthersParticleKindObjectIterator->second.SizeX, OthersParticleKindObjectIterator->second.SizeY, OthersParticleKindObjectIterator->second.SizeZ, OthersParticleKindObjectIterator->second.ParticleColor, OthersParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), OthersParticleKindObjectIterator->second.NameFromXML, "NAME_FROM_DATA_FILE" };
        }
        else
            ParticlesKindsManagerObject.ParticlesKinds[EntityId].GraphicData = ParticleKindGraphicData{ EntityId, ParticleKindObjectIterator->second.Visible, false, ParticleKindObjectIterator->second.SizeX, ParticleKindObjectIterator->second.SizeY, ParticleKindObjectIterator->second.SizeZ, ParticleKindObjectIterator->second.ParticleColor, ParticleKindObjectIterator->second.ParticleColor, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectIterator->second.NameFromXML, "NAME_FROM_DATA_FILE"};
    }
    CATCH("adding particle kind graphic data from xml config data")
}

vector3_Real32 CountCenterForListOfAtoms(const vector<CellEngineAtom>& ListOfAtoms)
{
    vmath::vec3 CenterOfParticle(0.0, 0.0, 0.0);

    try
    {
        for (const CellEngineAtom& AtomObject : ListOfAtoms)
            CenterOfParticle += AtomObject.Position();
    }
    CATCH("counting center for list of atoms")

    return { CenterOfParticle.X() / static_cast<RealType>(ListOfAtoms.size()), CenterOfParticle.Y() / static_cast<RealType>(ListOfAtoms.size()), CenterOfParticle.Z() / static_cast<RealType>(ListOfAtoms.size()) };
}

void CellEngineCIFDataFileReader::ReadDataFromCIFFile(const bool SetStartValuesBool)
{
    try
    {
        if (SetStartValuesBool == true)
            SetStartValues();

        ChainsNames.clear();
        TransformationsMatrixes.clear();

        string Line;
        vector<CellEngineAtom> LocalCellEngineAllAtomsObject;
        vector<CellEngineAtom> LocalCellEngineParticlesCentersObject;

        unordered_map<string, string> AutinIllinoisNamesMap;
        AssociateAutinNameWithIllinoisName(AutinIllinoisNamesMap);
        std::unordered_map<EntityIdInt, UniqueIdInt> ProteinIdFromGeneIdTranslator;
        std::unordered_map<EntityIdInt, string> ParticleAutinKindIdToAutinNameTranslator;

        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == true)
        {
            EntityIdInt LocalParticleKindId = max_element(ParticlesKindsManagerObject.ParticlesKinds.begin(), ParticlesKindsManagerObject.ParticlesKinds.end(), [](const pair<EntityIdInt, ParticleKind>& lhs, const pair<EntityIdInt, ParticleKind>& rhs){ return lhs.first < rhs.first; })->first + 1;
            ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::Lipid, LocalParticleKindId, "LipidOuterMembrane", "LipidOuterMembrane", "LipidOuterMembrane", -1, 0, "m", 1);
            ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::Lipid, LocalParticleKindId, "LipidInnerMembrane", "LipidInnerMembrane", "LipidInnerMembrane", -1, 0, "m", 1);
        }

        std::ifstream File(string(CellEngineConfigDataObject.CellStateFileName).c_str(), std::ios_base::in);

        LoggersManagerObject.Log(STREAM("STARTED READING FROM CIF FILE"));

        UnsignedInt NumberOfAtoms = 0;
        UnsignedInt NumberOfAtomsDNA = 0;
        UnsignedInt NumberOfNucleotidesInDNA = 0;
        UnsignedInt NumberOfParticles = 0;

        CellEngineColorsObject.SelectRandomEngineForColors();

        smatch SMatchObject;

        vector<UnsignedInt> AppliedMatrixesIds;
        regex RegexObject1(R"((\d+)-(\d+)|\((\d+)\))");

        vector<string> AppliedChainsNames;
        regex RegexObject2("([A-z]+[0-9]*)");

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == ". . ")
            {
                vector<string> AtomFields = split(Line, " ");

                EntityIdInt EntityId = stoi(AtomFields[2]);
                ParticlesKindsManagerObject.ParticlesKinds[EntityId].EntityId = EntityId;

                if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == true)
                {
                    const string StrToFind = "SYN3A_0";
                    auto Pos = AtomFields[5].find(StrToFind);
                    if (Pos != string::npos)
                        ProteinIdFromGeneIdTranslator[EntityId] = stoi(AtomFields[5].substr(Pos + StrToFind.length(), 3));
                    ParticleAutinKindIdToAutinNameTranslator[EntityId] = AtomFields[5];
                }
            }
            else
            if (Line.substr(0, 4) == "ATOM")
            {
                CellEngineAtom CellEngineAtomObject = ParseRecord(Line.c_str());
                ChainsNames[CellEngineAtomObject.Chain].emplace_back(CellEngineAtomObject);
            }
            else
            if (Line.find("point symmetry operation") != std::string::npos)
            {
                vector<string> MatrixFields = split(Line, " ");

                TransformationMatrix3x4 TransformationMatrix3x4Object{};

                UnsignedInt TransformationMatrix3x4ObjectId = stoi(MatrixFields[0]);

                CopyFieldsToMatrix(TransformationMatrix3x4Object, MatrixFields);

                TransformationsMatrixes[TransformationMatrix3x4ObjectId] = TransformationMatrix3x4Object;
            }
            else
            if (Line.substr(0, 4) == "1 '(")
            {
                AppliedMatrixesIds.clear();
                AppliedMatrixesIds.clear();
                auto pos = Line.cbegin();
                auto end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject1); pos = SMatchObject.suffix().first)
                    if (SMatchObject.str(1).empty() == false && SMatchObject.str(2).empty() == false)
                        for (UnsignedInt MatrixId = stoi(SMatchObject.str(1)); MatrixId <= stoi(SMatchObject.str(2)); MatrixId++)
                            AppliedMatrixesIds.emplace_back(MatrixId);
                    else
                        AppliedMatrixesIds.emplace_back(stoi(SMatchObject.str(3)));

                AppliedChainsNames.clear();
                pos = Line.cbegin();
                end = Line.cend();
                for ( ; regex_search(pos, end, SMatchObject, RegexObject2); pos = SMatchObject.suffix().first)
                    AppliedChainsNames.emplace_back(SMatchObject.str(1));

                vmath::vec3 ChainColor = CellEngineColorsObject.GetRandomColor();

                for (const auto& AppliedMatrixId : AppliedMatrixesIds)
                {
                    LocalCellEngineAllAtomsObject.clear();

                    if (AppliedChainsNames.empty() == true)
                        if (WriteErrorsInCIFFileToScreen == true)
                            LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE AppliedChainsNames EMPTY = " << to_string(AppliedMatrixId)));

                    for (const auto& AppliedChainName : AppliedChainsNames)
                    {
                        auto AtomsForChainNameIterator = ChainsNames.find(AppliedChainName);
                        if (AtomsForChainNameIterator == ChainsNames.end())
                        {
                            if (WriteErrorsInCIFFileToScreen == true)
                                LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE LACKS AppliedChainName: " << AppliedChainName));
                            continue;
                        }
                        else
                        {
                            if (AtomsForChainNameIterator->second.empty() == true)
                                if (WriteErrorsInCIFFileToScreen == true)
                                    LoggersManagerObject.Log(STREAM("ERROR IN CIF FILE - ChainsNames.find(AppliedChainName)->second.size() == 0 " << AppliedChainName));

                            vmath::vec3 UniqueParticleColor = CellEngineColorsObject.GetRandomColor();

                            EntityIdInt LocalEntityId = AtomsForChainNameIterator->second[0].EntityId;

                            if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == true)
                            {
                                LocalEntityId = GetParticleKindIdFromGeneIdOrName(ParticleAutinKindIdToAutinNameTranslator.find(LocalEntityId)->second, LocalEntityId, ProteinIdFromGeneIdTranslator, AutinIllinoisNamesMap);
                                AddParticleKindGraphicDataFromConfigXMLData(LocalEntityId);
                            }

                            if (CellEngineUseful::IsNucleotide(AppliedChainName))
                                NumberOfNucleotidesInDNA++;

                            ListOfAtomsType ListOfAtoms;

                            for (auto AppliedAtom : AtomsForChainNameIterator->second)
                            {
                                NumberOfAtoms++;

                                if (CellEngineUseful::IsNucleotide(AppliedChainName))
                                    NumberOfAtomsDNA++;

                                auto TransformationMatrixIterator = TransformationsMatrixes.find(AppliedMatrixId);

                                auto Result = CountResultPositionsFromTransformationMatrix(TransformationMatrixIterator, AppliedAtom);

                                AppliedAtom.SetAtomPositionsData(Result[0], Result[1], Result[2]);

                                ListOfAtoms.emplace_back(AppliedAtom);
                            }

                            vector3_Real32 Center = CountCenterForListOfAtoms(ListOfAtoms);

                            if (ListOfAtoms.empty() == false)
                            {
                                NumberOfParticles++;

                                SetCurrentSectorPos(CellEngineUseful::GetSectorPos(Center.X, Center.Y, Center.Z));

                                Particle ParticleToInsert(NumberOfParticles, LocalEntityId, CellEngineUseful::GetChainIdFromChainName(AppliedChainName), -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(UniqueParticleColor));
                                ParticleToInsert.ListOfAtoms = ListOfAtoms;
                                ParticleToInsert.Center = Center;
                                ParticleToInsert.ParticleColor = ParticlesKindsManagerObject.GetParticleKind(ParticleToInsert.EntityId).GraphicData.ParticleColor;
                                ParticleToInsert.UniqueParticleColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(UniqueParticleColor);
                                ParticleToInsert.RandomParticleKindColor = CellEngineUseful::GetVector3FormVMathVec3ForColor(ChainColor);
                                AddNewParticle(ParticleToInsert);
                            }
                        }
                    }
                }
            }
        }

        LoggersManagerObject.Log(STREAM("Number of erased particles from particle kinds = " << erase_if(ParticlesKindsManagerObject.ParticlesKinds, [](auto& P){ return P.first < CellEngineConfigDataObject.DNAIdentifier; })));

        PreprocessData(true);

        PrintStatistics();

        LoggersManagerObject.Log(STREAM("NumberOfAtoms = " << NumberOfAtoms << " NumberOfParticles = " << NumberOfParticles << " | LocalCellEngineParticlesCentersObject.size() = " << LocalCellEngineParticlesCentersObject.size() << " | NumberOfAtomsDNA = " << NumberOfAtomsDNA << " | NumberOfNucleotidesInDNA = " << NumberOfNucleotidesInDNA <<" | AtomsPositionsMatrixes.size() = " << TransformationsMatrixes.size() << " | " << endl));

        LoggersManagerObject.Log(STREAM("FINISHED READING FROM CIF FILE"));

        File.close();
    }
    CATCH("reading data from CIF file")
}