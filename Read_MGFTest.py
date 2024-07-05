import copy
import os
import random
import numpy as np
# import psycopg2
import math


# 保留最高峰度的mz值
def Fixed_Length(MSmz, MSRelaIntensity, length):
    mzlist = copy.deepcopy(MSmz)
    Relaintensitylist = copy.deepcopy(MSRelaIntensity)
    mz = []
    sortintensity = sorted(range(len(Relaintensitylist)), key=lambda k: Relaintensitylist[k], reverse=True)

    maxin = sortintensity[0 : length]

    for i in range(len(mzlist)):
        if i in maxin:
            mz.append(mzlist[i])

    if len(mz) != length:
        print("***********", len(mz))

    return mz

def Mz_Fix(mz, length):
    mzd = []
    mdl = length - len(mz)
    if mdl % 2 == 0:
        mzlen = mdl // 2
        for i in range(mzlen):
            mzd.append(0)
        for m in mz:
            mzd.append(m)
        for i in range(mzlen):
            mzd.append(0)

    else:
        mzlen = (mdl - 1) // 2
        for i in range(mzlen):
            mzd.append(0)
        mzd.append(0)
        for m in mz:
            mzd.append(m)
        for i in range(mzlen):
            mzd.append(0)

    return mzd



def Spectra_Ions(filepath, Connector, ConnectorFlag):
    print("Spectra_Ions")
    Homo_Sapiens = np.load("homo_sapiens.npy", allow_pickle=True).item()
    for root, dirs, files in os.walk(filepath):
        for file in files:
            if "pFind-Filtered.spectra" in file:
                Resultdata = {}
                Ionsdata = {}
                datalist = []
                ListName = {}
                f = open(root + Connector + file)
                print(root + Connector + file)
                lines = f.readlines()
                Name = lines[0].split('	')
                # print(len(Name))
                for i in range(len(Name)):
                    if 'File_Name' in Name[i]:
                        ListName["Title"] = i
                    elif 'Sequence' in Name[i]:
                        ListName["Sequence"] = i
                    elif 'Modification' in Name[i]:
                        ListName["Modification"] = i
                    elif 'Proteins' in Name[i]:
                        ListName["Proteins"] = i
                    elif 'Scan_No' in Name[i]:
                        ListName["ScanNum"] = i
                    elif 'Charge' in Name[i]:
                        ListName["Charge"] = i
                    elif 'Raw_Score' in Name[i]:
                        ListName["Raw_Score"] = i

                del lines[0]
                # print(len(lines))
                for i in range(len(lines)):
                    line = lines[i].split('	')
                    if len(line[ListName["Sequence"]].strip()) >= 9:
                        datalist.append(line)

                PXDNum = str(root).split(Connector)[ConnectorFlag]
                rootname = PXDNum + Connector
                print(rootname)

                for d in datalist:
                    protein = []
                    prot = d[ListName["Proteins"]].strip().split('/')
                    # print(prot)
                    for p in prot:
                        if p.strip() in Homo_Sapiens:
                            protein.append(p.strip())

                    if len(protein) > 0:
                        Resultdata[rootname + d[ListName["Title"]].strip()] = {}
                        Resultdata[rootname + d[ListName["Title"]].strip()]["FileName"] = str(d[ListName["Title"]].strip().split('.')[0])
                        Resultdata[rootname + d[ListName["Title"]].strip()]["ScanNum"] = int(d[ListName["ScanNum"]].strip())
                        Resultdata[rootname + d[ListName["Title"]].strip()]["Charge"] = int(d[ListName["Charge"]].strip())
                        Resultdata[rootname + d[ListName["Title"]].strip()]["Sequence"] = str(d[ListName["Sequence"]].strip())
                        Resultdata[rootname + d[ListName["Title"]].strip()]["Raw_Score"] = float(d[ListName["Raw_Score"]].strip())
                        Resultdata[rootname + d[ListName["Title"]].strip()]["Proteins"] = protein

                        modification = []
                        modi = d[ListName["Modification"]].strip().split(';')
                        for m in modi:
                            modification.append(m)

                        Resultdata[rootname + d[ListName["Title"]].strip()]["Modification"] = modification
                        Resultdata[rootname + d[ListName["Title"]].strip()]["USI"] = "mzspec:" + PXDNum + ":" + str(d[ListName["Title"]].strip().split('.')[0]) + ":scan:" + d[ListName["ScanNum"]].strip() + ":" + d[ListName["Sequence"]].strip() + "/" + d[ListName["Charge"]].strip() + ";"

                print(len(Resultdata))
                f = open(root[:-6] + "Ions")
                print(root[:-6] + "Ions")
                lines = f.readlines()
                index = 0
                bi = []
                ei = []
                while (index < len(lines)):
                    if "#" in lines[index] and ".dta" in lines[index]:
                        bi.append(index)
                    index += 1
                for i in range(len(bi) - 1):
                    ei.append(bi[i + 1] - 2)
                ei.append(len(lines) - 2)

                # 处理每一张谱图的峰匹配结果
                for i in range(len(bi)):
                    title = lines[bi[i]].split('#')[1].strip()
                    if rootname + title in Resultdata:
                        bymatch = {}
                        mzie = []
                        mzl = lines[bi[i] + 2: ei[i]]
                        # print(mzl)
                        for m in mzl:
                            mli = []
                            mzinten = m.strip().split("	")
                            if mzinten[2] != '0':
                                bymatch[mzinten[0].rstrip('+').strip()] = 1

                                mli.append(mzinten[0])
                                mli.append(float(mzinten[1]))
                                mli.append(float(mzinten[2]))
                                mli.append(float(mzinten[3]))
                                mzie.append(mli)

                        Ionsdata[rootname + title] = Resultdata[rootname + title]
                        Ionsdata[rootname + title]["byMatch"] = sorted(bymatch)
                        Ionsdata[rootname + title]["PeakMatch"] = mzie

                print(len(Ionsdata))
                np.save(root[:-6] + "Ions.npy", Ionsdata)

                IonsFilter = {}
                for key in Ionsdata:
                    if len(Ionsdata[key]["byMatch"]) >= math.ceil(len(Ionsdata[key]["Sequence"]) / 2):
                        IonsFilter[key] = Ionsdata[key]

                print(len(IonsFilter))
                np.save(root[:-6] + "Ions-Filtered.npy", IonsFilter)


def Seq_SatTest(filepath, FileName, Connector, SeqClass):
    print("Seq_Sat")
    Spectra_seq = {}
    for root, dirs, files in os.walk(filepath + FileName):
        for file in files:
            if "Ions-Filtered.npy" in file:
                IonsFilter = np.load(root + Connector + file, allow_pickle=True).item()
                print(root + Connector + file)

                for key in IonsFilter:
                    if "Phospho" not in str(IonsFilter[key]["Modification"]):
                        Seq_Class_Test = IonsFilter[key]["Sequence"]
                        if Seq_Class_Test in SeqClass:
                            Spectra_seq[key] = Seq_Class_Test

    np.save(filepath + FileName + Connector + FileName + Connector + FileName + "_Spectra_Sequence.npy", Spectra_seq)
    return Spectra_seq


def ReadMGFFileTest(filepath, FileName, SeqClass, Spectra_seq, Connector, ConnectorFlag, leng):
    MSData = {}
    MSData["xTrain"] = []
    MSData["yTrain"] = []
    MSData["MassTrain"] = []
    MSData["SeqTrain"] = []
    MSData["Title"] = []

    for root, dirs, files in os.walk(filepath + FileName):
        for file in files:
            if ".mgf" in file:
                print(os.path.join(root, file))
                f = open(root + Connector + file)
                lines = f.readlines()
                index = 0
                bi = []
                ei = []
                while (index < len(lines)):
                    if "BEGIN IONS" in lines[index]:
                        bi.append(index)
                    elif "END IONS" in lines[index]:
                        ei.append(index)
                    index += 1

                # print("=====>",len(bi), len(ei))
                rootname = str(root).split(Connector)[ConnectorFlag] + Connector
                for i in range(len(bi)):
                    title = lines[bi[i] + 1].split('=')[1].strip()
                    # print(title)
                    if rootname + title in Spectra_seq:
                        mz = []
                        intensity = []
                        pepmass = float(lines[bi[i] + 4].strip().split('=')[1])
                        MSData["MassTrain"].append(pepmass)
                        mzl = lines[bi[i] + 5 : ei[i]]

                        for m in mzl:
                            mzinten = m.strip().split(' ')
                            mz.append(float(mzinten[0]))
                            intensity.append(float(mzinten[1]))
                        if len(mz) > leng:
                            MSData["xTrain"].append(Fixed_Length(mz, intensity, leng))
                        elif len(mz) < leng:
                            MSData["xTrain"].append(Mz_Fix(mz, leng))
                        elif len(mz) == leng:
                            MSData["xTrain"].append(mz)

                        MSData["SeqTrain"].append(Spectra_seq[rootname + title])
                        MSData["yTrain"].append(SeqClass[Spectra_seq[rootname + title]])
                        MSData["Title"].append(rootname + title)

                f.close()
    print(len(MSData["SeqTrain"]))
    np.save(filepath + FileName + Connector + FileName + Connector + FileName + "_SpectraTestData2.npy", MSData)


def File_ProcessingTest(FileName, Parameter, SeqClass):
    print(Parameter["Test_file_path"])
    Spectra_seq = Seq_SatTest(Parameter["Test_file_path"], FileName, Parameter["Connector"], SeqClass)
    print("\n========================\n")

    # Spectra_seq = np.load(Parameter["Test_file_path"] + FileName + Parameter["Connector"] + FileName + Parameter["Connector"] + FileName + "_Spectra_Sequence.npy", allow_pickle=True).item()
    ReadMGFFileTest(Parameter["Test_file_path"], FileName, SeqClass, Spectra_seq, Parameter["Connector"], Parameter["ConnectorFlag"], 512)


def RawResult_ProcessingTest(FileName, Parameter):
    print(Parameter["Test_file_path"] + FileName)
    Spectra_Ions(Parameter["Test_file_path"] + FileName, Parameter["Connector"], Parameter["ConnectorFlag"])
    print("\n========================\n")


