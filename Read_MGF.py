import copy
import os
import random
import numpy as np
# import psycopg2
import math

def traverse_second_level_directories(root_dir):
    Second_level_dir = []
    for dir_name in os.listdir(root_dir):
        sub_dir_path = os.path.join(root_dir, dir_name)
        if os.path.isdir(sub_dir_path):
            for second_level_dir in os.listdir(sub_dir_path):
                second_level_dir_path = os.path.join(sub_dir_path, second_level_dir)
                if os.path.isdir(second_level_dir_path):
                    Second_level_dir.append(second_level_dir_path)

    return Second_level_dir

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
                # print(rootname)

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


def Seq_Sat(filepath, Connector):
    print("Seq_Sat")
    Resultdata = {}
    ResultProtein = {}
    ResultProtein2 = {}
    Spectra_seq = {}
    Spectra_seq2 = {}
    for root, dirs, files in os.walk(filepath):
        for file in files:
            if "Ions-Filtered.npy" in file:
                IonsFilter = np.load(root + Connector + file, allow_pickle=True).item()
                print(root + Connector + file, len(IonsFilter))

                for key in IonsFilter:
                    if "Phospho" not in str(IonsFilter[key]["Modification"]):
                        Seq_Class = IonsFilter[key]["Sequence"]
                        if Seq_Class not in Resultdata:
                            Resultdata[Seq_Class] = []
                        Resultdata[Seq_Class].append(key)
                        ResultProtein[Seq_Class] = IonsFilter[key]["Proteins"]

    seqindex = 0
    for key in Resultdata:
        if len(Resultdata[key]) >= 25000 and len(Resultdata[key]) <= 40000:
            Spectra_seq2[key] = seqindex
            ResultProtein2[key] = ResultProtein[key]
            seqindex += 1
            for p in Resultdata[key]:
                Spectra_seq[p] = key

    np.save(filepath + Connector + "Sequence_Class.npy", Spectra_seq2)
    np.save(filepath + Connector + "Spectra_Sequence-Filtered.npy", Spectra_seq)
    np.save(filepath + Connector + "Sequence_Protein.npy", ResultProtein2)
    print(len(Spectra_seq2), len(Spectra_seq), len(ResultProtein2))
    return Spectra_seq2, Spectra_seq


def ReadMGFFile(filepath, Spectra_seq, SeqClass, Connector, leng):
    MSData = {}
    MSData["xTrain"] = []
    MSData["yTrain"] = []
    MSData["MassTrain"] = []
    MSData["SeqTrain"] = []
    MSData["LenClass"] = len(SeqClass)

    for root, dirs, files in os.walk(filepath):
        for file in files:
            if "MGFDic.npy" in file:
                MGFDic = np.load(root + Connector + file, allow_pickle=True).item()
                print(os.path.join(root, file), len(MGFDic))

                for key in MGFDic:
                    if key in Spectra_seq:
                        MSData["MassTrain"].append(MGFDic[key]["MassTrain"])

                        if len(MGFDic[key]["MZ"]) > leng:
                            MSData["xTrain"].append(Fixed_Length(MGFDic[key]["MZ"], MGFDic[key]["Intensity"], leng))
                        elif len(MGFDic[key]["MZ"]) < leng:
                            MSData["xTrain"].append(Mz_Fix(MGFDic[key]["MZ"], leng))
                        elif len(MGFDic[key]["MZ"]) == leng:
                            MSData["xTrain"].append(MGFDic[key]["MZ"])

                        MSData["SeqTrain"].append(Spectra_seq[key])
                        lst1 = [0 for _ in range(len(SeqClass))]
                        lst1[SeqClass[Spectra_seq[key]]] = 1
                        MSData["yTrain"].append(lst1)

    print(len(MSData["SeqTrain"]))
    np.savez_compressed(filepath + Connector + "SpectraTrainData.npz", **MSData)
    np.save(filepath + Connector + "SpectraTrainData.npy", MSData)


def ProcessingMGFFile(filepath, Connector, ConnectorFlag):
    MSData = {}
    Ions_Filtered = {}
    PXDNum = str(filepath).split(Connector)[ConnectorFlag]
    rootname = PXDNum + Connector
    print(PXDNum)

    for root, dirs, files in os.walk(filepath):
        for file in files:
            if "Ions-Filtered.npy" in file:
                Ions = np.load(root + Connector + file, allow_pickle=True).item()
                for key in Ions:
                    Ions_Filtered[key] = Ions[key]

    for root, dirs, files in os.walk(filepath):
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

                for i in range(len(bi)):
                    title = lines[bi[i] + 1].split('=')[1].strip()
                    if rootname + title in Ions_Filtered:
                        MSData[rootname + title] = {}
                        mz = []
                        intensity = []
                        pepmass = float(lines[bi[i] + 4].strip().split('=')[1])
                        MSData[rootname + title]["MassTrain"] = pepmass
                        mzl = lines[bi[i] + 5 : ei[i]]

                        for m in mzl:
                            mzinten = m.strip().split(' ')
                            mz.append(float(mzinten[0]))
                            intensity.append(float(mzinten[1]))

                        MSData[rootname + title]["MZ"] = mz
                        MSData[rootname + title]["Intensity"] = intensity

                f.close()

    print(len(MSData))
    np.save(filepath + Connector + PXDNum + "_MGFDic.npy", MSData)


# 构造CNN数据格式
def Postgre_Data(MSData, conn):
    cursor = conn.cursor()
    for key in MSData:
        # Spectra_MZD
        sqlstr = "INSERT INTO public.\"Spectra_MZD\"(\"Title\", \"MZ_Denoising\", \"Modification\", \"MZLen\", \"Seq-Charge\") VALUES ("
        sqlstr = sqlstr + '\'' + str(key) + '\','
        sqlstr = sqlstr + 'ARRAY' + str(MSData[key]["MZ_Denoising"]) + ','
        sqlstr = sqlstr + 'ARRAY' + str(MSData[key]["Modification"]) + ','
        sqlstr = sqlstr + str(MSData[key]["MZLen"]) + ','
        sqlstr = sqlstr + '\'' + str(MSData[key]["Sequence"]) + " - " + str(MSData[key]["Charge"]) + '\');'

        cursor.execute(sqlstr)
        conn.commit()

        # Spectra_RawData
        sqlstr1 = "INSERT INTO public.\"Spectra_RawData\"(\"Title\", \"FileName\", \"ScanNum\", \"MZ\", \"Intensity\", \"Raw_Score\", \"RtinSeconds\", \"MZLen\") VALUES ("
        sqlstr1 = sqlstr1 + '\'' + str(key) + '\','
        sqlstr1 = sqlstr1 + '\'' + str(MSData[key]["FileName"]) + '\','
        sqlstr1 = sqlstr1 + str(MSData[key]["ScanNum"]) + ','
        sqlstr1 = sqlstr1 + 'ARRAY' + str(MSData[key]["MZ"]) + ','
        sqlstr1 = sqlstr1 + 'ARRAY' + str(MSData[key]["Intensity"]) + ','
        sqlstr1 = sqlstr1 + str(MSData[key]["Raw_Score"]) + ','
        sqlstr1 = sqlstr1 + str(MSData[key]["RtinSeconds"]) + ','
        sqlstr1 = sqlstr1 + str(MSData[key]["MZLen"]) + ');'

        cursor.execute(sqlstr1)
        conn.commit()

        # Spectra_Main
        sqlstr2 = "INSERT INTO public.\"Spectra_Main\"(\"Title\", \"Charge\", \"Sequence\", \"Seq_Charge\", \"PTM\", \"Proteins\", \"PepMass\", \"MZLen\") VALUES ("
        sqlstr2 = sqlstr2 + '\'' + str(key) + '\','
        sqlstr2 = sqlstr2 + str(MSData[key]["Charge"]) + ','
        sqlstr2 = sqlstr2 + '\'' + str(MSData[key]["Sequence"]) + '\','
        sqlstr2 = sqlstr2 + '\'' + str(MSData[key]["Sequence"]) + " - " + str(MSData[key]["Charge"]) + '\','
        sqlstr2 = sqlstr2 + 'ARRAY' + str(MSData[key]["PTM"]) + ','
        sqlstr2 = sqlstr2 + 'ARRAY' + str(MSData[key]["Proteins"]) + ','
        sqlstr2 = sqlstr2 + str(MSData[key]["PepMass"]) + ','
        sqlstr2 = sqlstr2 + str(MSData[key]["MZLen"]) + ');'

        cursor.execute(sqlstr2)
        conn.commit()

    cursor.close()

def File_Processing(Parameter):
    print(Parameter["file_path"])
    # Spectra_seq = Seq_Sat(Parameter["file_path"], Parameter["Connector"])
    # print("\n========================\n")

    Spectra_seq = np.load(Parameter["file_path"] + "Spectra_Sequence-Filtered.npy", allow_pickle=True).item()
    SeqClass = np.load(Parameter["file_path"] + "Spectra_Sequence.npy", allow_pickle=True).item()

    # ReadMGFFile(Parameter["file_path"], Spectra_seq, SeqClass, Parameter["Connector"], Parameter["ConnectorFlag"])

def File_Processing2(Parameter):
    print(Parameter["file_path"])
    # SeqClass, Spectra_seq = Seq_Sat(Parameter["file_path"], Parameter["Connector"])

    SeqClass = np.load(Parameter["file_path"] + "Sequence_Class.npy", allow_pickle=True).item()
    Spectra_seq = np.load(Parameter["file_path"] + "Spectra_Sequence-Filtered.npy", allow_pickle=True).item()
    print(len(SeqClass), len(Spectra_seq))

    ReadMGFFile(Parameter["file_path"], Spectra_seq, SeqClass, Parameter["Connector"],512)


def RawResult_Processing(FileName, Parameter):
    print(Parameter["file_path"] + FileName)
    Spectra_Ions(Parameter["file_path"] + FileName, Parameter["Connector"], Parameter["ConnectorFlag"])
    print("\n========================\n")




