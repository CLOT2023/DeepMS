import os
from time import time
import Read_MGF
import platform
import numpy as np
import tensorflow as tf


if __name__ == '__main__':
    system = platform.system()
    print(system)

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    FileName = ['HomoSapiens_10', 'HomoSapiens_20', 'HomoSapiens_30']
    if 'Linux' in str(system):
        FilePath = "/media/share/MSData_DDA/Human/"
        TestFilePath = "/media/share/MSData_DDA/Test/"
        OutputPath = "/media/share/NPYFile/"
        Connector = '/'
        ConnectorFlag = 7

    else:
        FilePath = "Z:\\MSData_DDA\\Human\\"
        TestFilePath = "Z:\\MSData_DDA\\Test\\"
        OutputPath = "Z:\\NPYFile\\"
        Connector = '\\'
        ConnectorFlag = 5

    Para = {"file_path" : FilePath,
            "Species" : "Homo_Sapiens",
            "Test_file_path" : TestFilePath,
            "Output_Path" : OutputPath,
            "IP_str" : "10.0.0.175",
            "Database" : "MassSpec_Human_DDA",
            "Connector" : Connector,
            "ConnectorFlag" : ConnectorFlag}

    #
    # for f in FileName:

    #     print(f)
    #     Second_level_dir = Read_MGF.traverse_second_level_directories(Para["file_path"] + f)
    #     print(Second_level_dir)
    #     begin_time = time()
    #     Read_MGF.RawResult_Processing(FileName=f, Parameter=Para)
    #     for s in Second_level_dir:
    #         Read_MGF.ProcessingMGFFile(s, Para["Connector"], Para["ConnectorFlag"])
    #     end_time = time()
    #     print("\n=======================\n\n")
    #     print("Train Time = ", end_time - begin_time)

    Read_MGF.File_Processing2(Parameter=Para)
    # SeqClass = np.load(Para["file_path"] + "Sequence_Class.npy", allow_pickle=True).item()
    # Spectra_seq = np.load(Para["file_path"] + "Spectra_Sequence-Filtered.npy", allow_pickle=True).item()
    #
    # SeqClass2 ={}
    # SeqClass3 = {}
    # Spectra_seq2 = {}
    # for key in Spectra_seq:
    #     if Spectra_seq[key] not in SeqClass2:
    #         SeqClass2[Spectra_seq[key]] = []
    #
    #     SeqClass2[Spectra_seq[key]].append(key)
    #
    # index = 0
    # for key in SeqClass2:
    #     if len(SeqClass2[key]) >= 40000:
    #         SeqClass3[key] = index
    #         index += 1
    #
    #
    # for key in Spectra_seq:
    #     if Spectra_seq[key] in SeqClass3:
    #         Spectra_seq2[key] = Spectra_seq[key]
    #
    #
    # print(len(SeqClass3), len(Spectra_seq2))
    #
    # MSSpectra = np.load(Para["file_path"] + "SpectraTrainData.npy", allow_pickle=True).item()
    # MSData = {}
    # MSData["xTrain"] = []
    # MSData["yTrain"] = []
    # MSData["MassTrain"] = []
    # MSData["SeqTrain"] = []
    # MSData["LenClass"] = len(SeqClass3)
    #
    # for i in range(len(MSSpectra["SeqTrain"])):
    #     if MSSpectra["SeqTrain"][i] in SeqClass3:
    #         MSData["xTrain"].append(MSSpectra["xTrain"][i])
    #         MSData["MassTrain"].append(MSSpectra["MassTrain"][i])
    #         MSData["SeqTrain"].append(MSSpectra["SeqTrain"][i])
    #
    #         lst1 = [0 for _ in range(len(SeqClass3))]
    #         lst1[SeqClass3[MSSpectra["SeqTrain"][i]]] = 1
    #         MSData["yTrain"].append(lst1)
    #
    #
    # print(len(MSData["xTrain"]))
    # np.save(Para["file_path"] + "Sequence_Class3.npy", SeqClass3)
    # np.save(Para["file_path"] + "Spectra_Sequence-Filtered3.npy", Spectra_seq2)
    # np.save(Para["file_path"] + "SpectraTrainData3.npy", MSData)



