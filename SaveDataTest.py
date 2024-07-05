import os
from time import time
import Read_MGFTest
import platform
import numpy as np

if __name__ == '__main__':
    system = platform.system()
    print(system)

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    TestFileName = ['PXD041292_97H_Beishida', 'PXD041292_97H_dalian', 'PXD041292_97H_Huagong', 'PXD014124_HeLa', 'PXD030166_U20S',
                    'PXD036650_K562', 'PXD034913_skeletal_muscle', 'PXD035780_saliva', 'PXD036800_Urine', 'PXD005371_PXD005372_HCT116',
                    'PXD030166_HEK293T', 'PXD022868']
    # TestFileName = ['PXD022868']
    if 'Linux' in str(system):
        FilePath = "/media/share/MSData_DDA/Human/"
        TestFilePath = "/media/share/MSData_DDA/Test/"
        OutputPath = "/media/share/NPYFile/"
        Connector = '/'
        ConnectorFlag = 5

    else:
        FilePath = "Z:\\MSData_DDA\\Human\\"
        TestFilePath = "Z:\\MSData_DDA\\Test\\"
        OutputPath = "Z:\\NPYFile\\"
        Connector = '\\'
        ConnectorFlag = 3

    Para = {"file_path" : FilePath,
            "Species" : "Homo_Sapiens",
            "Test_file_path" : TestFilePath,
            "Output_Path" : OutputPath,
            "IP_str" : "10.0.0.175",
            "Database" : "MassSpec_Human_DDA",
            "Connector" : Connector,
            "ConnectorFlag" : ConnectorFlag}

    SeqClass = np.load(Para["file_path"] + "Sequence_Class2.npy", allow_pickle=True).item()
    print(len(SeqClass))
    for f in TestFileName:
        print(f)
        begin_time = time()
        # Read_MGFTest.RawResult_ProcessingTest(FileName=f, Parameter=Para)
        Read_MGFTest.File_ProcessingTest(FileName=f, Parameter=Para, SeqClass=SeqClass)
        end_time = time()
        print("\n=======================\n\n")
        print("Train Time = ", end_time - begin_time)

    # Read_MGFTest_Phospho.File_ProcessingTest(Parameter=Para, SeqClass=SeqClass)