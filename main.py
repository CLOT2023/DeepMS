import os
import platform
import ReadPostgreSQL
import numpy as np

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    system = platform.system()
    print(system)

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    if 'Linux' in str(system):
        FilePath = "/media/share/MSData_DDA/Human/"
        OutputPath = "/media/share/NPYFile2/"
        Connector = '/'
        ModelPath = "/media/share/Train_Model/"

    else:
        FilePath = "Z:\\MSData_DDA\\Human\\"
        OutputPath = "Z:\\NPYFile2\\"
        Connector = '\\'
        ModelPath = "Z:\\Train_Model\\"

    # modelNN = ["AlexNet", "VGG16", "VGG19", "GoogleNet", "ResNet50", "ResNet101", "ResNet152", "DenseNet121",
    #            "DenseNet169", "DenseNet201", "DenseNet264"]
    modelNN = ["VGG16", "GoogleNet", "ResNet50", "DenseNet121"]

    Para = {"file_path": FilePath,
            "Species": "Human",
            "Output_Path": OutputPath,
            "IP_str": "10.0.0.175",
            "Database": "MassSpec_Human_DDA",
            "Connector": Connector,
            "ModelPath": ModelPath,
            "leng": 512,
            "Channel": 1,
            "Model": modelNN}

    # TrainData = ReadPostgreSQL.Data_Processing(Para, 4)
    TrainData = np.load(Para["file_path"] + "SpectraTrainData.npy", allow_pickle=True).item()
    trainX_tfcon, trainY_tfcon, lenMS, Sequence, value, strategy = ReadPostgreSQL.Train_Processing(Para, TrainData)
    ReadPostgreSQL.Train_Model(Para, trainX_tfcon, trainY_tfcon, lenMS, Sequence, value, strategy)
