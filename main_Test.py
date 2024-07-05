import os
import platform
import numpy as np
from tensorflow.keras import models
import tensorflow as tf
from time import time

def FDR_Filter(ytest, ytest_predict, fdr = 0.01):
    ytest_predict2 = []
    ytest_predict_prob = []
    for yt in ytest_predict:
        ytest_predict2.append(yt.index(max(yt)))
        ytest_predict_prob.append(max(yt))

    sort_ytest_predict_prob = sorted(range(len(ytest_predict_prob)), key=lambda k: ytest_predict_prob[k], reverse=True)
    maxin = sort_ytest_predict_prob[0: int(len(sort_ytest_predict_prob) * (1 - fdr))]

    # for s in range(len(maxin)):
    print(ytest_predict_prob)
    print(sort_ytest_predict_prob)
    print(maxin)
    print(ytest)



if __name__ == '__main__':
    ModelName = ["ResNet101"]
    FileName = ['PXD041292_97H_Beishida',
                'PXD041292_97H_dalian',
                'PXD041292_97H_Huagong',
                'PXD030166_U20S',
                'PXD014124_HeLa',
                'PXD036650_K562',
                'PXD034913_skeletal_muscle',
                'PXD035780_saliva',
                'PXD036800_Urine',
                'PXD005371_PXD005372_HCT116',
                'PXD030166_HEK293T',
                'PXD022868']

    system = platform.system()
    print(system)
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    Thre = 0.7

    if 'Linux' in str(system):
        FilePath = "/media/share/MSData_DDA/Human/"
        OutputPath = "/media/share/NPYFile/"
        Connector = '/'
        ModelPath = "/media/share/Train_Model/"

    else:
        FilePath = "Z:\\MSData_DDA\\Human\\"
        OutputPath = "Z:\\NPYFile\\"
        Connector = '\\'
        ModelPath = "Z:\\Train_Model\\"

    Para = {"file_path": FilePath,
            "Species": "Human",
            "Output_Path": OutputPath,
            "IP_str": "10.0.0.175",
            "Database": "MassSpec_Human_DDA",
            "Connector": Connector,
            "ModelPath": ModelPath,
            "leng": 512,
            "Channel": 1}

    for m in ModelName:
        model = models.load_model(Para["ModelPath"] + "HumanDDA_" + m + ".model")
        print(m)
        for f in FileName:
            print(f)
            TestData = np.load(Para["Output_Path"] + f + "_SpectraTestData2.npy", allow_pickle=True).item()
            xtest = TestData["xTrain"]
            ytest = TestData["yTrain"]
            Sequence = TestData["SeqTrain"]
            # PepMass = TestData["PepMass"]

            print(len(xtest), len(xtest[0]))

            testX_tfcon = tf.constant(xtest, shape=[len(xtest), Para["leng"], Para["Channel"]])
            print(testX_tfcon.shape)

            begin_time = time()
            testY_tfcon = model.predict(testX_tfcon)
            end_time = time()
            print("Test Time = ", end_time - begin_time, "s")

            ytest_predict = testY_tfcon.tolist()
            # print(ytest_predict[0])
            # ytest_predict2 = []
            # for yt in ytest_predict:
            #     ypre = []
            #     for y in yt:
            #         if y >= Thre:
            #             ypre.append(yt.index(y))
            #     ytest_predict2.append(ypre)
            #
            # sss = 0
            # for i in range(len(ytest)):
            #     if ytest[i] in ytest_predict2[i]:
            #         sss += 1
            #
            # print(sss/testX_tfcon.shape[0])

            FDR_Filter(ytest, ytest_predict, fdr=0.01)



            print("===============================\n")
        print("**********************************\n\n")
