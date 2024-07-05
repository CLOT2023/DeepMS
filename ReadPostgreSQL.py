import copy
import os
import random
# import psycopg2
import numpy as np
import tensorflow as tf
import DDA_DLAlgorithm as ddaalg
# import MLAlgorithms as mla

'''
def Data_Processing(Parameter, CharR):
    conn = psycopg2.connect(
        host=Parameter["IP_str"],  # IP地址
        port="5432",  # 端口号
        database=Parameter["Database"],  # 数据库名
        user="postgres",  # 用户名
        password="123456"  # 密码
    )

    cursor = conn.cursor()
    cursor.execute("SELECT \"PepClass\", \"PepMass\" FROM public.\"Peptide_Class\" where \"Num\">=25000 and \"Num\"<=30000;")
    conn.commit()
    result = cursor.fetchall()
    print(len(result), type(result[0]), result[0][0])

    PepClass = {}
    classindex = 0
    for r in result:
        if r[0] not in PepClass:
            PepClass[r[0]] = {}
        PepClass[r[0]]["PepMass"] = r[1]
        PepClass[r[0]]["Class"] = classindex

        classindex += 1

    np.save(Parameter["Output_Path"] + Parameter["Connector"] + "Human.npy", PepClass)
    print("Save PepClass Successful!\n")

    TrainData = {}
    TrainData["LenMS"] = len(PepClass)
    TrainData["xtrain"] = []
    TrainData["ytrain"] = []
    TrainData["Sequence"] = []
    XT = []
    YT = []
    SQ = []
    inputSQL = "SELECT \"MZ_Denoising\", \"PTM\", \"Charge\", \"Sequence\" FROM public.\"Spectra_MZD\" where "
    for key in PepClass:
        # inputSQL = inputSQL + '\"Sequence\"=\'' + key + "\' or "
        k = key.split(" - ")
        inputSQL = inputSQL + '(\"Sequence\"=\'' + str(k[0]) + "\' and \"Charge\"=" + str(k[1]) + ") or "

    inputSQL = inputSQL[:-4]
    inputSQL += ';'
    print(inputSQL)
    cursor.execute(inputSQL)
    conn.commit()
    result = cursor.fetchall()
    print(len(result), type(result[0]), result[0][3])

    for r in result:
        if 'Phospho' not in list(r[1]):
            try:
                array = np.array(r[0], dtype=int)
                XT.append(r[0])
                lst = [0] * len(PepClass)
                lst[PepClass[str(r[3]) + " - " + str(r[2])]["Class"]] = 1
                # lst[PepClass[str(r[3])]["Class"]] = 1
                YT.append(lst)
                SQ.append(str(r[3]))

            except Exception as e:
                print("Error MZD : ", r[0])
                continue

    DataIndex = [i for i in range(len(XT))]
    for j in range(5):
        random.shuffle(DataIndex)

    for d in DataIndex:
        TrainData["xtrain"].append(XT[d])
        TrainData["ytrain"].append(YT[d])
        TrainData["Sequence"].append(SQ[d])

    print(len(TrainData), len(TrainData["xtrain"]))
    np.save(Parameter["Output_Path"] + Parameter["Connector"] + "Human_DDATrain.npy", TrainData)
    print("Save TrainData Successful!")

    return TrainData
'''

def Train_Processing(Parameter, TrainData):
    xtrain = TrainData["xTrain"]
    ytrain = TrainData["yTrain"]
    lenMS = TrainData["LenClass"]
    Sequence = TrainData["SeqTrain"]

    device = ["/gpu:0", "/gpu:1", "/gpu:2", "/gpu:3"]
    strategy = tf.distribute.MirroredStrategy(devices=device)
    print('Number of device : {}'.format(strategy.num_replicas_in_sync))
    trainX_tfcon = tf.constant(xtrain, shape=[len(xtrain), Parameter["leng"], Parameter["Channel"]])
    trainY_tfcon = tf.constant(ytrain)
    print(trainX_tfcon.shape)
    print(trainY_tfcon.shape)
    print(lenMS)
    value = int(len(xtrain) * 0.8)

    return trainX_tfcon, trainY_tfcon, lenMS, Sequence, value, strategy


def Train_Model(Parameter, trainX_tfcon, trainY_tfcon, lenMS, Sequence, value, strategy):
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

    for m in Parameter["Model"]:
        model = ddaalg.Mutil_GPU_Model(m, Parameter["leng"], Parameter["Channel"], lenMS, strategy)
        model, Loss_Acc = ddaalg.Model_Train_Predict(model,
                                                    trainX_tfcon,
                                                    trainY_tfcon,
                                                    trainX_tfcon[value:],
                                                    trainY_tfcon[value:],
                                                    trainX_tfcon.shape[0])

        model.save(Parameter["ModelPath"] + "HumanDDA_" + m + ".model")
        np.save(Parameter["ModelPath"] + m + "_Train_Result.npy", Loss_Acc)
    # y_pred, accuracy = mla.MLClassifier('XGB', xtrain[:value], ytrain[:value], xtrain[value:], ytrain[value:], 'XGB.pkl')




