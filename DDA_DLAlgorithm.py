import math
from time import time
import ResNet
import DenseNet
import VGG_Google
import Transformer
from tensorflow.keras.optimizers import Adam, RMSprop, Nadam
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau, LearningRateScheduler


def Mutil_GPU_Model(namestr, shape1, shape2, Output, strategy):
    global model
    with strategy.scope():
        if namestr == "GoogleNet":
            model = VGG_Google.GoogLeNet_1D(shape1, shape2, Output)

        elif namestr == "VGG16":
            model = VGG_Google.VGG16_1D(shape1, shape2, Output)

        elif namestr == "VGG19":
            model = VGG_Google.VGG19_1D(shape1, shape2, Output)

        elif namestr == "AlexNet":
            model = VGG_Google.AlexNet_1D(shape1, shape2, Output)

        elif namestr == "LSTM":
            model = VGG_Google.MS_LSTM(shape1, shape2, Output)

        elif namestr == "GRU":
            model = VGG_Google.MS_GRU(shape1, shape2, Output)

        elif namestr == "ResNet50":
            model = ResNet.ResNet_1D(shape1, shape2, Output, 50)

        elif namestr == "ResNet101":
            model = ResNet.ResNet_1D(shape1, shape2, Output, 101)

        elif namestr == "ResNet152":
            model = ResNet.ResNet_1D(shape1, shape2, Output, 152)

        elif namestr == "DenseNet121":
            model = DenseNet.DenseNet_1D(shape1, shape2, Output, depth=121)

        elif namestr == "DenseNet169":
            model = DenseNet.DenseNet_1D(shape1, shape2, Output, depth=169)

        elif namestr == "DenseNet201":
            model = DenseNet.DenseNet_1D(shape1, shape2, Output, depth=201)

        elif namestr == "DenseNet264":
            model = DenseNet.DenseNet_1D(shape1, shape2, Output, depth=264)

        elif namestr == "Transformer":
            model = Transformer.Transformer_1D(shape1, shape2, Output)

        # 可以尝试 Adam, RMSprop 或 Nadam
        model.compile(loss='binary_crossentropy', optimizer=Nadam(learning_rate=0.001), metrics=['accuracy'])
        # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        print("Mutil GPU Model : " + namestr)
        return model

def Fit_Para(lenMS):
    if 10000 < lenMS <= 15000:
        return [1200, 40, 128]
    elif 15000 < lenMS <= 20000:
        return [1600, 40, 128]
    elif 20000 < lenMS <= 40000:
        return [1600, 50, 128]
    elif 40000 < lenMS <= 60000:
        return [1600, 70, 256]
    elif 60000 < lenMS <= 80000:
        return [1600, 80, 256]
    elif 80000 < lenMS <= 100000:
        return [1600, 100, 256]
    elif 100000 < lenMS <= 500000:
        return [1600, 120, 256]
    elif 500000 < lenMS <= 2000000:
        return [1800, 150, 256]
    elif 2000000 < lenMS <= 8000000:
        return [2000, 200, 256]

# 余弦退火学习率调度器
def cosine_annealing(epoch, lr_max=0.001, lr_min=0.00001, T_max=50):
    lr = lr_min + (lr_max - lr_min) * (1 + math.cos(math.pi * epoch / T_max)) / 2
    return lr


def Model_Train_Predict(model, trainX, trainY, testX, testY, lenxtrain):
    # 回调函数
    early_stopping = EarlyStopping(monitor='accuracy', patience=70, restore_best_weights=True)
    reduce_lr = ReduceLROnPlateau(monitor='accuracy', factor=0.25, patience=20, min_lr=0.00001)
    lr_scheduler = LearningRateScheduler(cosine_annealing)
    Fit_Parameter = Fit_Para(lenxtrain)
    # begin_time = time()
    DDA_Model = model.fit(x=trainX, y=trainY, validation_data=(testX, testY),
                epochs=Fit_Parameter[0], steps_per_epoch=Fit_Parameter[1], batch_size=Fit_Parameter[2],
                verbose=2, shuffle=True, callbacks=[early_stopping, reduce_lr])
    # end_time = time()
    print("\n=======================\n\n")
    # print("Train Time = ", end_time - begin_time, "s")

    Loss_Acc = {}
    Loss_Acc['loss'] = DDA_Model.history['loss']
    Loss_Acc['accuracy'] = DDA_Model.history['accuracy']

    return model, Loss_Acc




