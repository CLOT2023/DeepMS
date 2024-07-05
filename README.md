# DeepMS
The entire code involved in the article called 'DeepMS: Real-time MS Spectra Identification by Deep Learning' is here.

1.Configure the computer operating environment
Ubuntu 20.04 OR Ubuntu 18.04
CUDA v11.2
CUDNN v8.1.1
Python v3.8
keras v2.6.0
tensorflow-gpu v2.6.0
numpy v1.19.2
2.File description
Core code package installation command:pip install py_DDA_DL-1.0-py3-none-any.whl
Data preprocessing:
SaveData.py
SaveDataTest.py
Read_MGF.py
Read_MGFTest.py
Deep learning training code:
main.py
Test code
main_Test.py
3.Data description
Test dataset database:ProteomeXchange
PXD Number:PXD041292;PXD007863;PXD050120;PXD036650;PXD034913;PXD035780;PXD036800;

*:One can get all test data about this project through visiting database named ProteomeXchange.Then search number mentioned above as well as get data you are interested.Next,one must perform pFind v3.1 to deal with raw data.We provide a reference library for your reference.
Reference Library:(DB:Uniprot  DATA ID:UP000005640_9606  Species：Homo sapiens)
