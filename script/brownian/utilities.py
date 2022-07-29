

'''
utility function for tug-of-war project

'''
import os
import json

def getDefaultPrm():
    path = os.path.dirname(__file__)
    with open(path +"/prm.json", "r") as fp:
        default_parameter = json.load(fp)
    return default_parameter


def printPrm(prm):
    print("prm info:")
    for item, amount in prm.items():
        print("{}:\t{}".format(item, amount))
    return


def printProgress(n, N):
    percent = int(100.0*n/N)
    toPrint = "progress: "
    for i in range(percent//5):
        toPrint += '|'
    toPrint += "{:d}%".format(percent)
    print(toPrint, end='\r')
    return