from Window import GenomeWindow as WS
import pandas as pd

###############################################################
# Object for managing and handling Genome Sample Window
class WindowDataSet:

    dataSetLocationDirectory = "Dataset"
    featureSpreadsheetOne = "gam_feature_community.csv"
    # Public Class Properties
    SetCount = 0

    # Constructor takes a row in the file as an input
    def __init__(self):
        self.SetCount = 0
        self.WindowDataSet = []

    # Helper method for validating the window of the sample
    def isValidData(self, window):
        # print("Window Size : " + str(self.windowSize) + '\n')
        if self.window.oneCount >= 100:
            return False
        else:
            return True

    def addWindowByLine(self, line):
        # Plan on filtering Data here
        # No filtering done for windows
        NewWindow = WS(line)
        self.WindowDataSet.append(NewWindow)
        self.SetCount = self.SetCount + 1

    def GetSetSize(self):
        return len(self.WindowDataSet)

    def addWindowByObject(self, NewWindow):
        self.WindowDataSet.append(NewWindow)
        self.SetCount = self.SetCount + 1

    def getLastIndex(self):
        toReturn = (self.SetCount - 1)
        return toReturn

    def analyzeDataSetAndFindMappings(self):
        filePath = dataSetLocationDirectory + featureSpreadsheetOne
        DataSourceFile = pd.read_excel(filePath)
        for window in self.WindowDataSet:
            if DataSourceFile['']

