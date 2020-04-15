from Window import GenomeWindow as WS
import pandas as pd

###############################################################
# Object for managing and handling Genome Sample Window
class WindowDataSet:

    dataSetLocationDirectory = "Dataset"
    featureSpreadsheetOne = "gam_feature_community.csv"
    # Public Class Properties
    SetCount = 0

    LinkageDisequilibriumMatrixNormalised = []

    CentralityValues = []
    MeanCentralityOfWindows = -1
    MaxCentralityIndexOfWindows = -1
    MinCentralityIndexOfWindows = -1

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


    # Window Data Set class:
    # Used to fill out the 2D matrix item: LinkageDisequilibriumMatrix
    def CalculateLinkageMatricies(self, profileCount):
        currentRow = 0
        self.LinkageDisequilibriumMatrix = []
        for windowRow in self.WindowDataSet:
            newRow = []
            self.LinkageDisequilibriumMatrix.append(newRow)
            for windowColumn in self.WindowDataSet:
                self.LinkageDisequilibriumMatrix[currentRow].append(windowColumn.findCosegregation(windowRow, profileCount))
            currentRow = currentRow + 1


    def findAllCentralityValues(self):
        windowCount = len(self.WindowDataSet)
        self.MaxCentralityIndexOfWindows = -1
        self.MinCentralityIndexOfWindows = -1
        meanCentalitySum = 0
        x = 0
        while x < windowCount:
            currentCentalityValue = 0
            # Summing up all of the values
            for linkageValue in self.LinkageDisequilibriumMatrix[x]:
                currentCentalityValue = currentCentalityValue + linkageValue

            # Assigning the new linkage value to that window
            self.WindowDataSet[x].DegreeOfCentrality = linkageValue
            # Checking if value is a min or a max
            if self.MinCentralityIndexOfWindows == -1 or self.WindowDataSet[x].DegreeOfCentrality < self.WindowDataSet[self.MinCentralityIndexOfWindows].DegreeOfCentrality:
                self.MinCentralityIndexOfWindows = x
            if self.MaxCentralityIndexOfWindows == -1 or self.WindowDataSet[x].DegreeOfCentrality > self.WindowDataSet[self.MaxCentralityIndexOfWindows].DegreeOfCentrality:
                self.MaxCentralityIndexOfWindows = x
            meanCentalitySum = meanCentalitySum + self.WindowDataSet[x].DegreeOfCentrality
            #Incrementing X value
            x = x + 1
        self.MeanCentralityOfWindows = meanCentalitySum / windowCount


    # Helper Method to get an array of all of the window sample names:
    def GetAllWindowNames(self):
        newNameArray = []
        for window in self.WindowDataSet:
            newNameArray.append(window.sampleID)
        return  newNameArray

    # Helper Method to get an array of all of the window sample names:
    def GetAllWindowIndecies(self):
        newIndexArray = []
        index = 0
        for window in self.WindowDataSet:
            newIndexArray.append(index)
            index = index + 1
        return newIndexArray

    def printSortedCentralityList(self):
        newlist = sorted(self.WindowDataSet, key=lambda x: x.DegreeOfCentrality, reverse=True)
        for window in newlist:
            print ("Centrality: " + str(window.DegreeOfCentrality) + "  |  Window Start/End: " + str(window.rowStart) + "/" + str(window.rowEnd))

    def analyzeDataSetAndFindMappings(self):
        print("Mapping Windows to occurance in external data file...")
        dataFile = pd.read_csv('DataSets/gam_feature_community.csv')

        windowName = dataFile['name'].apply(lambda x: x.split(':')[0])
        windowStart = dataFile['name'].apply(lambda x: x.split(':')[1].split('-')[0])
        windowStop= dataFile['name'].apply(lambda x: x.split('-')[1])
        isInHistOne = dataFile['Hist1']
        isInLAD = dataFile['LAD']
        for x in range(len(isInHistOne)):
            if isInHistOne[x] == 1:
                # print("found lad occurance in " + str(x))
                self.WindowDataSet[x].LADPresence = 1
            if isInLAD[x] == 1:
                # print("found hist1 occurance in " + str(x))
                self.WindowDataSet[x].HIST1_Presence = 1


