from Window import GenomeWindow as WS
import pandas as pd
from DataVisualiser import DataVisualiser as DataVis

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

    @staticmethod
    def print2D_arry(array):
        print("")
        print("")
        for sA in array:
            rowS = ""
            for item in sA:
                rowS = rowS + str(item) + " "
            print(rowS)

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

    def getSortedCentralityList(self):
        newlist = sorted(self.WindowDataSet, key=lambda x: x.DegreeOfCentrality, reverse=True)
        return newlist

    def getSetCommunityAmountAndDisplayResults(self, communityCount):
        print("")
        print("")
        print("Generating " + str(communityCount) + " communities....")
        orderedList = self.getSortedCentralityList()
        currentCom = 1
        GeneratedComs = []
        nextCommunity = []
        while currentCom <= communityCount:
            # print("Generating Community # : " + str(currentCom + 1))
            nextCommunity = self.getCommunitiesSetByCenterNodeAndPrintSummaryData(orderedList[currentCom])
            # print("Generated community length: " + str(len(nextCommunity)))
            # print("Displaying Community 2D visualisation....")
            titleString = "Community #" + str(currentCom) + ":"
            currentCom = currentCom + 1
            DataVis.GenerateBinaryCommunityHeatMap(titleString, "Windows", "Windows", nextCommunity, self.GetAllWindowIndecies())
        # Resets current coms to


    def findNeighborsFromLinkageTable(self):
        gridDimensions = len(self.LinkageDisequilibriumMatrix)
        x = 0
        while x < gridDimensions:
            y = 0
            while y < gridDimensions:
                if self.LinkageDisequilibriumMatrix[x][y] > self.MeanCentralityOfWindows:
                    self.WindowDataSet[x].NeighborIndecies.append(y)
                y = y + 1
            x = x + 1



    def getCommunitiesSetByCenterNodeAndPrintSummaryData(self, centerNode):
        self.findNeighborsFromLinkageTable()
        communityNodeMembers = []
        # Adds the new node to the list of current members
        # Defines the grid used to show the community
        newCommunity = []

        # Summary data tracked:
        communitySize = 0
        communityHistOneCount = 0
        communityLADCount = 0

        # Calculates the dimensions of the community grid
        communityDimensions = len(self.WindowDataSet)
        x = 1
        while x < communityDimensions:
            newCommunityRow = []
            y = 1
            while y < communityDimensions:
                # Checks if the current node has a shared occurance of a community node member
                if y in centerNode.NeighborIndecies or x in centerNode.NeighborIndecies:
                    newCommunityRow.append(1)
                    # Updates summary totals:
                    communitySize = communitySize + 1
                    if self.WindowDataSet[y].LADPresence == 1:
                        communityLADCount = communityLADCount + 1
                    if self.WindowDataSet[y].HIST1_Presence == 1:
                        communityHistOneCount = communityHistOneCount + 1
                    if y not in communityNodeMembers:
                        communityNodeMembers.append(y)
                else:
                    newCommunityRow.append(0)
                y = y + 1
            newCommunity.append(newCommunityRow)
            x = x + 1
        # Here we print out summary results before returning the 2d visualsation grid:
        # WindowDataSet.print2D_arry(newCommunity)
        print("Community Summary Data:")
        print("Community Size - " + str(len(newCommunity)) + "X" + str(len(newCommunity[0])))
        print("Community Count - " + str(communitySize))
        print("Community Percentage Hist1 - " + str(communityHistOneCount/communitySize))
        print("Community Percentage LAD - " + str(communityLADCount/communitySize))
        print("Community Member List: ")
        comListString = ""
        memberNum = 1
        z = 0
        while z < len(communityNodeMembers):
            comListString = comListString + "Member #" + str(memberNum) + ":" + str(self.WindowDataSet[communityNodeMembers[z]].rowStart) + "/" + str(self.WindowDataSet[communityNodeMembers[z]].rowEnd) + "; "
            z = z + 1
            memberNum = memberNum + 1
        print(comListString)
        return newCommunity




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


