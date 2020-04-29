###############################################################
# Object for managing and handling Genome Sample Window
class GenomeWindow:

    # Public Class Properties
    sampleID = ""
    windowSize = 0
    rank = -1
    rowStart = -1
    rowEnd = -1
    zeroCount = 0
    oneCount = 0
    HistOnePresence = False
    LADPresence = 0
    HIST1_Presence = 0
    profilePresenceIndecies = []
    DegreeOfCentrality = 0
    NeighborIndecies = []


    # Constructor takes a row in the file as an input
    def __init__(self, line):
        self.profilePresenceIndecies = []
        self.LADPresence = 0
        self.windowSize = 0
        self.HIST1_Presence = 0

        self.NeighborIndecies = []

        currentProfileIndex = 0
        sampleColumn = 0
        # Handles the first three columns respective data
        for value in line.split():
            if sampleColumn == 0:
                self.sampleID = value
            else:
                if sampleColumn == 1:
                    self.rowStart = value
                elif sampleColumn == 2:
                    self.rowEnd = value
                else: # Starts the actual occurrence intake of the window
                    if value == "1":
                        self.oneCount = self.oneCount + 1
                        self.windowSize = self.windowSize + 1
                        self.profilePresenceIndecies.append(currentProfileIndex)
                        currentProfileIndex = currentProfileIndex + 1
                        # print("Occur Found on index: " + str(currentProfileIndex))
                    elif value == "0":
                        self.zeroCount = self.zeroCount + 1
                        currentProfileIndex = currentProfileIndex + 1
                    else:
                        print ("Error! unexpected Value: " + value)
            sampleColumn = sampleColumn + 1

    # Method for easy display of Sample Window Properties
    def printSampleWindowSummary(self):
        print("ID: " + self.sampleID)
        print("Sample Start/End: " + str(self.rowStart) + "/" + str(self.rowEnd))
        print("Sample Window Size: " + str(self.windowSize))
        # print("Rank: " + str(self.findRank()))
        print("Degree of Centrality: " + str(self.DegreeOfCentrality))
        print("")

    # Helper method for validating the window of the sample
    def isValidData(self):
        if self.windowSize >= 100 or self.windowSize == 0:
            return False
        else:
            return True

    def printBothOccuranceLists(self, otherList):
        biggestList = self.profilePresenceIndecies
        smallestList = otherList
        if len(otherList) > len(self.profilePresenceIndecies):
            smallestList = self.profilePresenceIndecies
            biggestList = otherList
        index = 0
        print("Occurance List One:       Occurance list Two: ")
        while index < len(biggestList):
            currentSmallListString = ""
            if index < len(smallestList):
                currentSmallListString = "   " + str(smallestList[index])
            print(str(biggestList[index]) + currentSmallListString)
            index = index + 1


    # Helper method defining if the nodes are neighbors with one another (Returns a 1 for yes 0 for no):
    # Checks if the current node has a shared occurance of a community node member
    def isNeighbor(self, centerNode):
        isNeighbor = 0
        # self.printBothOccuranceLists(validOccurances)
        # cycle = 0
        if str(round(self.DegreeOfCentrality, 1)) == str(round(centerNode.DegreeOfCentrality, 1)):
            isNeighbor = 1
            # print("Neighbor found! In: " + self.sampleID + " my centrality " + str(round(self.DegreeOfCentrality, 1)) + " center centrality: " + str(round(centerNode.DegreeOfCentrality)))
        # for occurance in validOccurances:
        #     for myOccurance in self.profilePresenceIndecies:
        #         if occurance == myOccurance:
        #             # print("Neighbor found! In: " + self.sampleID + " my occurance: " + str(myOccurance) + " center occurances: " + str(occurance) + " CYCLE: " + str(cycle))
        #             isNeighbor = 1
        #             return isNeighbor
        #         cycle = cycle + 1
        return isNeighbor


    # Window Class:
    # Helper method that returns the cosegregation value
    def findCosegregation(self, otherWindow, ProfileCount):
        occuranceMatches = 0
        maxProfileOccurances = self.profileOA()
        if otherWindow.profileOA() > maxProfileOccurances:
            maxProfileOccurances = otherWindow.profileOA()
        for profile in self.profilePresenceIndecies:
            if profile in otherWindow.profilePresenceIndecies:
                occuranceMatches = occuranceMatches + 1
        # FORMULA FOR LINKAGE MATRICIES:
        result = (occuranceMatches/maxProfileOccurances)#  - (self.profileOA() * otherWindow.profileOA())
        return result

    # Helper method to get window length (Profile occurance ammount)
    def profileOA(self):
        return len(self.profilePresenceIndecies)


    # def FindCentrality(self, otherWindow, WindowCount):
    #     occuranceMatches = 0
    #     maxProfileOccurances = self.profileOA()
    #     if otherWindow.profileOA() > maxProfileOccurances:
    #         maxProfileOccurances = otherWindow.profileOA()
    #     for profile in self.profilePresenceIndecies:
    #         if profile in otherWindow.profilePresenceIndecies:
    #             occuranceMatches = occuranceMatches + 1
    #     # FORMULA FOR Normalized Linkage LINKAGE MATRICIES:
    #     self.DegreeOfCentrality = (occuranceMatches / maxProfileOccurances)/ (WindowCount - 1)





    # Boolean operators for comparing the windows to one another:
    def __gt__(self, other):
        if self.windowSize > other.windowSize:
            return True
        else:
            return False

    # Boolean operators for comparing the windows to one another:
    def __lt__(self, other):
        if self.windowSize < other.windowSize:
            return True
        else:
            return False
