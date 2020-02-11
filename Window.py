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
    profilePresenceIndecies = []

    # Constructor takes a row in the file as an input
    def __init__(self, line):
        self.profilePresenceIndecies = []
        self.windowSize = 0

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
        # self.findRank()
        # print("Indecies Array length = " + str(len(self.profilePresenceIndecies)))
        # print("Window Row Size Count: " + str(self.oneCount))

    # Method for easy display of Sample Window Properties
    def printSampleWindowSummary(self):
        print("ID: " + self.sampleID)
        print("Sample Start/End: " + str(self.rowStart) + "/" + str(self.rowEnd))
        print("Sample Window Size: " + str(self.windowSize))
        print("Rank: " + str(self.findRank()))
        print("")

    # Helper method for validating the window of the sample
    def isValidData(self):
        if self.windowSize >= 100:
            return False
        else:
            return True

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
