from Window import GenomeWindow as WS
###############################################################
# Object for managing and handling Genome Sample Window
class WindowDataSet:

    # Public Class Properties
    WindowDataSet = []
    SetCount = 0


    # Constructor takes a row in the file as an input
    def __init__(self):
        self.SetCount = 0

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

    def addWindowByObject(self, NewWindow):
        self.WindowDataSet.append(NewWindow)
        self.SetCount = self.SetCount + 1

    def getLastIndex(self):
        toReturn = (self.SetCount - 1)
        return toReturn
