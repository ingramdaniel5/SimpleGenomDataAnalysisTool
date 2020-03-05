###############################################################
# Object for managing and handling Genome Profiles
class GenomeProfile:

    # Public Class Properties
    profileID = ""
    occuranceAmount = 0
    IndeciesOfValidWindowSamples = []
    HistOneWindowOccurancePercentage = 0
    LADWindowOccurancePercentage = 0
    rank = 0
    # Constructor called as first line of the input set is read (Only available data is the name)
    def __init__(self, profileID):
        # print("Making Profile: " + profileID)
        self.profileID = profileID
        self.IndeciesOfValidWindowSamples = []
        self.HistOneWindowOccurancePercentage = 0
        self.LADWindowOccurancePercentage = 0
        self.rank = 0;

    # Helper method called to handle the detection of a 1 in a window for this column
    def occuranceInWindowFound(self, windowIndex):
        self.occuranceAmount = self.occuranceAmount + 1
        self.IndeciesOfValidWindowSamples.append(windowIndex)

    def GetProfileSize(self):
        return self.occuranceAmount

    @staticmethod
    def ProfileIsEmpty(self):
        if len(self.IndeciesOfValidWindowSamples) != 0:
            return True
        else:
            return False


    def JaccardSimlarityIndexNormalized(self, otherNP):
        # Determines which of the two Profiles to use as the baseline based on which is smallest:
        if len(self.IndeciesOfValidWindowSamples) != 0 and len(otherNP.IndeciesOfValidWindowSamples) != 0:
            if len(self.IndeciesOfValidWindowSamples) >= len(otherNP.IndeciesOfValidWindowSamples):
                chosenSampleBase = otherNP.IndeciesOfValidWindowSamples
                chosenSampleCompare = self.IndeciesOfValidWindowSamples
            else:
                chosenSampleBase = self.IndeciesOfValidWindowSamples
                chosenSampleCompare = otherNP.IndeciesOfValidWindowSamples

            # Generates order of Union off of that:
            intersection = len(list(set(chosenSampleBase).intersection(chosenSampleCompare)))
            # union = (len(chosenSampleBase) + len(chosenSampleCompare)) - intersection
            return float(intersection) / len(chosenSampleBase)
        elif self == otherNP:
            return 1
        else:
            return 0

     # Boolean operators for comparing the profiles to one another:
    def __gt__(self, other):
        if self.occuranceAmount > other.occuranceAmount:
            return True
        else:
            return False

    # Boolean operators for comparing the profiles to one another:
    def __lt__(self, other):
        if self.occuranceAmount < other.occuranceAmount:
            return True
        else:
            return False
