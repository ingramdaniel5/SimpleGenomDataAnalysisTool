###############################################################
# Object for managing and handling Genome Profiles
class GenomeProfile:

    # Public Class Properties
    profileID = ""
    occuranceAmount = 0
    IndeciesOfValidWindowSamples = []

    # Constructor called as first line of the input set is read (Only available data is the name)
    def __init__(self, profileID):
        # print("Making Profile: " + profileID)
        self.profileID = profileID
        self.IndeciesOfValidWindowSamples = []

    # Helper method called to handle the detection of a 1 in a window for this column
    def occuranceInWindowFound(self, windowIndex):
        self.occuranceAmount = self.occuranceAmount + 1
        self.IndeciesOfValidWindowSamples.append(windowIndex)

    def GetProfileSize(self):
        return self.occuranceAmount

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
