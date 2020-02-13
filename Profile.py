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

    def JaccardSimlarityIndexUnNormalized(self, otherNP):
        matches = 0
        for windowIndex in self.IndeciesOfValidWindowSamples:
            for otherProfileWindowIndex in otherNP.IndeciesOfValidWindowSamples:
                if windowIndex == otherProfileWindowIndex:
                    matches = matches + 1
                    break
                elif otherProfileWindowIndex > windowIndex:
                    break
        print("Jaccard Score: " + str(matches/len(self.IndeciesOfValidWindowSamples)))
        return matches/len(self.IndeciesOfValidWindowSamples)

    def JaccardSimlarityIndexNormalized(self, otherNP):
        matches = 0
        if self != otherNP:
            if len(self.IndeciesOfValidWindowSamples) != 0 and len(otherNP.IndeciesOfValidWindowSamples) != 0:
                for windowIndex in self.IndeciesOfValidWindowSamples:
                    for otherProfileWindowIndex in self.IndeciesOfValidWindowSamples:
                        if windowIndex == otherProfileWindowIndex:
                            matches = matches + 1
                            break
                        elif otherProfileWindowIndex > windowIndex:
                            break
                # print("Jaccard Score: " + str(matches/len(self.IndeciesOfValidWindowSamples)))
                if len(self.IndeciesOfValidWindowSamples) < len(otherNP.IndeciesOfValidWindowSamples):
                    return matches/len(self.IndeciesOfValidWindowSamples)
                else:
                    return matches/len(otherNP.IndeciesOfValidWindowSamples)
            else:
                return 0
        else:
            return 1

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
