from Profile import GenomeProfile as Profile
from Window import GenomeWindow as Window
from WindowDataSet import WindowDataSet as WSet

###############################################################
# Object for managing and handling a group od Genome Profiles
class ProfileDataSet:

    # Public Class Properties
    DataSetName = "DEFAULT"
    ProfileDataSet = []

    AverageOccurrence = -1
    LargestProfileIndex = -1
    SmallestProfileIndex = -1

    DefaultRanksAmount = 5
    ProfileRankIndecies = []


    # Constructor called as first line of the input set is read (Only available data is the name)
    def __init__(self, SetName):
        print("Making dataset by name: " + SetName)
        self.DataSetName = SetName
        self.WindowSamples = WSet()

    def AddProfileByObject(self, Profile):
        self.ProfileDataSet.append(Profile)

    def DivideLineAndCreateProfiles(self, line):
        print("Dividing initial Line and making profiles")
        # Column count utilized for managing the first few columns
        columnCount = 0
        for ProfileID in line.split():
            # Excludes the first three columns that describe nomenclature
            if columnCount >= 3:
                newProfile = Profile(ProfileID)
                self.AddProfileByObject(newProfile)
            columnCount = columnCount + 1

    def AddWindowByLine(self, line, windowIndex):
        # print("Adding new potential window: " + str(windowIndex))
        newWindow = Window(line)
        if newWindow.isValidData():
            self.WindowSamples.addWindowByObject(newWindow)
            currentWindowIndex = self.WindowSamples.getLastIndex()
            # Loops through all of the matching profile indecies found and adds them to the profile object.
            for occuranceIndex in newWindow.profilePresenceIndecies:
                self.ProfileDataSet[occuranceIndex].occuranceInWindowFound(currentWindowIndex)
        else:
            print(self.DataSetName + ": Invalid Sample Window Found in data set!")

    def AddWindowByLineWithRangeFilter(self, line, windowIndex, rangeStart, RangeEnd):
        # print("Adding new potential window: " + str(windowIndex))
        newWindow = Window(line)
        if newWindow.isValidData():
            self.WindowSamples.addWindowByObject(newWindow)
            currentWindowIndex = self.WindowSamples.getLastIndex()
            # Loops through all of the matching profile indecies found and adds them to the profile object.
            for occuranceIndex in newWindow.profilePresenceIndecies:
                self.ProfileDataSet[occuranceIndex].occuranceInWindowFound(currentWindowIndex)
        else:
            print(self.DataSetName + ": Invalid Sample Window Found in data set!")

    def CalculateSummaryValues(self):
        self.AverageOccurrence = self.CalculateAverageOccurrence()
        self.LargestProfileIndex = self.FindLargestProfileIndex()
        self.SmallestProfileIndex = self.FindSmallestProfileIndex()
        self.CalculateProfileRanks(5)

    # WARNING! Screws up index tracking for sub components
    def SortSetBySize(self):
        self.ProfileDataSet.sort(key=Profile.GetProfileSize)

    def GetSetSize(self):
        return len(self.ProfileDataSet)

    def CalculateAverageOccurrence(self):
        occurrenceSum = 0
        for Profile in self.ProfileDataSet:
            occurrenceSum = occurrenceSum + Profile.GetProfileSize()
        return occurrenceSum /(len(self.ProfileDataSet))


    def CalculateProfileRanks(self, RankCount):
        # Find smallest index up here because it is used twice:
        smallestIndex = self.FindSmallestProfileIndex()
        largestIndex = self.FindLargestProfileIndex()
        # Finds Range of data set Profile Size and divides by RankAmount
        ZonesSize = (self.ProfileDataSet[largestIndex].occuranceAmount)/RankCount
        self.ProfileRankIndecies = []
        x = 0
        # Loops through and appends empty arrays to the rank counts array
        while x < RankCount:
            newRankArray = []
            self.ProfileRankIndecies.append(newRankArray)
            x = x + 1
        # Loops through all of the Profiles Stored and assigns them ranks
        CurrentPIndex = 0
        SetCount = self.GetSetSize()
        while CurrentPIndex < SetCount:
            # For each profile, loop through until it qualifies for a certain rank
            rank = 1
            while (rank * ZonesSize) < self.ProfileDataSet[CurrentPIndex].GetProfileSize() and rank <= len(self.ProfileRankIndecies):
                # Increment through all the ranks until the occurance amount is less then the current rank's limit
                rank = rank + 1
            #print("Rank of profile found:  " + str(rank))
            self.ProfileRankIndecies[rank - 1].append(CurrentPIndex)
            CurrentPIndex = CurrentPIndex + 1

    def FindLargestProfileIndex(self):
        largestIndex = 0
        x = 0
        size = self.GetSetSize()
        while x < size:
            if self.ProfileDataSet[x].GetProfileSize() > self.ProfileDataSet[largestIndex].GetProfileSize():
                largestIndex = x
            x = x + 1
        return largestIndex

    def FindSmallestProfileIndex(self):
        smallestIndex = 0
        x = 0
        size = self.GetSetSize()
        while x < size:
            if self.ProfileDataSet[x].GetProfileSize() < self.ProfileDataSet[smallestIndex].GetProfileSize():
                smallestIndex = x
            x = x + 1
        return smallestIndex






    # Method for easy display of Sample Window Properties
    def printProfileDataSetSummaryInTerminal(self):
        if self.LargestProfileIndex == -1 or self.LargestProfileIndex == -1 or self.AverageOccurrence == 0:
            self.CalculateSummaryValues()
        print("")
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print("Data Summary for " + self.DataSetName)
        print("Profiles Available: " + str(self.GetSetSize()))
        print("Window Sample Count: " + str(self.WindowSamples.GetSetSize()))
        print("Average occurrences Detected Per NP: " + str(self.AverageOccurrence))
        print("Largest Profile Found: " + self.ProfileDataSet[self.LargestProfileIndex].profileID + " with: " + str(self.ProfileDataSet[self.LargestProfileIndex].occuranceAmount) + " sample occurrences")
        print("Smallest profile found: " + self.ProfileDataSet[self.SmallestProfileIndex].profileID + " with: " + str(self.ProfileDataSet[self.SmallestProfileIndex].occuranceAmount) + " sample occurrences")

        # Loops through each rank and prints the count of that rank in the terminal
        print("Profile Rank Count:___________")
        currentRank = 1
        for rankArray in self.ProfileRankIndecies:
            print("Rank " + str(currentRank) + ": " + str(len(rankArray)))
            currentRank = currentRank + 1
        print("______________________________________________")
        print("")

        # Old data printing from week one and two:
        # print("Samples Rank 1/2/3/4/5: " + str(SamplesRank1) + "/" + str(SamplesRank2) + "/" + str(SamplesRank3) + "/" + str(SamplesRank4) + "/" + str(SamplesRank5))
        # validSamples[smallestValidSampleIndex].printSample()
        # Find window for 21.7Mb -> 24.1Mb
