from Profile import GenomeProfile as Profile
from Window import GenomeWindow as Window
from WindowDataSet import WindowDataSet as WSet
from HeatMapWindow import heatmap
import matplotlib.pyplot as plt
import random
# importing copy module
import copy

import numpy as np

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
    ProfileJaccardSimilarityMarix = []
    ProfileJaccardDifferenceMarix = []
    CurrentOptimalNPClusters = []


    # Constructor called as first line of the input set is read (Only available data is the name)
    def __init__(self, SetName):
        print("Making dataset by name: " + SetName)
        self.DataSetName = SetName
        self.WindowSamples = WSet()
        self.ProfileDataSet = []
        self.ProfileJaccardSimilarityMarix = []
        self.ProfileJaccardDifferenceMarix = []
        self.CurrentOptimalNPClusters = []

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

    def FindAllJaccardMatriciesNormalized(self):
        currentRow = 0
        self.ProfileJaccardSimilarityMarix = []
        self.ProfileJaccardDifferenceMarix = []
        for profileRow in self.ProfileDataSet:
            newSimRow = []
            newDifRow = []
            self.ProfileJaccardSimilarityMarix.append(newSimRow)
            self.ProfileJaccardDifferenceMarix.append(newDifRow)
            for profileColumn in self.ProfileDataSet:
                self.ProfileJaccardSimilarityMarix[currentRow].append(profileRow.JaccardSimlarityIndexNormalized(profileColumn))
                self.ProfileJaccardDifferenceMarix[currentRow].append(1 - profileRow.JaccardSimlarityIndexNormalized(profileColumn))
            currentRow = currentRow + 1
        print("Normalised Jaccard Similarities and Difference Matrix Found!")

    @staticmethod
    def getGroupAverage(group):
        summation = 0
        clusterSize = len(group)
        for profile in group:
            summation = summation + profile.GetProfileSize()
        return summation/clusterSize

    @staticmethod
    def findGroupMediodValue(cluster, meanValue):
        return

    @staticmethod
    def findGroupRange(cluster):
        largestProfileIndex = 0
        smallestProfileIndex = 0
        currentClusterIndex = 0
        clusterSize = len(cluster)
        while currentClusterIndex < clusterSize:
            if cluster[largestProfileIndex].GetProfileSize() <= cluster[currentClusterIndex].GetProfileSize():
                largestProfileIndex = currentClusterIndex
            if cluster[smallestProfileIndex].GetProfileSize() >= cluster[currentClusterIndex].GetProfileSize():
                smallestProfileIndex = currentClusterIndex
            cluster = cluster + 1
        return cluster[largestProfileIndex].GetProfileSize() - cluster[smallestProfileIndex].GetProfileSize()

    # Returns the index of the closest number in the array
    @staticmethod
    def getClosestIndex(group, TargetValue):
        group = np.asarray(group)
        idx = (np.abs(group - TargetValue)).argmin()
        return idx

    @staticmethod
    def getClosestValue(group, TargetValue):
        return group[min(range(len(group)), key=lambda i: abs(group[i] - TargetValue))]

    def findOptimalClustersInDataset(self, clusterAmount, RepAmount):
            print("Beginning Search for Optimal clustering based on given parameters...")
            # Initialization stuff
            self.CurrentOptimalNPClusters = []
            currentClustersAdded = 0
            # Adds the desired amount of clusters to the chart
            while currentClustersAdded < clusterAmount:
                newSubClusterArray = []
                self.CurrentOptimalNPClusters.append(newSubClusterArray)
                currentClustersAdded = currentClustersAdded + 1

            # Gets range for which to generate random starting values:
            smallestValuePossible = self.ProfileDataSet[self.FindSmallestProfileIndex()].GetProfileSize()
            largestValuePossible = self.ProfileDataSet[self.FindLargestProfileIndex()].GetProfileSize()

            # Gets the middle points at random initially:
            middleValuesGenerated = 0
            currentMiddleValues = []
            while middleValuesGenerated < clusterAmount:
                #  through appending random numbers:
                currentMiddleValues.append(random.randint(smallestValuePossible, largestValuePossible))
                middleValuesGenerated = middleValuesGenerated + 1

            # Beginning of optimal cluster finding:

            # Loops through desired times to try and find an optimal break down given the parameters
            currentRep = 0
            while currentRep < RepAmount:
                # Loops through appending to closest cluster AKA Makes Groups...
                for profile in self.ProfileDataSet:
                    optimalClusterIndex = ProfileDataSet.getClosestValue(currentMiddleValues, profile.GetProfileSize())
                    self.CurrentOptimalNPClusters[optimalClusterIndex].append(copy.deepcopy(profile))

                # Goes through and finds the middles of all clusters
                currentClustersAdded = 0
                while currentClustersAdded < clusterAmount:
                    cAverage = ProfileDataSet.getGroupAverage(self.CurrentOptimalNPClusters[currentClustersAdded])
                    currentMiddleValues[currentClustersAdded] = ProfileDataSet.findGroupMediodValue(self.CurrentOptimalNPClusters[currentClustersAdded], cAverage)
                    currentClustersAdded = currentClustersAdded + 1
                if currentRep + 1 != RepAmount:
                    # Before restarting clears the existing clusters if not on exit pass
                    self.CurrentOptimalNPClusters = []
                    clustersAdded = 0
                    # Adds the desired amount of clusters to the chart
                    while clustersAdded < clusterAmount:
                        newSubClusterArray = []
                        self.CurrentOptimalNPClusters.append(newSubClusterArray)
                        clusterAmount = clusterAmount + 1
                currentRep = currentRep + 1
            print("Attempt amount to find optimal clusters reached!")



    @staticmethod
    def print2DMatrix(matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                print(matrix[i][j], end=' ')
            print()


    def FilterEmptyProfiles(self):
        newDataSet = list(filter(Profile.ProfileIsEmpty, self.ProfileDataSet))
        self.ProfileDataSet = newDataSet.copy()
        print("Filtered empty profiles from dataset: " + self.DataSetName)


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

    def AddWindowByLineWithRangeFilter(self, line, windowIndex, RangeStart, RangeEnd, chrName):
        # print("Adding new potential window: " + str(windowIndex))
        newWindow = Window(line)
        if newWindow.isValidData() and int(newWindow.rowStart) >= RangeStart and int(newWindow.rowEnd) <= RangeEnd and newWindow.sampleID == chrName:
            self.WindowSamples.addWindowByObject(newWindow)
            currentWindowIndex = self.WindowSamples.getLastIndex()
            # Loops through all of the matching profile indecies found and adds them to the profile object.
            for occuranceIndex in newWindow.profilePresenceIndecies:
                self.ProfileDataSet[occuranceIndex].occuranceInWindowFound(currentWindowIndex)
        # else:
        #    print(self.DataSetName + ": Invalid Sample Window Found in data set!")

    def CalculateSummaryValues(self):
        self.AverageOccurrence = self.CalculateAverageOccurrence()
        self.LargestProfileIndex = self.FindLargestProfileIndex()
        self.SmallestProfileIndex = self.FindSmallestProfileIndex()
        self.CalculateProfileRanks(5)
        self.FindAllJaccardMatriciesNormalized()

    # WARNING! Screws up index tracking for sub components
    def SortSetBySize(self):
        self.ProfileDataSet.sort(key=Profile.GetProfileSize)

    def GetSetSize(self):
        return len(self.ProfileDataSet)

    def CalculateAverageOccurrence(self):
        occurrenceSum = 0
        for profile in self.ProfileDataSet:
            occurrenceSum = occurrenceSum + profile.GetProfileSize()
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

        # Displays all of the data for each of the calculated 2D Matrices:
    def printAllSummaryPlots(self):
        if len(self.ProfileJaccardSimilarityMarix) != 0:
            JS_display = plt.figure(num=1, figsize=(6, 6), dpi=150)
            JS_plot = JS_display.subplots()
            JS_plot.set_title('Jaccard Similarity')
            JS_plot.set_aspect('auto')
            JS_plot.set_autoscalex_on(True)
            plt.imshow(self.ProfileJaccardSimilarityMarix)
            plt.colorbar(orientation='vertical')

        if len(self.ProfileJaccardDifferenceMarix) != 0:
            JD_display = plt.figure(num=2, figsize=(6, 6), dpi=150)
            JD_plot = JD_display.subplots()
            JD_plot.set_title('Jaccard Difference')
            JD_plot.set_aspect('auto')
            JD_plot.set_autoscalex_on(True)
            plt.imshow(self.ProfileJaccardDifferenceMarix)
            plt.colorbar(orientation='vertical')
        plt.show()