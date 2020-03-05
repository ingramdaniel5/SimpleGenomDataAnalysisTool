from Profile import GenomeProfile as Profile
from Window import GenomeWindow as Window
from WindowDataSet import WindowDataSet as WSet
from HeatMapWindow import heatmap
import matplotlib.pyplot as plt
import random
# importing copy module
import copy
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
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


    def getGroupAverage(self, group):
        summation = 0
        clusterSize = len(group)
        for profileIndex in group:
            summation = summation + self.ProfileDataSet[profileIndex].GetProfileSize()
        return summation/clusterSize

    def findGroupMediodValue(self, cluster, meanValue):
        smallestProfileSizeDeviationIndex = 0
        currentClusterProfileIndex = 0
        groupSize = len(cluster)
        while currentClusterProfileIndex < groupSize:
            if abs(self.ProfileDataSet[cluster[currentClusterProfileIndex]].GetProfileSize() - meanValue) < abs(
                    self.ProfileDataSet[cluster[smallestProfileSizeDeviationIndex]].GetProfileSize() - meanValue):
                smallestProfileSizeDeviationIndex = currentClusterProfileIndex
            currentClusterProfileIndex = currentClusterProfileIndex + 1
        return self.ProfileDataSet[cluster[smallestProfileSizeDeviationIndex]].GetProfileSize()

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
        smallestProfileSizeDeviationIndex = 0
        currentClusterProfileIndex = 0
        groupSize = len(group)
        while currentClusterProfileIndex < groupSize:
            if abs(group[currentClusterProfileIndex].GetProfileSize() - TargetValue) < abs(group[smallestProfileSizeDeviationIndex].GetProfileSize() - TargetValue):
                smallestProfileSizeDeviationIndex = currentClusterProfileIndex
            currentClusterProfileIndex = currentClusterProfileIndex + 1
        return smallestProfileSizeDeviationIndex

    @staticmethod
    def getClosestIndex_NON_PROFILE_GROUP(group, TargetValue):
        smallestProfileSizeDeviationIndex = 0
        currentIndex = 0
        groupSize = len(group)
        while currentIndex < groupSize:
            if abs(group[currentIndex] - TargetValue) < abs(group[smallestProfileSizeDeviationIndex] - TargetValue):
                smallestProfileSizeDeviationIndex = currentIndex
            currentIndex = currentIndex + 1
        return smallestProfileSizeDeviationIndex

    @staticmethod
    def getClosestValue(group, TargetValue):
        smallestProfileSizeDeviationIndex = 0
        currentClusterProfileIndex = 0
        groupSize = len(group)
        while currentClusterProfileIndex < groupSize:
            if abs(group[currentClusterProfileIndex].GetProfileSize() - TargetValue) < abs(group[smallestProfileSizeDeviationIndex].GetProfileSize() - TargetValue):
                smallestProfileSizeDeviationIndex = currentClusterProfileIndex
            currentClusterProfileIndex = currentClusterProfileIndex + 1
        return group[smallestProfileSizeDeviationIndex]

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
                print("Cluster Optimizaition: " + str(currentRep) + "/" + str(RepAmount))
                # Loops through appending to closest cluster AKA Makes Groups...
                currentProfile = 0
                for profile in self.ProfileDataSet:
                    targetValue = profile.GetProfileSize()
                    optimalClusterIndex = ProfileDataSet.getClosestIndex_NON_PROFILE_GROUP(currentMiddleValues, targetValue)
                    self.CurrentOptimalNPClusters[optimalClusterIndex].append(currentProfile)
                    currentProfile = currentProfile + 1

                # Goes through and finds the middles of all clusters
                currentMiddleValues = []
                currentMiddleValueReplaced = 0
                while currentMiddleValueReplaced < clusterAmount:
                    cAverage = self.getGroupAverage(self.CurrentOptimalNPClusters[currentMiddleValueReplaced])
                    newCenterMediodValue = self.findGroupMediodValue(self.CurrentOptimalNPClusters[currentMiddleValueReplaced], cAverage)
                    currentMiddleValues.append(newCenterMediodValue)
                    currentMiddleValueReplaced = currentMiddleValueReplaced + 1

                    # Resets all of the clusters for another optomization:
                if currentRep + 1 != RepAmount:
                    # Before restarting clears the existing clusters if not on exit pass
                    self.CurrentOptimalNPClusters = []
                    clustersAdded = 0
                    # Adds the desired amount of clusters to the chart
                    while clustersAdded < clusterAmount:
                        newSubClusterArray = []
                        self.CurrentOptimalNPClusters.append(newSubClusterArray)
                        clustersAdded = clustersAdded + 1
                currentRep = currentRep + 1
            print("Attempt amount to find optimal clusters reached!")

    @staticmethod
    def print2DMatrix(matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                print(matrix[i][j].profileID, end=' ')
            print('\n\n')


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
        self.WindowSamples.analyzeDataSetAndFindMappings()
        self.findAllProfileSamplePercentagesForExternalFile()

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
            # print("Rank of profile found:  " + str(rank))
            self.ProfileRankIndecies[rank - 1].append(CurrentPIndex)
            self.ProfileDataSet[CurrentPIndex].rank = rank
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

    def findAllProfileSamplePercentagesForExternalFile(self):
        for x in range(len(self.ProfileDataSet)):
            HistSum = 0
            LadSum = 0
            for windowIndex in self.ProfileDataSet[x].IndeciesOfValidWindowSamples:
                if self.WindowSamples.WindowDataSet[windowIndex].LADPresence == 1:
                    LadSum = LadSum + 1
                if self.WindowSamples.WindowDataSet[windowIndex].HIST1_Presence == 1:
                    HistSum = HistSum + 1
            currentSetWindowTotal = len(self.ProfileDataSet[x].IndeciesOfValidWindowSamples)
            self.ProfileDataSet[x].HistOneWindowOccurancePercentage = HistSum / currentSetWindowTotal
            self.ProfileDataSet[x].LADWindowOccurancePercentage = LadSum / currentSetWindowTotal
            # print ("Hist percent = " + str(self.ProfileDataSet[x].HistOneWindowOccurancePercentage))
            # print("LAD percent = " + str(self.ProfileDataSet[x].LADWindowOccurancePercentage))

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
        IndeciesByClusterOrder = []
        allNamesInClusters = []
        # LAD DATA
        LADallClusterPercentages = []
        LADcombinedClusterPercentages = []
        #HIST1 DATA
        HISTallClusterPercentages = []
        HISTcombinedClusterPercentages = []

        # Radial to Cluster
        ClusterAllRadialRankings = []

        #ClusterJacquardMatricies:
        ClusterAllJacquardMatricies = []

        AllNames = []
        AllJacqs = []

        for cluster in self.CurrentOptimalNPClusters:
            newHISTClusterPercentage = []
            newLADClusterPercentage = []
            newClusterRadRank = []
            newJacquardMatricies = []
            newClusterNames = []

            for profileIndex in cluster:
                newHISTClusterPercentage.append(self.ProfileDataSet[profileIndex].HistOneWindowOccurancePercentage)
                HISTcombinedClusterPercentages.append(self.ProfileDataSet[profileIndex].HistOneWindowOccurancePercentage)
                newLADClusterPercentage.append(self.ProfileDataSet[profileIndex].LADWindowOccurancePercentage)
                LADcombinedClusterPercentages.append(self.ProfileDataSet[profileIndex].LADWindowOccurancePercentage)
                IndeciesByClusterOrder.append(profileIndex)
                newClusterRadRank.append(self.ProfileDataSet[profileIndex].rank)
                newJacquardMatricies.append(self.ProfileJaccardSimilarityMarix[profileIndex])
                newClusterNames.append(self.ProfileDataSet[profileIndex].profileID)
                AllNames.append(self.ProfileDataSet[profileIndex].profileID)
                AllJacqs.append(self.ProfileJaccardSimilarityMarix[profileIndex])
            HISTallClusterPercentages.append(newHISTClusterPercentage)
            LADallClusterPercentages.append(newLADClusterPercentage)
            ClusterAllRadialRankings.append(newClusterRadRank)
            ClusterAllJacquardMatricies.append(newJacquardMatricies)
            allNamesInClusters.append(newClusterNames)

        # Cluster Heatmap plotting individual:
        fig3 = make_subplots(rows=3, cols=1, vertical_spacing=0.05, shared_xaxes=False,
                            specs=[[{}],
                                   [{}],
                                   [{}]])
        fig3.add_trace(go.Heatmap( y=allNamesInClusters[0],x=allNamesInClusters[0], z=ClusterAllJacquardMatricies[0], colorscale='Spectral',reversescale=True, showscale=False), 1, 1)
        fig3.add_trace(go.Heatmap(y=allNamesInClusters[1], x=allNamesInClusters[1], z=ClusterAllJacquardMatricies[1], colorscale='Spectral', reversescale=True, showscale=False), 2, 1)
        fig3.add_trace(go.Heatmap(y=allNamesInClusters[2], x=allNamesInClusters[2], z=ClusterAllJacquardMatricies[2], colorscale='Spectral', reversescale=True, showscale=False), 3, 1)
        fig3.show()

        # Cluster Heatmap plotting individual:
        fig4 = make_subplots(rows=1, cols=1, vertical_spacing=0.05, shared_xaxes=False,
                             specs=[[{}]])
        fig4.add_trace(go.Heatmap(y=AllNames, x=AllNames, z=AllJacqs,
                                  colorscale='Spectral', reversescale=True, showscale=False), 1, 1)
        fig4.show()
        # HIST1 PLOTTING
        fig = make_subplots(rows=2, cols=1, vertical_spacing=0.1,
                            subplot_titles=('NPs vs Feature by Cluster for HIST1 and LAD',),
                            specs=[[{}], [{}]])
        clrs = px.colors.sequential.Viridis
        fig.add_trace(go.Box(y=HISTcombinedClusterPercentages, name='Hist1 Total', marker_color='royalblue', boxpoints='all'), 1, 1)
        fig.add_trace(
        go.Box(y=HISTallClusterPercentages[0], name='Hist1 (cluster 1)', marker_color=clrs[1],
               boxpoints='all'), 1, 1)
        fig.add_trace(
         go.Box(y=HISTallClusterPercentages[1], name='Hist1 (cluster 2)', marker_color=clrs[3],
                boxpoints='all'), 1, 1)
        fig.add_trace(
         go.Box(y=HISTallClusterPercentages[2], name='Hist1 (cluster 3)', marker_color=clrs[5],
                boxpoints='all'), 1, 1)

        # LAD1 PLOTTING
        fig.add_trace(go.Box(y=LADcombinedClusterPercentages, name='LAD Total', marker_color='royalblue', boxpoints='all'), 2, 1)
        fig.add_trace(go.Box(y=LADallClusterPercentages[0], name='LAD (cluster 1)', marker_color=clrs[1],
                             boxpoints='all'), 2, 1)
        fig.add_trace(go.Box(y=LADallClusterPercentages[1], name='LAD (cluster 2)', marker_color=clrs[3],
                             boxpoints='all'), 2, 1)
        fig.add_trace(go.Box(y=LADallClusterPercentages[2], name='LAD (cluster 3)', marker_color=clrs[5],
                             boxpoints='all'), 2, 1)
        fig.update_layout(height=600, title_text="<b>Box Plots: Feature Table</b>", showlegend=False,
                          yaxis_tickformat='%')
        fig.update_yaxes(title_text="Percent Hist1 Content", row=1, col=1)

        fig.update_yaxes(title_text="Percent LAD Content", row=2, col=1, tickformat='%')
        fig.update_xaxes(title_text="Clusters", row=2, col=1)
        fig.show()

        # Radial plotting based on clusters:
        fig2 = go.Figure()
        fig2.add_trace(go.Violin(y=ClusterAllRadialRankings[0], name='Cluster 1', marker_color=clrs[1],
                                    points='all', box_visible=True, meanline_visible=True))
        fig2.add_trace(go.Violin(y=ClusterAllRadialRankings[1], name='Cluster 2', marker_color=clrs[3],
                                points='all', box_visible=True, meanline_visible=True))
        fig2.add_trace(go.Violin(y=ClusterAllRadialRankings[2], name='Cluster 3', marker_color=clrs[5],
                                points='all', box_visible=True, meanline_visible=True))
        fig2.update_layout(yaxis_title='Radial score', title_text="<b>Cluster Radials: Continuous</b>")
        fig2.show()


        # Demo one plots (OLD)
        # if len(self.ProfileJaccardSimilarityMarix) != 0:
        #     JS_display = plt.figure(num=1, figsize=(6, 6), dpi=150)
        #     JS_plot = JS_display.subplots()
        #     JS_plot.set_title('Jaccard Similarity')
        #     JS_plot.set_aspect('auto')
        #     JS_plot.set_autoscalex_on(True)
        #     plt.imshow(self.ProfileJaccardSimilarityMarix)
        #     plt.colorbar(orientation='vertical')
        #
        # if len(self.ProfileJaccardDifferenceMarix) != 0:
        #     JD_display = plt.figure(num=2, figsize=(6, 6), dpi=150)
        #     JD_plot = JD_display.subplots()
        #     JD_plot.set_title('Jaccard Difference')
        #     JD_plot.set_aspect('auto')
        #     JD_plot.set_autoscalex_on(True)
        #     plt.imshow(self.ProfileJaccardDifferenceMarix)
        #     plt.colorbar(orientation='vertical')
        #
        #     gX_Values = []
        #     gY_Values = []
        #     clusterColors = ["red", "green", "blue", "orange", "purple", "brown", "pink", "gray", "olive", "cyan", "yellow"]
        #     clusterTitles = []
        #     currentCluster = 1
        #     for Cluster in self.CurrentOptimalNPClusters:
        #         clusterTitles.append("Cluster " + str(currentCluster))
        #         new_gx = []
        #         new_gy = []
        #         # Makeup Distance of scatter plots:
        #         # while len(new_gx) < len(gY_Values):
        #         #    new_gx.append(0)
        #         for profileIndex in Cluster:
        #             for windowIndex in self.ProfileDataSet[profileIndex].IndeciesOfValidWindowSamples:
        #                 new_gx.append(windowIndex)
        #                 new_gy.append(profileIndex)
        #         currentCluster = currentCluster + 1
        #         gX_Values.append(new_gx)
        #         gY_Values.append(new_gy)
        #
        #     # Compute the clusters 2d arrays
        #     fig = plt.figure(num=3, figsize=(6, 6), dpi=150)
        #     ax = fig.add_subplot(111)
        #
        #     # Counting values to increment through cluster properties:
        #     plottedClusters = 0
        #     currentColor = 0
        #     clusterAmmount = len(self.CurrentOptimalNPClusters)
        #     colorArrayLen = len(clusterColors)
        #     # While loop to plot clusters:
        #     while plottedClusters < clusterAmmount:
        #         if currentColor >= colorArrayLen:
        #             currentColor = 0
        #         ax.scatter(gX_Values[plottedClusters], gY_Values[plottedClusters], alpha=.8, c=clusterColors[currentColor], edgecolors='none', s=30, label=clusterTitles[plottedClusters])
        #         plottedClusters = plottedClusters + 1
        #         currentColor = currentColor + 1
        #     plt.title('Optimal Cluster Results V1.0')
        #     plt.legend(loc=2)
        # plt.show()
