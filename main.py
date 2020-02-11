from ProfileDataSet import ProfileDataSet as ProfileSet

###############################################################
# Main loop of the application
def main():
    DataSetProfiles = ProfileSet("Valid NP Samples")
    print("Opening File and extracting Data..." + '\n')
    # Data path on desktop:
    # C:/Users/danny/OneDrive/Documents/SchoolWork/DataAnalasysAndGenomeResearch/GSE64881_segmentation_at_30000bp.passqc.multibam.txt
    # Data path on surface pro:
    # C:/Users/Daniel Ingram/Documents/Schoolwork/GenomeResearch/DataSets/GSE64881_segmentation_at_30000bp.passqc.multibam.txt
    # with open('C:/Users/Daniel Ingram/Documents/Schoolwork/GenomeResearch/DataSets/GSE64881_segmentation_at_30000bp.passqc.multibam.txt') as f:
    with open('./DataSets/GSE64881_segmentation_at_30000bp.passqc.multibam.txt') as f:
        lineCount = 0
        windowCount = 0
        for line in f.readlines():
            # Handles the first Line being the profiles description (4628 Bytes or 0.004628 Mb)
            if lineCount == 0:
                # In the newest revision, the whole first line is passed to the dataset object
                DataSetProfiles.DivideLineAndCreateProfiles(line)
            # In the Event it is not the starting line (Begins the genome window intake) (1656 Bytes or 0.001656 Mb) @ line
            else:
                DataSetProfiles.AddWindowByLine(line, windowCount)
                windowCount = windowCount + 1
            lineCount = lineCount + 1
            # print("Reading Line: " + str(lineCount))
        print("Completed Data Extraction of file. Computing Data Summaries...")
        DataSetProfiles.CalculateSummaryValues()
        DataSetProfiles.printProfileDataSetSummaryInTerminal()


if __name__ == "__main__":
    main()
