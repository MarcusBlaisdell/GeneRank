###

import random
import numpy as np

fileSize = {"0157":500000, "2016C":500000, "CH611":500000, "Co6114":500000, "ED1a":500000, "EDL933-1":500000, "FAP1":500000, "_isolate102":500000, "RS76":500000, "UMN026":500000}
#fileSize = {"0157":39530625, "2016C":37059519, "CH611":34837752, "Co6114":35684248, "ED1a":38008400, "EDL933-1":38957663, "FAP1":35662001, "_isolate102":33888563, "RS76":34293140, "UMN026":38092060}


class rClass():
    weight = {}     # weights are stored as a dict to save space
    target = 0.5    # the target precision
    tau = 50        # tau value for top tau
    lam = 1e-3  # lam for regularization coefficient

    def __init__ (self):
        pass

    ### training function: doTrain:

    def doTrain (self, gene, fileName, sampleSize):
        randomList = []
        randomList = self.getRandomList(gene, sampleSize)

        ### sort the random list in numerical order
        ### to read lines from the file in order:

        randomList.sort()

        ### read the randomly selected lines from the file and use for training:

        ### The format of trainList will be:
        ### [score, [(index0, value), (index1, value), ...(indexn, value)], label]
        trainList = []

        readFile = open(fileName, 'r')

        j = 0
        for lineNum in randomList:

            ### don't process for duplicates:
            if lineNum != j:
                theScore = 0

                for k in range (lineNum - j - 1):
                    ### increment to the line we need to read
                    readString = readFile.readline ()
                j = lineNum

                readString = readFile.readline ()

                ### parse the input string:
                theDataPre = readString.split(" ")[:-1]
                theData = []

                for element in theDataPre:
                    theDataOne = int(element.split(":")[0])
                    theDataTwo = int(element.split(":")[-1])
                    theDataTuple = (theDataOne, theDataTwo)
                    theData.append(theDataTuple)

                try:
                    labelVal = int(readString.split (" ")[-1])
                except:
                    labelVal = 0

                theLabel = labelVal

                trainList.append([theScore, theData, theLabel])

        readFile.close ()

        ### Evaluate the lines using current weight:

        ### store the scores of all kmers with a positive label:
        positive_scores = []

        '''
        for kmer in trainList:
            kmer[0] = self.scoreKmer (kmer)
            if kmer[2] == 1:
                positive_scores.append(kmer[0])
        '''

        for kmer in trainList:
            kmer[0] = self.l2regularization (kmer)
            if kmer[2] == 1:
                positive_scores.append(kmer[0])


        #trainList.sort ()

        ### update if there is a positive score below the threshold

        positive_scores = np.array(positive_scores)
        threshold = np.percentile(positive_scores, self.target * 100)

        for kmer in trainList:
            if kmer[2] == 1:
                if kmer[0] <= threshold:
                    self.updateWeight (kmer)
            else:
                if kmer[0] > threshold:
                    self.updateWeight (kmer)

    ### end doTrain function

    ### doTest function
    def doTest (self, gene, fileName, writeFile):
        ### store the scores and labels in a list:
        scoreList = []
        positive_scores = []
        tpfb = 0    # true positives
        fpfb = 0    # false positives
        yPlus = 0   # count of actual positives
        cTp = 0     # classification true positive
        cFp = 0     # classification false positive
        cTn = 0     # classification true negative
        cFn = 0     # classification false negative
        cTotal = 0  # classification - total number of records

        readFile = open (fileName, 'r')
        readString = readFile.readline ()
        ### skip the first line
        readString = readFile.readline ()

        while readString:
            ### parse the input string:
            theDataPre = readString.split(" ")[:-1]
            theData = []

            for element in theDataPre:
                theDataOne = int(element.split(":")[0])
                theDataTwo = int(element.split(":")[-1])
                theDataTuple = (theDataOne, theDataTwo)
                theData.append(theDataTuple)

            try:
                labelVal = int(readString.split (" ")[-1])
            except:
                labelVal = 0

            theLabel = labelVal

                ### evaluate classification to contrast with ranking:

            thePrediction = self.classifyKmer ([0, theData, theLabel])
            cTotal += 1

            ### if the prediction is good,
            ### if the label is good, increment true positive
            ### else, increment false positive
            if thePrediction == 1:
                if theLabel == 1:
                    cTp += 1
                else:
                    cFp += 1
            ### but, if the prediction is negative,
            ### if the label is negative, increment true negative
            ### else, increment false negative
            else:
                if theLabel == -1:
                    cTn += 1
                else:
                    cFn += 1

            theScore = self.scoreKmer ([0, theData, theLabel])

            scoreList.append([theScore, theLabel])
            if theLabel == 1:
                positive_scores.append(theScore)
                yPlus += 1

            readString = readFile.readline ()

        readFile.close ()

        scoreList.sort (reverse=True)

        positive_scores = np.array (positive_scores)
        threshold = np.percentile(positive_scores, self.target * 100)
        print '(int(len(scoreList) * (float (self.tau) / 100))): ', (int(len(scoreList) * (float (self.tau) / 100)))
        print 'len(scoreList): ', len(scoreList)
        print 'float (self.tau): ', float(self.tau)
        print 'float (self.tau) / 100', float (self.tau) / 100

        for i in range (int(len(scoreList) * (float (self.tau) / 100))):
            if scoreList[i][0] > threshold:
                if scoreList[i][1] == 1:
                    tpfb += 1
                else:
                    fpfb += 1

        '''
        print 'tpfb: ', tpfb
        print 'fpfb: ', fpfb
        print 'yPlus: ',
        '''
        '''
        writeFile.write ( 'tpfb: ' + str(tpfb)  + '\n')
        writeFile.write ( 'fpfb: ' + str (fpfb)  + '\n')
        writeFile.write ( 'yPlus: ' + str (yPlus)  + '\n')
        '''
        if (tpfb + fpfb) > 0:
            rPrecision = float (tpfb) / (tpfb + fpfb)
        else:
            rPrecision = 0
        if (yPlus) > 0:
            rRecall = float(tpfb) / (yPlus)
        else:
            rRecall = 0
        if (rPrecision + rRecall) > 0:
            rF1 = 2 * ( (rPrecision * rRecall) / (rPrecision + rRecall) )
        else:
            rF1 = 0

        if (cTotal > 0):
            cAccuracy = float (cTp) / cTotal
        else:
            cAccuracy = 0
        if (cTp + cFp) > 0:
            cPrecision = float (cTp) / (cTp + cFp)
        else:
            cPrecision = 0
        if (cTp + cFn) > 0:
            cRecall = float(cTp) / (cTp + cFn)
        else:
            cRecall = 0
        if (cPrecision + cRecall) > 0:
            cF1 = 2 * ( (cPrecision * cRecall) / (cPrecision + cRecall) )
        else:
            cF1 = 0

        writeFile.write (', ' + str(tpfb) + ', ' + str (fpfb) + ', ' + str (yPlus) )
        writeFile.write (', ' + str(rPrecision) + ', ' + str(rRecall) + ', ' + str(rF1) )
        writeFile.write (', ' + str(cTp) + ', ' + str(cFp) + ', ' + str(cTn) + ', ' + str(cFn) )
        writeFile.write (', ' + str(cTotal) + ', ')
        writeFile.write (', ' + str(cAccuracy) + ', ' + str(cPrecision) + ', ' + str(cRecall) )
        writeFile.write (', ' + str(cF1) + '\n')

    ### end doTest function

    ### updateWeight function
    def updateWeight (self, kmer):
        for index in kmer[1]:
            if self.weight.get(index[0], '--') == '--':
                self.weight[index[0]] = index[1]
            else:
                self.weight[index[0]] += index[1] * kmer[2]
    ### end updateWeight function

    ### l1regularization function
    def l1regularization (self, kmer):
        theScore = 0.0
        Y = kmer[2]
        l1SumOne = 0.0
        l1SumTwo = 0.0

        for pair in kmer[1]:
            if self.weight.get (pair[0], '--') == '--':
                self.weight[pair[0]] = 0
            else:
                l1SumOne += pair[1]
                l1SumTwo += abs(self.weight.get(pair[0]))
        l1SumTwo = l1SumTwo * self.lam

        theScore = ((Y - l1SumOne)**2) + l1SumTwo

        return theScore

    ### end l1regularization function

    ### l2regularization function
    def l2regularization (self, kmer):
        theScore = 0.0
        Y = kmer[2]
        l1SumOne = 0.0
        l1SumTwo = 0.0

        for pair in kmer[1]:
            if self.weight.get (pair[0], '--') == '--':
                self.weight[pair[0]] = 0
            else:
                l1SumOne += pair[1]
                l1SumTwo += ((self.weight.get(pair[0]))**2)
        l1SumTwo = l1SumTwo * self.lam

        theScore = ((Y - l1SumOne)**2) + l1SumTwo

        return theScore

    ### end l2regularization function

    ### scoreKmer function

    def scoreKmer (self, kmer):
        theScore = 0.0
        for index in kmer[1]:
            if self.weight.get (index[0], '--') != '--':
                theScore += self.weight.get (index[0]) * index[1]
            else:
                self.weight[index[0]] = 0
        return theScore

    ### end scoreKmer function

    ### classifyKmer function

    def classifyKmer (self, kmer):
        theScore = 0.0
        thePrediction = 0
        theLabel = kmer[2]

            ### index[0] is the index in the weight dictionary (vector)
            ### index[1] is the weight value at that index in the sample
            ### sample label is kmer[2]
            ### first, calculate dot product:
        for index in kmer[1]:
            if self.weight.get (index[0], '--') != '--':
                theScore += self.weight.get (index[0]) * index[1]
            else:
                self.weight[index[0]] = 0

        ### multiple by the value of the label to make a prediction:
        theScore = theScore * theLabel

        if theScore >= 0:
            thePrediction = 1
        else:
            thePrediction = -1

        return thePrediction

    ### end classifyKmer function

    ### getRandomList function

    def getRandomList (self, gene, sampleSize):
        ### store random numbers in a local list that we will return
        returnRandomList = []

        ### get a list of random line numbers in the gene file

        for i in range(sampleSize):
            returnRandomList.append(random.randrange(0, fileSize.get(gene)))

        return returnRandomList

    ### end getRandomList function
