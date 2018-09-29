###

from rClass import *
import time

myRClass = rClass ()
T = 3   # Number of iterations
n = 11  # number of N-grams to use
sampleSize = 100000

#genePath = "/data/doppa/users/mblaisdell/genome-assembly-datasets/"
genePath = '/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/'
#geneList = ['0157', '2016C', 'CH611', 'Co6114']
geneList = ['2016C', 'CH611', 'Co6114', 'ED1a', 'EDL933-1', 'FAP1', '_isolate102', 'RS76', 'UMN026']

geneTest = '0157'

#outFile = 'outfile-l2_1e-3.csv'
outFile = '/Users/MarcusBlaisdell/Documents/REU/RankRedo/outfile-l2_1e-3_50-test.csv'
writeFile = open (outFile, 'w')

writeFile.write ('notes, rank-Tp, rank-Fp, rank-yPlus, rank-Precision, rank-Recall, rank-F1, ')
writeFile.write ('class-Tp, class-Fp, class-Tn, class-Fn, class-Total, class-Accuracy, ' )
writeFile.write ('class-Precision, class-Recall, class-F1 \n')

totalTimeStart = time.time()
for t in range(T):
    writeFile.write ('iteration: ' + str(t)  + '\n')

    trainStart = time.time ()

    for gene in geneList:
        fileName = genePath + gene + '-' + str(n)

        myRClass.doTrain (gene, fileName, sampleSize)

    trainEnd = time.time ()

    writeFile.write ('trainTime: ' + str (trainEnd - trainStart)  + '\n')

    testStart = time.time ()

    myRClass.doTest (geneTest, fileName, writeFile)

    testEnd = time.time ()

    writeFile.write ('testTime: ' + str (testEnd - testStart) + '\n' )

totalTimeEnd = time.time()

writeFile.write ('total Time: ' + str(totalTimeEnd - totalTimeStart) + '\n')

writeFile.close ()
