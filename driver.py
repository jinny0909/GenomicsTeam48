import proteinIndex
import sys

if len(sys.argv) == 2: #load from file
    programName = sys.argv[0]
    inputFile = str(sys.argv[1]) # input file, an encoded FM index
    pIndex = proteinIndex.FMIndex(inputFile)
if len(sys.argv) == 5:
    programName = sys.argv[0]
    inputFile = str(sys.argv[1]) # input file, a list of protein IDs
    cpStep = int(sys.argv[2])
    saStep = int(sys.argv[3])
    saveFile = str(sys.argv[4])
    pIndex = proteinIndex.buildProteinIndex(inputFile, cpStep, saStep)
else:
    #print("Enter valid args: \na) inputFile (encoding a previously saved FM index) \nb) inputFile (listing protein IDs), cpStep, saStep, saveFile")
    #TESTING
    pIndex = proteinIndex.FMIndex("abcdefghijdefklmnopqdefrstuv$", 1, 1)

print(pIndex.bwt_t)
print(pIndex.query("def"))