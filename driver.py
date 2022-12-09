import proteinIndex
import sys

if len(sys.argv) == 2: #Build by loading from file
    programName = sys.argv[0]
    inputFile = str(sys.argv[1]) # input file, an encoded FM index
    pIndex = proteinIndex.FMIndex(readFile=inputFile)
elif len(sys.argv) == 4: #Build from list of IDs
    programName = sys.argv[0]
    inputFile = str(sys.argv[1]) # input file, a list of protein IDs
    cpStep = int(sys.argv[2])
    saStep = int(sys.argv[3])
    pIndex = proteinIndex.buildProteinIndex(inputFile, cpStep, saStep)
elif len(sys.argv) == 5: #Same as previous but save to file
    programName = sys.argv[0]
    inputFile = str(sys.argv[1]) # input file, a list of protein IDs
    cpStep = int(sys.argv[2])
    saStep = int(sys.argv[3])
    saveFile = str(sys.argv[4])
    pIndex = proteinIndex.buildProteinIndex(inputFile, cpStep, saStep)
    proteinIndex.saveProteinIndex(saveFile, pIndex)
else:
    print("Enter valid args: \na) inputFile (encoding a previously saved FM index) \nb) inputFile (listing protein IDs), cpStep, saStep, saveFile")
    #TESTING
    pIndex = proteinIndex.FMIndex("MTVSSHRLELLSPARDAAIAREAILHGADAVYIGGPGFGARHNASNSLKDIAELVPFAHRYGAKIFVTLNTILHDDELEPAQRLITDLYQTGVDALIVQDMGILELDIPPIELHASTQCDIRTVEKAKFLSDVGFTQIVLARELNLDQIRAIHQATDATIEFFIHGALCVAYSGQCYISHAQTGRSANRGDCSQACRLPYTLKDDQGRVVSYEKHLLSMKDNDQTANLGALIDAGVRSFKIEGRYKDMSYVKNITAHYRQMLDAIIEERGDLARASSGRTEHFFVPSTEKTFHRGSTDYFVNARKGDIGAFDSPKFIGLPVGEVVKVAKDHLDVAVTEPLANGDGLNVLIKREVVGFRANTVEKTGENQYRVWPNEMPADLHKIRPHHPLNRNLDHNWQQALTKTSSERRVAVDIELGGWQEQLILTLTSEEGVSITHTLDGQFDEANNAEKAMNNLKDGLAKLGQTLYYARDVQINLPGALFVPNSLLNQFRREAADMLDAARLASYQRGSRKPVADPAPVYPQTHLSFLANVYNQKAREFYHRYGVQLIDAAYEAHEEKGEVPVMITKHCLRFAFNLCPKQAKGNIKSWKATPMQLVNGDEVLTLKFDCRPCEMHVIGKIKNHILKMPLPGSVVASVSPDELLKTLPKRKG$", 7, 20)
    print(pIndex.query('QIVLARELNLDQIRAIHQATDATIEF'))
    proteinIndex.saveProteinIndex("index.txt", pIndex)
