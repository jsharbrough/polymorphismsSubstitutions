import sys
# -*- coding: utf-8 -*-
def polymorphismsSubstitutionsHaploid(alignment):
    speciesDict = {'MH':[">NC_012920.1 Homo sapiens mitochondrion, complete genome",	">AF346963.1 Homo sapiens mitochondrion, complete genome",	">AF346964.1 Homo sapiens mitochondrion, complete genome",	">AF346965.1 Homo sapiens mitochondrion, complete genome",	">AF346966.1 Homo sapiens mitochondrion, complete genome",	">AF346967.1 Homo sapiens mitochondrion, complete genome",	">AF346968.1 Homo sapiens mitochondrion, complete genome",	">AF346969.1 Homo sapiens mitochondrion, complete genome",	">AF346970.1 Homo sapiens mitochondrion, complete genome",	">AF346971.1 Homo sapiens mitochondrion, complete genome",	">AF346972.1 Homo sapiens mitochondrion, complete genome",	">AF346973.1 Homo sapiens mitochondrion, complete genome",	">AF346974.1 Homo sapiens mitochondrion, complete genome",	">AF346975.1 Homo sapiens mitochondrion, complete genome",	">AF346976.1 Homo sapiens mitochondrion, complete genome",	">AF346977.1 Homo sapiens mitochondrion, complete genome",	">AF346978.1 Homo sapiens mitochondrion, complete genome",	">AF346979.1 Homo sapiens mitochondrion, complete genome",	">AF346980.1 Homo sapiens mitochondrion, complete genome",	">AF346981.1 Homo sapiens mitochondrion, complete genome",	">AF346982.1 Homo sapiens mitochondrion, complete genome",	">AF346983.1 Homo sapiens mitochondrion, complete genome",	">AF346984.1 Homo sapiens mitochondrion, complete genome",	">AF346985.1 Homo sapiens mitochondrion, complete genome",	">AF346986.1 Homo sapiens mitochondrion, complete genome",	">AF346987.1 Homo sapiens mitochondrion, complete genome",	">AF346988.1 Homo sapiens mitochondrion, complete genome",	">AF346989.1 Homo sapiens mitochondrion, complete genome",	">AF346990.1 Homo sapiens mitochondrion, complete genome",	">AF346991.1 Homo sapiens mitochondrion, complete genome",	">AF346992.1 Homo sapiens mitochondrion, complete genome",	">AF346993.1 Homo sapiens mitochondrion, complete genome",	">AF346994.1 Homo sapiens mitochondrion, complete genome",	">AF346995.1 Homo sapiens mitochondrion, complete genome",	">AF346996.1 Homo sapiens mitochondrion, complete genome",	">AF346997.1 Homo sapiens mitochondrion, complete genome",	">AF346998.1 Homo sapiens mitochondrion, complete genome",	">AF346999.1 Homo sapiens mitochondrion, complete genome",	">AF347000.1 Homo sapiens mitochondrion, complete genome",	">AF347001.1 Homo sapiens mitochondrion, complete genome",	">AF347002.1 Homo sapiens mitochondrion, complete genome",	">AF347003.1 Homo sapiens mitochondrion, complete genome",	">AF347004.1 Homo sapiens mitochondrion, complete genome",	">AF347005.1 Homo sapiens mitochondrion, complete genome",	">AF347006.1 Homo sapiens mitochondrion, complete genome",	">AF347007.1 Homo sapiens mitochondrion, complete genome",	">AF347008.1 Homo sapiens mitochondrion, complete genome",	">AF347009.1 Homo sapiens mitochondrion, complete genome",	">AF347010.1 Homo sapiens mitochondrion, complete genome",	">AF347011.1 Homo sapiens mitochondrion, complete genome",	">AF347012.1 Homo sapiens mitochondrion, complete genome",	">AF347013.1 Homo sapiens mitochondrion, complete genome",	">AF347014.1 Homo sapiens mitochondrion, complete genome",	">AF347015.1 Homo sapiens mitochondrion, complete genome"],'N':[">FM865407.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 1",	">FM865408.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 2",	">FM865409.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate El Sidron 1253",	">FM865410.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Vindija 33.25",	">FM865411.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Mezmaiskaya 1",	">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome"],'D':[">FR695060.1 Homo sp. Altai complete mitochondrial genome, isolate Denisova molar",	">KT780370.1 Homo sapiens ssp. Denisova isolate Denisova8 mitochondrion, complete genome",	">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia"],'PoB':[">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"]}
    referenceDict = {'MH':">NC_012920.1 Homo sapiens mitochondrion, complete genome",'N':">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome",'D':">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia",'PoB':">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"}
    speciesList = ['MH','N','D','PoB']
    geneticCode = {"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"} #standard code
    startCodons = ['TTG','CTG','ATG']
    regions = {"D-Loop":{"type":"noncoding","region":"16049-592","direction":"forward"},	"TRNF":{"type":"tRNA","region":"593-663","direction":"forward"},	"RNR1":{"type":"rRNA","region":"664-1617","direction":"forward"},	"TRNV":{"type":"tRNA","region":"1618-1686","direction":"forward"},	"RNR2":{"type":"rRNA","region":"1687-3247","direction":"forward"},	"TRNL1":{"type":"tRNA","region":"3248-3322","direction":"forward"},	"ND1":{"type":"protein","region":"3325-4280","direction":"forward"},	"TRNI":{"type":"tRNA","region":"4281-4349","direction":"forward"},	"TRNQ":{"type":"tRNA","region":"4347-4418","direction":"complement"},	"TRNM":{"type":"tRNA","region":"4420-4487","direction":"forward"},	"ND2":{"type":"protein","region":"4488-5531","direction":"forward"},	"TRNW":{"type":"tRNA","region":"5530-5597","direction":"forward"},	"TRNA":{"type":"tRNA","region":"5605-5673","direction":"complement"},	"TRNN":{"type":"tRNA","region":"5675-5747","direction":"complement"},	"TRNC":{"type":"tRNA","region":"5779-5845","direction":"complement"},	"TRNY":{"type":"tRNA","region":"5845-5910","direction":"complement"},	"COI":{"type":"protein","region":"5929-7470","direction":"forward"},	"TRNS1":{"type":"tRNA","region":"7471-7539","direction":"complement"},	"TRND":{"type":"tRNA","region":"7543-7610","direction":"forward"},	"COII":{"type":"protein","region":"7611-8294","direction":"forward"},	"TRNK":{"type":"tRNA","region":"8320-8389","direction":"forward"},	"ATP8":{"type":"protein","region":"8391-8597","direction":"forward"},	"ATP6":{"type":"protein","region":"8552-9232","direction":"forward"},	"COIII":{"type":"protein","region":"9232-10015","direction":"forward"},	"TRNG":{"type":"tRNA","region":"10016-10083","direction":"forward"},	"ND3":{"type":"protein","region":"10084-10429","direction":"forward"},	"TRNR":{"type":"tRNA","region":"10430-10494","direction":"forward"},	"ND4L":{"type":"protein","region":"10495-10791","direction":"forward"},	"ND4":{"type":"protein","region":"10785-12162","direction":"forward"},	"TRNH":{"type":"tRNA","region":"12163-12231","direction":"forward"},	"TRNS2":{"type":"tRNA","region":"12232-12290","direction":"forward"},	"TRNL2":{"type":"tRNA","region":"12291-12361","direction":"forward"},	"ND5":{"type":"protein","region":"12362-14173","direction":"forward"},	"ND6":{"type":"protein","region":"14174-14698","direction":"complement"},	"TRNE":{"type":"tRNA","region":"14699-14767","direction":"complement"},	"CYTB":{"type":"protein","region":"14772-15912","direction":"forward"},	"TRNT":{"type":"tRNA","region":"15913-15978","direction":"forward"},	"TRNP":{"type":"tRNA","region":"15981-16048","direction":"complement"}, 'whole mt genome':{"type":"noncoding","region":"1-16599","direction":"forward"}}
    regionList = ["D-Loop",	"TRNF",	"RNR1",	"TRNV",	"RNR2",	"TRNL1",	"ND1",	"TRNI",	"TRNQ",	"TRNM",	"ND2",	"TRNW",	"TRNA",	"TRNN",	"TRNC",	"TRNY",	"COI",	"TRNS1",	"TRND",	"COII",	"TRNK",	"ATP8",	"ATP6",	"COIII",	"TRNG",	"ND3",	"TRNR",	"ND4L",	"ND4",	"TRNH",	"TRNS2",	"TRNL2",	"ND5",	"ND6",	"TRNE",	"CYTB",	"TRNT",	"TRNP",	"whole mt genome"]
    seqDict = {}
    seqList = []
    currSeq = ''
    seqName = ''
    alignFile = open(alignment,'r')
    for line in alignFile:
        if line[0] == '>':
            if seqName != '':
                seqDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
    	    seqList.append(seqName)
    	    currSeq = ''
        else:	
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    seqDict[seqName] = currSeq
    alignFile.close()
    featureLengths = [1141]
    featurePolyDict = {}
    for feature in regionList:
        logfile = open('logfile_test.txt','a')
        logfile.write('~~~~~~~~~~~' + feature + ' – Polymorphisms ~~~~~~~~~~~\n')
        logfile.close()
        items = regions[feature]
        direction = items['direction']
        region = items['region']
        regionSplit = region.split('-')
        regionStart = int(regionSplit[0]) - 1
        regionStop = int(regionSplit[1]) - 1
        currSeqDict = {}
        if feature == 'D-Loop':
            for seq in seqList:
                currSeq = seqDict[seq]
                regionSeq = currSeq[regionStart:] + currSeq[0:regionStop]
                currSeqDict[seq] = regionSeq
        else:    
            featureLengths.append(regionStop - regionStart + 1)
            for seq in seqList:
                if direction == 'forward':
                    currSeq = seqDict[seq]
                    regionSeq = currSeq[regionStart:regionStop + 1]
                    currSeqDict[seq] = regionSeq
                else:
                    currSeq = seqDict[seq]
                    regionSeq = currSeq[regionStart:regionStop + 1]
                    newSeq = reverseComplementHaploid(regionSeq)
                    currSeqDict[seq] = newSeq
        regionType = items['type']
        speciesNum = 0
        polyDict = {}
        for species in speciesList:
            currSpeciesList = speciesDict[species]
            if regionType == 'protein':
                codonDict,AADict = buildCodonDictHaploid(currSeqDict)
                outfile = open(feature + '_aaSeq.fasta','w')
                for seq in seqList:
                    outfile.write(seq + '\n')
                    outfile.write(AADict[seq] + '\n')
                outfile.close()
                if len(currSpeciesList) > 1:
                    invariantSites,nonsynPolymorphisms,synPolymorphisms = polymorphicCodonsHaploid(currSeqDict,currSpeciesList)
                    totalPolymorphisms = nonsynPolymorphisms + synPolymorphisms
                    polyDict[species] = [invariantSites, nonsynPolymorphisms,synPolymorphisms, totalPolymorphisms]
                else:
                    invariantSites = range(len(seqDict[currSpeciesList[0]]))
                    polyDict[species] = [invariantSites]
            else:      
                if len(currSpeciesList) > 1:
                    invariantSites,totalPolymorphisms = polymorphicSitesHaploid(currSeqDict,currSpeciesList)
                    polyDict[species] = [invariantSites,totalPolymorphisms]
                else:
                    invariantSites = range(len(seqDict[currSpeciesList[0]]))
                    polyDict[species] = [invariantSites]
        featurePolyDict[feature] = polyDict
    featureSubDict = {}
    for feature in regionList:
        logfile = open('logfile_test.txt','a')
        logfile.write('~~~~~~~~~~~' + feature + ' – Fixed Differences ~~~~~~~~~~~\n')
        logfile.close()
        polyDict = featurePolyDict[feature]
        items = regions[feature]
        direction = items['direction']
        region = items['region']
        regionSplit = region.split('-')
        regionStart = int(regionSplit[0]) - 1
        regionStop = int(regionSplit[1]) - 1
        currSeqDict = {}
        if feature == 'D-Loop':
            for seq in seqList:
                currSeq = seqDict[seq]
                regionSeq = currSeq[regionStart:] + currSeq[0:regionStop]
                currSeqDict[seq] = regionSeq
        else:    
            featureLengths.append(regionStop - regionStart)
            for seq in seqList:
                if direction == 'forward':
                    currSeq = seqDict[seq]
                    regionSeq = currSeq[regionStart:regionStop +1]
                    currSeqDict[seq] = regionSeq
                else:
                    currSeq = seqDict[seq]
                    regionSeq = currSeq[regionStart:regionStop + 1]
                    newSeq = reverseComplementHaploid(regionSeq)
                    currSeqDict[seq] = newSeq
        speciesNum = 1
        subDict = {}
        for speciesA in speciesList:
            speciesAPolyDict = polyDict[speciesA]
            invariantSitesA = speciesAPolyDict[0]
            for speciesB in speciesList[speciesNum:]:
                comparison = speciesA + '-' + speciesB
                speciesBPolyDict = polyDict[speciesB]
                invariantSitesB = speciesBPolyDict[0]
                if items['type'] == 'protein':
                    nonsynFixedDifferences,synFixedDifferences = fixedDifferencesCodonsHaploid(currSeqDict,invariantSitesA,invariantSitesB,speciesA,speciesB)
                    totalFixedDifferences = synFixedDifferences + nonsynFixedDifferences
                    subDict[comparison] = [totalFixedDifferences,nonsynFixedDifferences,synFixedDifferences]
                else:
                    totalFixedDifferences = fixedDifferencesSitesHaploid(currSeqDict,invariantSitesA,invariantSitesB,speciesA,speciesB)
                    subDict[comparison] = [totalFixedDifferences]   
            speciesNum += 1
        featureSubDict[feature] = subDict
    outfile = open('test.out','w')
    outfile.write('\t\t\tMH-N\t\t\tMH-D\t\t\tMH-SdlH\t\t\tN-D\t\t\tN-SdlH\t\t\tD-SdlH\t\t\tIntraspecific MH\t\t\tIntraspecific N\t\t\tIntraspecfici D\t\t\nRegion\tLength\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Polymorphisms\tSynonymous Polymorphisms\tTotal Polymorphisms\tNonsynonymous Polymorphisms\tSynonymous Polymorphisms\tTotal Polymorphisms\tNonsynonymous Polymorphisms\tSynonymous Polymorphisms\tTotal Polymorphisms\n')
    geneNum = 0
    for feature in regionList:
        outfile.write(feature + '\t')
        outfile.write(str(featureLengths[geneNum]))
        items = regions[feature]
        polyDict = featurePolyDict[feature]
        subDict = featureSubDict[feature]
        regionType = items['type']
        speciesNum = 1
        for speciesA in speciesList:
            for speciesB in speciesList[speciesNum:]:
                comparison = speciesA + '-' + speciesB
                compNumbers = subDict[comparison]
                if regionType == 'protein':
                    outfile.write('\t' + str(compNumbers[1]) + '\t' + str(compNumbers[2]) + '\t' + str(compNumbers[0]))
                else:
                    outfile.write('\t-\t-\t' + str(compNumbers[0]))
            speciesNum += 1
        for species in speciesList:
            if len(speciesDict[species]) > 1:
                polNumbers = polyDict[species]
                if regionType == 'protein':
                    outfile.write('\t' + str(polNumbers[1]) + '\t' + str(polNumbers[2]) + '\t' + str(polNumbers[3]))
                else:
                    outfile.write('\t-\t-\t' + str(polNumbers[1]))
        outfile.write('\n')
        geneNum += 1

    outfile.close()
        

                 
def polymorphicSitesHaploid(seqDict,seqList):
    speciesDict = {">NC_012920.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346963.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346964.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346965.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346966.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346967.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346968.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346969.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346970.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346971.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346972.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346973.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346974.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346975.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346976.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346977.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346978.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346979.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346980.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346981.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346982.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346983.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346984.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346985.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346986.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346987.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346988.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346989.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346990.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346991.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346992.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346993.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346994.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346995.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346996.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346997.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346998.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346999.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347000.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347001.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347002.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347003.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347004.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347005.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347006.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347007.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347008.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347009.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347010.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347011.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347012.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347013.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347014.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347015.1 Homo sapiens mitochondrion, complete genome":"MH",	">FM865407.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 1":"N",	">FM865408.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 2":"N",	">FM865409.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate El Sidron 1253":"N",	">FM865410.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Vindija 33.25":"N",	">FM865411.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Mezmaiskaya 1":"N",	">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome":"N",	">FR695060.1 Homo sp. Altai complete mitochondrial genome, isolate Denisova molar":"D",	">KT780370.1 Homo sapiens ssp. Denisova isolate Denisova8 mitochondrion, complete genome":"D",	">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia":"D",	">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome":"PoB"}
    regionDict = {"D-Loop":(1,592),"TRNF":(593,663),"RNR1":(664,1617),"TRNV":(1618,1686),"RNR2":(1687,3247),"TRNL1":(3248,3322),"ND1":(3325,4280),"TRNI":(4281,4349),"TRNQ":(4347,4418),"TRNM":(4420,4487),"ND2":(4488,5531),"TRNW":(5530,5597),"TRNA":(5605,5673),"TRNN":(5675,5747),"TRNC":(5779,5845),"TRNY":(5845,5910),"COI":(5929,7470),"TRNS1":(7471,7539),"TRND":(7543,7610),"COII":(7611,8294),"TRNK":(8320,8389),"ATP8":(8391,8597),"ATP6":(8552,9232),"COIII":(9232,10015),"TRNG":(10016,10083),"ND3":(10084,10429),"TRNR":(10430,10494),"ND4L":(10495,10791),"ND4":(10785,12162),"TRNH":(12163,12231),"TRNS2":(12232,12290),"TRNL2":(12291,12361),"ND5":(12362,14173),"ND6":(14174,14698),"TRNE":(14699,14767),"CYTB":(14772,15912),"TRNT":(15913,15978),"TRNP":(15981,16048), "D-Loop2":(16048,16599)}
    regionList = ["D-Loop",	"TRNF",	"RNR1",	"TRNV",	"RNR2",	"TRNL1",	"ND1",	"TRNI",	"TRNQ",	"TRNM",	"ND2",	"TRNW",	"TRNA",	"TRNN",	"TRNC",	"TRNY",	"COI",	"TRNS1",	"TRND",	"COII",	"TRNK",	"ATP8",	"ATP6",	"COIII",	"TRNG",	"ND3",	"TRNR",	"ND4L",	"ND4",	"TRNH",	"TRNS2",	"TRNL2",	"ND5",	"ND6",	"TRNE",	"CYTB",	"TRNT",	"TRNP","D-Loop2"]
    seqLength = len(seqDict[seqList[0]])
    currSpecies = speciesDict[seqList[0]]
    i = 0
    numPolymorphisms = 0
    invariantSites = []
    polymorphicSites = []
    while i < seqLength:
        currSite = []
        for seq in seqList:
            currSeq = seqDict[seq]
            if currSeq[i] not in currSite:
                currSite.append(currSeq[i])
        if len(currSite) > 1 and '-' not in currSite and 'N' not in currSite:
            numPolymorphisms += len(currSite) - 1
            polymorphicSites.append(i)
            logfile = open('logfile_test.txt','a')
            logfile.write(currSpecies + ': Site Num – ' + str(i + 1) + ' (' + currSite[0])
            for allele in currSite[1:]:
                logfile.write(', ' + allele)
            logfile.write(') -- Polymorphism -- ')
            regions = []
            for region in regionList:
                regionPositions = regionDict[region]
                if region == 'TRNQ' or region == 'TRNA' or region == 'TRNN' or region == 'TRNC' or region == 'TRNY' or region == 'TRNS1' or region == 'ND6' or region == 'TRNE' or region == 'TRNP':
                    regionStart = regionPositions[0]
                    regionStop = regionPositions[1] + 1
                else:
                    regionStart = regionPositions[0] - 1
                    regionStop = regionPositions[1]
                if i >= regionStart and i <= regionStop:
                    regions.append(region)
            if len(regions) > 0:
                logfile.write(str(regions) + '\n')
            else:
                logfile.write('intergenic\n')
            logfile.close() 
        elif len(currSite) == 1 and '-' not in currSite and 'N' not in currSite:
            invariantSites.append(i)
        i += 1
    return invariantSites, numPolymorphisms
          
def fixedDifferencesSitesHaploid(seqDict,invariantSitesA,invariantSitesB,species,outgroup):
    referenceDict = {'MH':">NC_012920.1 Homo sapiens mitochondrion, complete genome",'N':'>AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome','D':">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia",'PoB':">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"}
    regionDict = {"D-Loop":(1,592),"TRNF":(593,663),"RNR1":(664,1617),"TRNV":(1618,1686),"RNR2":(1687,3247),"TRNL1":(3248,3322),"ND1":(3325,4280),"TRNI":(4281,4349),"TRNQ":(4347,4418),"TRNM":(4420,4487),"ND2":(4488,5531),"TRNW":(5530,5597),"TRNA":(5605,5673),"TRNN":(5675,5747),"TRNC":(5779,5845),"TRNY":(5845,5910),"COI":(5929,7470),"TRNS1":(7471,7539),"TRND":(7543,7610),"COII":(7611,8294),"TRNK":(8320,8389),"ATP8":(8391,8597),"ATP6":(8552,9232),"COIII":(9232,10015),"TRNG":(10016,10083),"ND3":(10084,10429),"TRNR":(10430,10494),"ND4L":(10495,10791),"ND4":(10785,12162),"TRNH":(12163,12231),"TRNS2":(12232,12290),"TRNL2":(12291,12361),"ND5":(12362,14173),"ND6":(14174,14698),"TRNE":(14699,14767),"CYTB":(14772,15912),"TRNT":(15913,15978),"TRNP":(15981,16048), "D-Loop2":(16048,16599)}
    regionList = ["D-Loop",	"TRNF",	"RNR1",	"TRNV",	"RNR2",	"TRNL1",	"ND1",	"TRNI",	"TRNQ",	"TRNM",	"ND2",	"TRNW",	"TRNA",	"TRNN",	"TRNC",	"TRNY",	"COI",	"TRNS1",	"TRND",	"COII",	"TRNK",	"ATP8",	"ATP6",	"COIII",	"TRNG",	"ND3",	"TRNR",	"ND4L",	"ND4",	"TRNH",	"TRNS2",	"TRNL2",	"ND5",	"ND6",	"TRNE",	"CYTB",	"TRNT",	"TRNP","D-Loop2"]
    speciesA = referenceDict[species]
    speciesB = referenceDict[outgroup]
    speciesASeq = seqDict[speciesA]
    speciesBSeq = seqDict[speciesB] 
    numFixedDifferences = 0
    if invariantSitesA == False:
        invariantSitesA = range(len(speciesASeq))
    if invariantSitesB == False:
        invariantSitesB = range(len(speciesBSeq))
    invariantSites = []
    for siteA in invariantSitesA:
        for siteB in invariantSitesB:
            if siteA == siteB and siteA not in invariantSites:
                invariantSites.append(siteA)
    for site in invariantSites:
        currSite = [speciesASeq[site],speciesBSeq[site]]
        if currSite[0] != '-' and currSite[0] != 'N' and currSite[1] != '-' and currSite[1] != 'N':
            if currSite[0] != currSite[1]:
                numFixedDifferences += 1  
                comparison = species + '-' + outgroup
                logfile = open('logfile_test.txt','a')
                logfile.write(comparison + ': Site Num – ' + str(site + 1) + ' (' + currSite[0])
                for allele in currSite[1:]:
                    logfile.write(', ' + allele)
                logfile.write(') -- Fixed Difference -- ')
                regions = []
                for region in regionList:
                    regionPositions = regionDict[region]
                    if region == 'TRNQ' or region == 'TRNA' or region == 'TRNN' or region == 'TRNC' or region == 'TRNY' or region == 'TRNS1' or region == 'ND6' or region == 'TRNE' or region == 'TRNP':
                        regionStart = regionPositions[0]
                        regionStop = regionPositions[1] + 1
                    else:
                        regionStart = regionPositions[0] - 1
                        regionStop = regionPositions[1]
                    if site >= regionStart and site <= regionStop:
                        regions.append(region)
                if len(regions) > 0:
                    logfile.write(str(regions) + '\n')
                else:
                    logfile.write('intergenic\n')
                    logfile.close()
    return numFixedDifferences
        
def polymorphicCodonsHaploid(seqDict,seqList):
    speciesDict = {">NC_012920.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346963.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346964.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346965.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346966.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346967.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346968.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346969.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346970.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346971.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346972.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346973.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346974.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346975.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346976.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346977.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346978.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346979.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346980.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346981.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346982.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346983.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346984.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346985.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346986.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346987.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346988.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346989.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346990.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346991.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346992.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346993.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346994.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346995.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346996.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346997.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346998.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF346999.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347000.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347001.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347002.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347003.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347004.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347005.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347006.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347007.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347008.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347009.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347010.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347011.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347012.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347013.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347014.1 Homo sapiens mitochondrion, complete genome":"MH",	">AF347015.1 Homo sapiens mitochondrion, complete genome":"MH",	">FM865407.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 1":"N",	">FM865408.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 2":"N",	">FM865409.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate El Sidron 1253":"N",	">FM865410.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Vindija 33.25":"N",	">FM865411.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Mezmaiskaya 1":"N",	">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome":"N",	">FR695060.1 Homo sp. Altai complete mitochondrial genome, isolate Denisova molar":"D",	">KT780370.1 Homo sapiens ssp. Denisova isolate Denisova8 mitochondrion, complete genome":"D",	">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia":"D",	">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome":"PoB"}
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    startCodons = {'ATT':'M','ATC':'M','ATA':'M','ATG':'M','GTG':'M'}
    codonDict,AADict = buildCodonDictHaploid(seqDict) 
    currSpecies = speciesDict[seqList[0]]
    seqLength = len(seqDict[seqList[0]])
    numCodons = len(codonDict[seqList[0]])
    i = 0
    numPolymorphisms = 0
    numNonsynPolymorphisms = 0
    numSynPolymorphisms = 0
    invariantSites = []
    polymorphicCodons = []
    while i/3 < numCodons:
        currSite = []
        for seq in seqList:
            currSeq = seqDict[seq]
            if currSeq[i] not in currSite:
                currSite.append(currSeq[i])
        if len(currSite) > 1 and '-' not in currSite and 'N' not in currSite:
            numPolymorphisms += len(currSite)
            codonNum = i/3
            if codonNum not in polymorphicCodons:
                polymorphicCodons.append(codonNum)
        elif len(currSite) == 1 and '-' not in currSite and 'N' not in currSite:
            invariantSites.append(i)
        i += 1 
    for codonNumber in polymorphicCodons:
        currCodon = []
        for seq in seqList:
            codons = codonDict[seq]
            if codons[codonNumber] not in currCodon:
                currCodon.append(codons[codonNumber])
        if len(currCodon) > 1:
            poly = True
            j = 0
            while poly == True and j < len(currCodon):
                if '-' in currCodon[j] or 'N' in currCodon[j]:
                    poly = False
                j += 1
            if poly == True:
                AAs = []
                site1 = []
                site2 = []
                site3 = []
                for polyCodon in currCodon:
                    if polyCodon[0] not in site1:
                        site1.append(polyCodon[0])
                    if polyCodon[1] not in site2:
                        site2.append(polyCodon[1])
                    if polyCodon[2] not in site3:
                        site3.append(polyCodon[2])
                    if codonNumber == 0:
                        if polyCodon in startCodons:
                            currAA = startCodons[polyCodon]
                        else:
                            currAA = geneticCode[polyCodon]
                    else:
                        currAA = geneticCode[polyCodon]
                    if currAA not in AAs:
                        AAs.append(currAA)
                totalChanges = (len(site1) - 1) + (len(site2) - 1) + (len(site3) - 1)
                if totalChanges == 1:
                    if len(AAs) == 1:
                        numSynPolymorphisms += 1
                        logfile = open('logfile_test.txt','a')
                        logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                        for allele in currCodon[1:]:
                            logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                        logfile.write(') -- ' + str(totalChanges) + ' Synonymous Polymorphism\n')
                        logfile.close()
                    else:
                        numNonsynPolymorphisms += 1
                        logfile = open('logfile_test.txt','a')
                        logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                        for allele in currCodon[1:]:
                            logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                        logfile.write(') -- ' + str(totalChanges) + ' Nonsynonymous Polymorphism\n')
                        logfile.close()
                elif totalChanges == 2:
                    numSynDiff = 0
                    numNonSynDiff = 0
                    if len(currCodon) == 3:
                        codon1 = currCodon[0]
                        codon2 = currCodon[1]
                        codon3 = currCodon[2]
                        if codonNumber == 0:
                            if codon1 in startCodons:
                                aa1 = startCodons[codon1]
                            else:
                                aa1 = geneticCode[codon1]
                            if codon2 in startCodons:
                                aa2 = startCodons[codon2]
                            else:
                                aa2 = geneticCode[codon2]
                            if codon3 in startCodons:
                                aa3 = startCodons[codon3]
                            else:
                                aa3 = geneticCode[codon3]
                        else:  
                            aa1 = geneticCode[codon1]
                            aa2 = geneticCode[codon2]
                            aa3 = geneticCode[codon3]
                        codon1_2 = 0
                        codon1_3 = 0
                        codon2_3 = 0
                        k = 0
                        while k < 3:
                            if codon1[k] != codon2[k]:
                                codon1_2 += 1
                            if codon1[k] != codon3[k]:
                                codon1_3 += 1
                            if codon2[k] != codon3[k]:
                                codon2_3 += 1
                            k += 1
                        if codon1_2 > codon1_3 and codon1_2 > codon2_3:
                            if aa1 != aa3:
                                numNonSynDiff += 1
                            else:
                                numSynDiff += 1 
                            if aa2 != aa3:
                                numNonSynDiff += 1
                            else:
                                numSynDiff += 1
                        elif codon1_3 > codon1_2 and codon1_3 > codon2_3:
                            if aa2 != aa1:
                                numNonSynDiff += 1
                            else:
                                numSynDiff += 1 
                            if aa2 != aa3:
                                numNonSynDiff += 1
                            else:
                                numSynDiff += 1
                        elif codon2_3 > codon1_3 and codon2_3 > codon1_2:
                            if aa1 != aa2:
                                numNonSynDiff += 1
                            else:
                                numSynDiff += 1 
                            if aa1 != aa3:
                                numNonSynDiff += 1
                            else:
                                numSynDiff += 1
                        else:
                            logfile = open('logfile_test.txt','a')
                            logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' has a complex evolutionary history --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                            for allele in currCodon[1:]:
                                logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                            logfile.write(') -- Polymorphism\n')
                            logfile.close()                            
                    elif len(currCodon) == 2:
                        codon1 = currCodon[0]
                        codon2 = currCodon[1]
                        if codonNumber == 0:
                            if codon1 in startCodons:
                                aa1 = startCodons[codon1]
                            else:
                                aa1 = geneticCode[codon1]
                            if codon2 in startCodons:
                                aa2 = startCodons[codon2]
                            else:
                                aa2 = geneticCode[codon2]
                        else:  
                            aa1 = geneticCode[codon1]
                            aa2 = geneticCode[codon2]
                        sites = [site1,site2,site3]
                        int1 = ''
                        int2 = ''
                        changeNum = 0
                        siteNum = 0
                        for site in sites:
                            if len(site) == 1:
                                int1 += site[0]
                                int2 += site[0]
                            elif changeNum == 0:
                                int1 += codon1[siteNum]
                                int2 += codon2[siteNum]
                                changeNum +=1
                            elif changeNum == 1:
                                int1 += codon2[siteNum]
                                int2 += codon1[siteNum]
                            siteNum += 1
                        int1AA = geneticCode[int1]
                        int2AA = geneticCode[int2]
                        numSynDiff_a = 0
                        numSynDiff_b = 0
                        numNonSynDiff_a = 0
                        numNonSynDiff_b = 0
                        if aa1 != int1AA:
                            numNonSynDiff_a += 1
                            if aa2 != int1AA:
                                numNonSynDiff_a += 1
                            else:
                                numSynDiff_a += 1
                        else:
                            numSynDiff_a += 1
                            if aa2 != int1AA:
                                numNonSynDiff_a += 1
                            else:
                                numSynDiff_a += 1
                        if aa1 != int2AA:
                            numNonSynDiff_a += 1
                            if aa2 != int2AA:
                                numNonSynDiff_b += 1
                            else:
                                numSynDiff_b += 1
                        else:
                            numSynDiff_b += 1
                            if aa2 != int2AA:
                                numNonSynDiff_b += 1
                            else:
                                numSynDiff_b += 1
                        if int2AA != '*' and int1AA != '*':
                            if numSynDiff_a > numSynDiff_b:
                                numSynDiff = numSynDiff_a
                                numNonSynDiff = numNonSynDiff_a
                            elif numSynDiff_b > numSynDiff_a:
                                numSynDiff = numSynDiff_b
                                numNonSynDiff = numNonSynDiff_b
                            else:
                                numSynDiff = numSynDiff_a
                                numNonSynDiff = numNonSynDiff_a
                        elif int1AA == '*':
                            if int2AA == '*':
                                numSynDiff = numSynDiff_a
                                numNonSynDiff = numNonSynDiff_a
                            else:
                                numSynDiff = numSynDiff_b
                                numNonSynDiff = numNonSynDiff_b
                        else:
                            numSynDiff = numSynDiff_a
                            numNonSynDiff = numNonSynDiff_a
                    logfile = open('logfile_test.txt','a')
                    logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                    for allele in currCodon[1:]:
                        logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                    logfile.write(') -- ' + str(numSynDiff)  +' Synonymous Polymorphisms, ' + str(numNonSynDiff) + ' Nonsynyonmous Polymorphisms\n')
                    logfile.close()
                    numSynPolymorphisms += numSynDiff
                    numNonsynPolymorphisms += numNonSynDiff
                    numPolymorphisms += 2
                elif totalChanges == 3:
                    logfile = open('logfile_test.txt','a')
                    logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                    for allele in currCodon[1:]:
                        logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                    logfile.write(') -- has a complex evolutionary history\n')
                    logfile.close()         
    return invariantSites, numNonsynPolymorphisms, numSynPolymorphisms
    
    
def fixedDifferencesCodonsHaploid(seqDict,invariantSitesA,invariantSitesB,species,outgroup):
    referenceDict = {'MH':">NC_012920.1 Homo sapiens mitochondrion, complete genome",'N':'>AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome','D':">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia",'PoB':">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"}
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    startCodons = {'ATT':'M','ATC':'M','ATA':'M','ATG':'M','GTG':'M'}
    nonsynFixedDifferences = 0
    synFixedDifferences = 0
    codonDict,AADict = buildCodonDictHaploid(seqDict)
    speciesA = referenceDict[species]
    speciesB = referenceDict[outgroup]
    speciesASeq = seqDict[speciesA]
    speciesBSeq = seqDict[speciesB] 
    speciesACodons = codonDict[speciesA]
    speciesBCodons = codonDict[speciesB]
    if invariantSitesA == False:
        invariantSitesA = range(len(speciesASeq))
    if invariantSitesB == False:
        invariantSitesB = range(len(speciesBSeq))
    invariantSites = []
    for siteA in invariantSitesA:
        for siteB in invariantSitesB:
            if siteA == siteB and siteA not in invariantSites:
                invariantSites.append(siteA)
    comparison = species + '-' + outgroup
    variableCodons = []
    for site in invariantSites:
        currSite = [speciesASeq[site],speciesBSeq[site]]
        if currSite[0] != currSite[1] and '-' not in currSite and 'N' not in currSite:
            codonValue = site/3
            if codonValue not in variableCodons:
                variableCodons.append(site/3)
    for codonNumber in variableCodons:
        if codonNumber < len(speciesACodons):
            currCodon = [speciesACodons[codonNumber],speciesBCodons[codonNumber]]
            codon1 = speciesACodons[codonNumber]
            codon2 = speciesBCodons[codonNumber]
            if '-' not in codon1 and 'N' not in codon1 and '-' not in codon2 and 'N' not in codon2:
                if codonNumber == 0:
                    if codon1 in startCodons:
                        aa1 = startCodons[codon1]
                    else:
                        aa1 = geneticCode[codon1]
                    if codon2 in startCodons:
                        aa2 = startCodons[codon2]
                    else:
                        aa2 = geneticCode[codon2]
                else:
                    aa1 = geneticCode[codon1]
                    aa2 = geneticCode[codon2]
                AAs = [aa1,aa2]
                site1 = [codon1[0],codon2[0]]
                site2 = [codon1[1],codon2[1]]
                site3 = [codon1[2],codon2[2]]
                totalChanges = 0
                sites = [site1,site2,site3]
                for site in sites:
                    if site[0] != site[1]:
                        totalChanges += 1
                if totalChanges == 1:
                    if aa1 != aa2:
                        nonsynFixedDifferences += 1
                        logfile = open('logfile_test.txt','a')
                        logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + codon1 + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- ' + str(totalChanges) + ' Nonsynonymous Fixed Difference\n')
                        logfile.close()
                    else:
                        synFixedDifferences += 1    
                        logfile = open('logfile_test.txt','a')
                        logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + codon1 + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- ' + str(totalChanges) + ' Synonymous Fixed Difference\n') 
                        logfile.close()
                elif totalChanges > 1:
                    site1Num = codonNumber*3
                    site2Num = (codonNumber*3) + 1
                    site3Num = (codonNumber*3) + 2
                    if site1Num in invariantSites and site2Num in invariantSites and site3Num in invariantSites:
                        if totalChanges == 2:
                            if site1[0] == site1[1]:
                                int1 = site1[0] + site2[1] + site2[0]
                                int2 = site1[0] + site2[0] + site2[1]
                            elif site2[0] == site2[1]:
                                int1 = site1[0] + site2[0] + site3[1]
                                int2 = site1[1] + site2[0] + site3[0]
                            else:
                                int1 = site1[0] + site2[1] + site3[0]
                                int2 = site1[1] + site2[0] + site3[0]
                            if codonNumber == 0:
                                if int1 in startCodons:
                                    int1AA = startCodons[int1]
                                else:
                                    int1AA = geneticCode[int1]
                                if int2 in startCodons:
                                    int2AA = startCodons[int2]
                                else:
                                    int2AA = geneticCode[int2]
                            else:
                                int1AA = geneticCode[int1]
                                int2AA = geneticCode[int2]
                            numSynDiff_a = 0
                            numSynDiff_b = 0
                            numNonSynDiff_a = 0
                            numNonSynDiff_b = 0
                            if aa1 != int1AA:
                                numNonSynDiff_a += 1
                                if aa2 != int1AA:
                                    numNonSynDiff_a += 1
                                else:
                                    numSynDiff_a += 1
                            else:
                                numSynDiff_a += 1
                                if aa2 != int1AA:
                                    numNonSynDiff_a += 1
                                else:
                                    numSynDiff_a += 1
                            if aa1 != int2AA:
                                numNonSynDiff_b += 1
                                if aa2 != int2AA:
                                    numNonSynDiff_b += 1
                                else:
                                    numSynDiff_b += 1
                            else:
                                numSynDiff_b += 1
                                if aa2 != int2AA:
                                    numNonSynDiff_b += 1
                                else:
                                    numSynDiff_b += 1
                            if int2AA != '*' and int1AA != '*':
                                if numSynDiff_a > numSynDiff_b:
                                    numSynDiff = numSynDiff_a
                                    numNonSynDiff = numNonSynDiff_a
                                elif numSynDiff_b > numSynDiff_a:
                                    numSynDiff = numSynDiff_b
                                    numNonSynDiff = numNonSynDiff_b
                                else:
                                    numSynDiff = numSynDiff_a
                                    numNonSynDiff = numNonSynDiff_a
                            elif int1AA == '*':
                                if int2AA == '*':
                                    numSynDiff = numSynDiff_a
                                    numNonSynDiff = numNonSynDiff_a
                                else:
                                    numSynDiff = numSynDiff_b
                                    numNonSynDiff = numNonSynDiff_b
                            else:
                                numSynDiff = numSynDiff_a
                                numNonSynDiff = numNonSynDiff_a
                            logfile = open('logfile_test.txt','a')
                            logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + codon1 + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- ' + str(numSynDiff) + ' Synonymous Fixed Differences, ' + str(numNonSynDiff) + ' Nonsynonymous Fixed Differences\n')
                            logfile.close()
                            synFixedDifferences += numSynDiff
                            nonsynFixedDifferences += numNonSynDiff
                        elif totalChanges == 3:
                            logfile = open('logfile_test.txt','a')
                            logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + str(codon1) + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- has a complex evolutionary history')
                            logfile.close()        
                    else:
                        logfile = open('logfile_test.txt','a')
                        logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + str(codon1) + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- some sites in codon are polymorphic in 1 or more of compared species\n')  
                        logfile.close()
    return nonsynFixedDifferences,synFixedDifferences
        

def buildCodonDictHaploid(seqDict):
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    startCodons = ['ATT','ATC','ATA','ATG','GTG']
    codonDict = {}
    AADict = {}
    for seq in seqDict:
        nucleotideSeq = seqDict[seq]
        codonList = []
        i = 2
        while i < len(nucleotideSeq):
            currCodon = nucleotideSeq[i-2] + nucleotideSeq[i-1] + nucleotideSeq[i]
            codonList.append(currCodon)
            i += 3
        codonDict[seq] = codonList
        AAseq = ''
        codonNum = 1
        for codon in codonList:
            if codonNum == 1:
                if codon in startCodons:
                    aa = 'M'
                else:
                    aa = geneticCode[codon]
            elif codon in geneticCode:
                aa = geneticCode[codon]
            else:
                aa = 'X'
            AAseq += aa
            codonNum += 1
        if AAseq[-1] == '*':
            AAseq = AAseq[0:-1]
        AADict[seq] = AAseq
    return codonDict,AADict
    
def reverseComplementHaploid(seq):
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc
    
def seqDictGeneratorHaploid(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict, scaffoldList

def extracttRNASeqsHaploid(fasta):
    speciesDict = {'MH':[">NC_012920.1 Homo sapiens mitochondrion, complete genome",	">AF346963.1 Homo sapiens mitochondrion, complete genome",	">AF346964.1 Homo sapiens mitochondrion, complete genome",	">AF346965.1 Homo sapiens mitochondrion, complete genome",	">AF346966.1 Homo sapiens mitochondrion, complete genome",	">AF346967.1 Homo sapiens mitochondrion, complete genome",	">AF346968.1 Homo sapiens mitochondrion, complete genome",	">AF346969.1 Homo sapiens mitochondrion, complete genome",	">AF346970.1 Homo sapiens mitochondrion, complete genome",	">AF346971.1 Homo sapiens mitochondrion, complete genome",	">AF346972.1 Homo sapiens mitochondrion, complete genome",	">AF346973.1 Homo sapiens mitochondrion, complete genome",	">AF346974.1 Homo sapiens mitochondrion, complete genome",	">AF346975.1 Homo sapiens mitochondrion, complete genome",	">AF346976.1 Homo sapiens mitochondrion, complete genome",	">AF346977.1 Homo sapiens mitochondrion, complete genome",	">AF346978.1 Homo sapiens mitochondrion, complete genome",	">AF346979.1 Homo sapiens mitochondrion, complete genome",	">AF346980.1 Homo sapiens mitochondrion, complete genome",	">AF346981.1 Homo sapiens mitochondrion, complete genome",	">AF346982.1 Homo sapiens mitochondrion, complete genome",	">AF346983.1 Homo sapiens mitochondrion, complete genome",	">AF346984.1 Homo sapiens mitochondrion, complete genome",	">AF346985.1 Homo sapiens mitochondrion, complete genome",	">AF346986.1 Homo sapiens mitochondrion, complete genome",	">AF346987.1 Homo sapiens mitochondrion, complete genome",	">AF346988.1 Homo sapiens mitochondrion, complete genome",	">AF346989.1 Homo sapiens mitochondrion, complete genome",	">AF346990.1 Homo sapiens mitochondrion, complete genome",	">AF346991.1 Homo sapiens mitochondrion, complete genome",	">AF346992.1 Homo sapiens mitochondrion, complete genome",	">AF346993.1 Homo sapiens mitochondrion, complete genome",	">AF346994.1 Homo sapiens mitochondrion, complete genome",	">AF346995.1 Homo sapiens mitochondrion, complete genome",	">AF346996.1 Homo sapiens mitochondrion, complete genome",	">AF346997.1 Homo sapiens mitochondrion, complete genome",	">AF346998.1 Homo sapiens mitochondrion, complete genome",	">AF346999.1 Homo sapiens mitochondrion, complete genome",	">AF347000.1 Homo sapiens mitochondrion, complete genome",	">AF347001.1 Homo sapiens mitochondrion, complete genome",	">AF347002.1 Homo sapiens mitochondrion, complete genome",	">AF347003.1 Homo sapiens mitochondrion, complete genome",	">AF347004.1 Homo sapiens mitochondrion, complete genome",	">AF347005.1 Homo sapiens mitochondrion, complete genome",	">AF347006.1 Homo sapiens mitochondrion, complete genome",	">AF347007.1 Homo sapiens mitochondrion, complete genome",	">AF347008.1 Homo sapiens mitochondrion, complete genome",	">AF347009.1 Homo sapiens mitochondrion, complete genome",	">AF347010.1 Homo sapiens mitochondrion, complete genome",	">AF347011.1 Homo sapiens mitochondrion, complete genome",	">AF347012.1 Homo sapiens mitochondrion, complete genome",	">AF347013.1 Homo sapiens mitochondrion, complete genome",	">AF347014.1 Homo sapiens mitochondrion, complete genome",	">AF347015.1 Homo sapiens mitochondrion, complete genome"],'N':[">FM865407.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 1",	">FM865408.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 2",	">FM865409.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate El Sidron 1253",	">FM865410.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Vindija 33.25",	">FM865411.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Mezmaiskaya 1",	">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome"],'D':[">FR695060.1 Homo sp. Altai complete mitochondrial genome, isolate Denisova molar",	">KT780370.1 Homo sapiens ssp. Denisova isolate Denisova8 mitochondrion, complete genome",	">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia"],'PoB':[">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"]}
    referenceDict = {'MH':">NC_012920.1 Homo sapiens mitochondrion, complete genome",'N':">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome",'D':">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia",'PoB':">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"}
    speciesList = ['MH','N','D','PoB']
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    regions = {"D-Loop":{"type":"noncoding","region":"16049-592","direction":"forward"},	"TRNF":{"type":"tRNA","region":"593-663","direction":"forward"},	"RNR1":{"type":"rRNA","region":"664-1617","direction":"forward"},	"TRNV":{"type":"tRNA","region":"1618-1686","direction":"forward"},	"RNR2":{"type":"rRNA","region":"1687-3247","direction":"forward"},	"TRNL1":{"type":"tRNA","region":"3248-3322","direction":"forward"},	"ND1":{"type":"protein","region":"3325-4280","direction":"forward"},	"TRNI":{"type":"tRNA","region":"4281-4349","direction":"forward"},	"TRNQ":{"type":"tRNA","region":"4347-4418","direction":"complement"},	"TRNM":{"type":"tRNA","region":"4420-4487","direction":"forward"},	"ND2":{"type":"protein","region":"4488-5531","direction":"forward"},	"TRNW":{"type":"tRNA","region":"5530-5597","direction":"forward"},	"TRNA":{"type":"tRNA","region":"5605-5673","direction":"complement"},	"TRNN":{"type":"tRNA","region":"5675-5747","direction":"complement"},	"TRNC":{"type":"tRNA","region":"5779-5845","direction":"complement"},	"TRNY":{"type":"tRNA","region":"5845-5910","direction":"complement"},	"COI":{"type":"protein","region":"5929-7470","direction":"forward"},	"TRNS1":{"type":"tRNA","region":"7471-7539","direction":"complement"},	"TRND":{"type":"tRNA","region":"7543-7610","direction":"forward"},	"COII":{"type":"protein","region":"7611-8294","direction":"forward"},	"TRNK":{"type":"tRNA","region":"8320-8389","direction":"forward"},	"ATP8":{"type":"protein","region":"8391-8597","direction":"forward"},	"ATP6":{"type":"protein","region":"8552-9232","direction":"forward"},	"COIII":{"type":"protein","region":"9232-10015","direction":"forward"},	"TRNG":{"type":"tRNA","region":"10016-10083","direction":"forward"},	"ND3":{"type":"protein","region":"10084-10429","direction":"forward"},	"TRNR":{"type":"tRNA","region":"10430-10494","direction":"forward"},	"ND4L":{"type":"protein","region":"10495-10791","direction":"forward"},	"ND4":{"type":"protein","region":"10785-12162","direction":"forward"},	"TRNH":{"type":"tRNA","region":"12163-12231","direction":"forward"},	"TRNS2":{"type":"tRNA","region":"12232-12290","direction":"forward"},	"TRNL2":{"type":"tRNA","region":"12291-12361","direction":"forward"},	"ND5":{"type":"protein","region":"12362-14173","direction":"forward"},	"ND6":{"type":"protein","region":"14174-14698","direction":"complement"},	"TRNE":{"type":"tRNA","region":"14699-14767","direction":"complement"},	"CYTB":{"type":"protein","region":"14772-15912","direction":"forward"},	"TRNT":{"type":"tRNA","region":"15913-15978","direction":"forward"},	"TRNP":{"type":"tRNA","region":"15981-16048","direction":"complement"}, 'whole mt genome':{"type":"noncoding","region":"1-16599","direction":"forward"}}
    regionList = ["D-Loop",	"TRNF",	"RNR1",	"TRNV",	"RNR2",	"TRNL1",	"ND1",	"TRNI",	"TRNQ",	"TRNM",	"ND2",	"TRNW",	"TRNA",	"TRNN",	"TRNC",	"TRNY",	"COI",	"TRNS1",	"TRND",	"COII",	"TRNK",	"ATP8",	"ATP6",	"COIII",	"TRNG",	"ND3",	"TRNR",	"ND4L",	"ND4",	"TRNH",	"TRNS2",	"TRNL2",	"ND5",	"ND6",	"TRNE",	"CYTB",	"TRNT",	"TRNP",	"whole mt genome"]
    seqDict = {}
    seqList = []
    currSeq = ''
    seqName = ''
    alignFile = open(fasta,'r')
    for line in alignFile:
        if line[0] == '>':
            if seqName != '':
                seqDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
    	    seqList.append(seqName)
    	    currSeq = ''
        else:	
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    seqDict[seqName] = currSeq
    alignFile.close()
    for feature in regionList:
        items = regions[feature]
        direction = items['direction']
        region = items['region']
        regionSplit = region.split('-')
        regionStart = int(regionSplit[0]) - 1
        regionStop = int(regionSplit[1]) - 1
        currSeqDict = {}    
        for seq in seqList:
            if direction == 'forward':
                currSeq = seqDict[seq]
                regionSeq = currSeq[regionStart:regionStop]
                currSeqDict[seq] = regionSeq
            else:
                currSeq = seqDict[seq]
                regionSeq = currSeq[regionStart:regionStop + 1]
                newSeq = reverseComplementHaploid(regionSeq)
                currSeqDict[seq] = newSeq
        regionType = items['type']
        if regionType == 'tRNA':
            outfile = open(feature + '_seqs.fasta','w')
            for seq in seqList:
                outfile.write(seq + '\n')
                outfile.write(currSeqDict[seq] + '\n')
            outfile.close()

def extractrRNASeqsHaploid(fasta):
    speciesDict = {'MH':[">NC_012920.1 Homo sapiens mitochondrion, complete genome",	">AF346963.1 Homo sapiens mitochondrion, complete genome",	">AF346964.1 Homo sapiens mitochondrion, complete genome",	">AF346965.1 Homo sapiens mitochondrion, complete genome",	">AF346966.1 Homo sapiens mitochondrion, complete genome",	">AF346967.1 Homo sapiens mitochondrion, complete genome",	">AF346968.1 Homo sapiens mitochondrion, complete genome",	">AF346969.1 Homo sapiens mitochondrion, complete genome",	">AF346970.1 Homo sapiens mitochondrion, complete genome",	">AF346971.1 Homo sapiens mitochondrion, complete genome",	">AF346972.1 Homo sapiens mitochondrion, complete genome",	">AF346973.1 Homo sapiens mitochondrion, complete genome",	">AF346974.1 Homo sapiens mitochondrion, complete genome",	">AF346975.1 Homo sapiens mitochondrion, complete genome",	">AF346976.1 Homo sapiens mitochondrion, complete genome",	">AF346977.1 Homo sapiens mitochondrion, complete genome",	">AF346978.1 Homo sapiens mitochondrion, complete genome",	">AF346979.1 Homo sapiens mitochondrion, complete genome",	">AF346980.1 Homo sapiens mitochondrion, complete genome",	">AF346981.1 Homo sapiens mitochondrion, complete genome",	">AF346982.1 Homo sapiens mitochondrion, complete genome",	">AF346983.1 Homo sapiens mitochondrion, complete genome",	">AF346984.1 Homo sapiens mitochondrion, complete genome",	">AF346985.1 Homo sapiens mitochondrion, complete genome",	">AF346986.1 Homo sapiens mitochondrion, complete genome",	">AF346987.1 Homo sapiens mitochondrion, complete genome",	">AF346988.1 Homo sapiens mitochondrion, complete genome",	">AF346989.1 Homo sapiens mitochondrion, complete genome",	">AF346990.1 Homo sapiens mitochondrion, complete genome",	">AF346991.1 Homo sapiens mitochondrion, complete genome",	">AF346992.1 Homo sapiens mitochondrion, complete genome",	">AF346993.1 Homo sapiens mitochondrion, complete genome",	">AF346994.1 Homo sapiens mitochondrion, complete genome",	">AF346995.1 Homo sapiens mitochondrion, complete genome",	">AF346996.1 Homo sapiens mitochondrion, complete genome",	">AF346997.1 Homo sapiens mitochondrion, complete genome",	">AF346998.1 Homo sapiens mitochondrion, complete genome",	">AF346999.1 Homo sapiens mitochondrion, complete genome",	">AF347000.1 Homo sapiens mitochondrion, complete genome",	">AF347001.1 Homo sapiens mitochondrion, complete genome",	">AF347002.1 Homo sapiens mitochondrion, complete genome",	">AF347003.1 Homo sapiens mitochondrion, complete genome",	">AF347004.1 Homo sapiens mitochondrion, complete genome",	">AF347005.1 Homo sapiens mitochondrion, complete genome",	">AF347006.1 Homo sapiens mitochondrion, complete genome",	">AF347007.1 Homo sapiens mitochondrion, complete genome",	">AF347008.1 Homo sapiens mitochondrion, complete genome",	">AF347009.1 Homo sapiens mitochondrion, complete genome",	">AF347010.1 Homo sapiens mitochondrion, complete genome",	">AF347011.1 Homo sapiens mitochondrion, complete genome",	">AF347012.1 Homo sapiens mitochondrion, complete genome",	">AF347013.1 Homo sapiens mitochondrion, complete genome",	">AF347014.1 Homo sapiens mitochondrion, complete genome",	">AF347015.1 Homo sapiens mitochondrion, complete genome"],'N':[">FM865407.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 1",	">FM865408.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Feldhofer 2",	">FM865409.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate El Sidron 1253",	">FM865410.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Vindija 33.25",	">FM865411.1 Homo sapiens neanderthalensis complete mitochondrial genome, isolate Mezmaiskaya 1",	">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome"],'D':[">FR695060.1 Homo sp. Altai complete mitochondrial genome, isolate Denisova molar",	">KT780370.1 Homo sapiens ssp. Denisova isolate Denisova8 mitochondrion, complete genome",	">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia"],'PoB':[">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"]}
    referenceDict = {'MH':">NC_012920.1 Homo sapiens mitochondrion, complete genome",'N':">AM948965.1 Homo sapiens neanderthalensis complete mitochondrial genome",'D':">FN673705.1 Homo sp. Altai complete mitochondrial genome sequence from Denisova, Altai Russia",'PoB':">NC_023100.1 Homo heidelbergensis mitochondrion, complete genome"}
    speciesList = ['MH','N','D','PoB']
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    regions = {"D-Loop":{"type":"noncoding","region":"16049-592","direction":"forward"},	"TRNF":{"type":"tRNA","region":"593-663","direction":"forward"},	"RNR1":{"type":"rRNA","region":"664-1617","direction":"forward"},	"TRNV":{"type":"tRNA","region":"1618-1686","direction":"forward"},	"RNR2":{"type":"rRNA","region":"1687-3247","direction":"forward"},	"TRNL1":{"type":"tRNA","region":"3248-3322","direction":"forward"},	"ND1":{"type":"protein","region":"3325-4280","direction":"forward"},	"TRNI":{"type":"tRNA","region":"4281-4349","direction":"forward"},	"TRNQ":{"type":"tRNA","region":"4347-4418","direction":"complement"},	"TRNM":{"type":"tRNA","region":"4420-4487","direction":"forward"},	"ND2":{"type":"protein","region":"4488-5531","direction":"forward"},	"TRNW":{"type":"tRNA","region":"5530-5597","direction":"forward"},	"TRNA":{"type":"tRNA","region":"5605-5673","direction":"complement"},	"TRNN":{"type":"tRNA","region":"5675-5747","direction":"complement"},	"TRNC":{"type":"tRNA","region":"5779-5845","direction":"complement"},	"TRNY":{"type":"tRNA","region":"5845-5910","direction":"complement"},	"COI":{"type":"protein","region":"5929-7470","direction":"forward"},	"TRNS1":{"type":"tRNA","region":"7471-7539","direction":"complement"},	"TRND":{"type":"tRNA","region":"7543-7610","direction":"forward"},	"COII":{"type":"protein","region":"7611-8294","direction":"forward"},	"TRNK":{"type":"tRNA","region":"8320-8389","direction":"forward"},	"ATP8":{"type":"protein","region":"8391-8597","direction":"forward"},	"ATP6":{"type":"protein","region":"8552-9232","direction":"forward"},	"COIII":{"type":"protein","region":"9232-10015","direction":"forward"},	"TRNG":{"type":"tRNA","region":"10016-10083","direction":"forward"},	"ND3":{"type":"protein","region":"10084-10429","direction":"forward"},	"TRNR":{"type":"tRNA","region":"10430-10494","direction":"forward"},	"ND4L":{"type":"protein","region":"10495-10791","direction":"forward"},	"ND4":{"type":"protein","region":"10785-12162","direction":"forward"},	"TRNH":{"type":"tRNA","region":"12163-12231","direction":"forward"},	"TRNS2":{"type":"tRNA","region":"12232-12290","direction":"forward"},	"TRNL2":{"type":"tRNA","region":"12291-12361","direction":"forward"},	"ND5":{"type":"protein","region":"12362-14173","direction":"forward"},	"ND6":{"type":"protein","region":"14174-14698","direction":"complement"},	"TRNE":{"type":"tRNA","region":"14699-14767","direction":"complement"},	"CYTB":{"type":"protein","region":"14772-15912","direction":"forward"},	"TRNT":{"type":"tRNA","region":"15913-15978","direction":"forward"},	"TRNP":{"type":"tRNA","region":"15981-16048","direction":"complement"}, 'whole mt genome':{"type":"noncoding","region":"1-16599","direction":"forward"}}
    regionList = ["D-Loop",	"TRNF",	"RNR1",	"TRNV",	"RNR2",	"TRNL1",	"ND1",	"TRNI",	"TRNQ",	"TRNM",	"ND2",	"TRNW",	"TRNA",	"TRNN",	"TRNC",	"TRNY",	"COI",	"TRNS1",	"TRND",	"COII",	"TRNK",	"ATP8",	"ATP6",	"COIII",	"TRNG",	"ND3",	"TRNR",	"ND4L",	"ND4",	"TRNH",	"TRNS2",	"TRNL2",	"ND5",	"ND6",	"TRNE",	"CYTB",	"TRNT",	"TRNP",	"whole mt genome"]
    seqDict = {}
    seqList = []
    currSeq = ''
    seqName = ''
    alignFile = open(fasta,'r')
    for line in alignFile:
        if line[0] == '>':
            if seqName != '':
                seqDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
    	    seqList.append(seqName)
    	    currSeq = ''
        else:	
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    seqDict[seqName] = currSeq
    alignFile.close()
    for feature in regionList:
        items = regions[feature]
        direction = items['direction']
        region = items['region']
        regionSplit = region.split('-')
        regionStart = int(regionSplit[0]) - 1
        regionStop = int(regionSplit[1]) - 1
        currSeqDict = {}    
        for seq in seqList:
            if direction == 'forward':
                currSeq = seqDict[seq]
                regionSeq = currSeq[regionStart:regionStop]
                currSeqDict[seq] = regionSeq
            else:
                currSeq = seqDict[seq]
                regionSeq = currSeq[regionStart:regionStop + 1]
                newSeq = reverseComplementHaploid(regionSeq)
                currSeqDict[seq] = newSeq
        regionType = items['type']
        if regionType == 'rRNA':
            outfile = open(feature + '_seqs.fasta','w')
            for seq in seqList:
                outfile.write(seq + '\n')
                outfile.write(currSeqDict[seq] + '\n')
            outfile.close()


def readLogfileHaploid(logfile):
    infile  = open(logfile,'r')
    polyDict = {}
    fdDict = {}
    speciesList = ['MH','N','D','PoB']
    regionDict = {"D-Loop":(1,592),"TRNF":(593,663),"RNR1":(664,1617),"TRNV":(1618,1686),"RNR2":(1687,3247),"TRNL1":(3248,3322),"ND1":(3325,4280),"TRNI":(4281,4349),"TRNQ":(4347,4418),"TRNM":(4420,4487),"ND2":(4488,5531),"TRNW":(5530,5597),"TRNA":(5605,5673),"TRNN":(5675,5747),"TRNC":(5779,5845),"TRNY":(5845,5910),"COI":(5929,7470),"TRNS1":(7471,7539),"TRND":(7543,7610),"COII":(7611,8294),"TRNK":(8320,8389),"ATP8":(8391,8597),"ATP6":(8552,9232),"COIII":(9232,10015),"TRNG":(10016,10083),"ND3":(10084,10429),"TRNR":(10430,10494),"ND4L":(10495,10791),"ND4":(10785,12162),"TRNH":(12163,12231),"TRNS2":(12232,12290),"TRNL2":(12291,12361),"ND5":(12362,14173),"ND6":(14174,14698),"TRNE":(14699,14767),"CYTB":(14772,15912),"TRNT":(15913,15978),"TRNP":(15981,16048), "D-Loop2":(16048,16599)}
    regionList = ["D-Loop",	"TRNF",	"RNR1",	"TRNV",	"RNR2",	"TRNL1",	"ND1",	"TRNI",	"TRNQ",	"TRNM",	"ND2",	"TRNW",	"TRNA",	"TRNN",	"TRNC",	"TRNY",	"COI",	"TRNS1",	"TRND",	"COII",	"TRNK",	"ATP8",	"ATP6",	"COIII",	"TRNG",	"ND3",	"TRNR",	"ND4L",	"ND4",	"TRNH",	"TRNS2",	"TRNL2",	"ND5",	"ND6",	"TRNE",	"CYTB",	"TRNT",	"TRNP","D-Loop2"]
    for region in regionList:
        polyDict[region] = {}
        fdDict[region] = {}
        speciesNum = 1
        for speciesA in speciesList:
            currDict = polyDict[region]
            currDict[speciesA] = 0
            for speciesB in speciesList[speciesNum:]:
                comp = speciesA + '-' + speciesB
                currDict = fdDict[region]
                currDict[comp] = 0
            speciesNum += 1
    for line in infile:
        if line[0] != '~':
            realLine = line
            while realLine[-1] == '\n' or realLine[-1] == '\t' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            lineSplit = realLine.split(' -- ')
            changeType = lineSplit[1]
            mutationInfo = lineSplit[0]
            locations = lineSplit[2]
            locations = locations[2:-2]
            locations = locations.split(',')
            if len(locations) == 2:
                newLocations = []
                for item in locations:
                    if item[0] == "'":
                        newLocations.append(item[1:])
                    elif item[-1] == "'":
                        newLocations.append(item[0:-1])
                    else:
                        newLocations.append(item)
                locations = newLocations
            mutationInfoSplit = mutationInfo.split(' ')
            compType = mutationInfoSplit[0]
            compType = compType[0:-1]
            alleles = mutationInfoSplit[5:]
            changeCount = -1
            for allele in alleles:
                changeCount += 1
            if changeType == 'Polymorphism' and locations[0] != 'tergen':
                for location in locations:
                    if location in polyDict:
                        currDict = polyDict[location]
                        currDict[compType] += changeCount
                        polyDict[location] = currDict
                    else:
                        print location, changeCount
            elif changeType == 'Fixed Difference' and locations[0] != 'tergen':
                for location in locations:
                    if location in fdDict:
                        currDict = fdDict[location]
                        currDict[compType] += changeCount
                        fdDict[location] = currDict
    infile.close()
    outfile = open('polymorphicSites_Locations.txt','w')
    for region in regionList:
        outfile.write(region)
        regionPolyDict = polyDict[region]
        regionfdDict = fdDict[region]
        speciesNum = 1
        for speciesA in speciesList:
            for speciesB in speciesList[speciesNum:]:
                numFDs = regionfdDict[speciesA + '-' + speciesB]
                outfile.write('\t' + str(numFDs))
            speciesNum += 1
        for species in speciesList:
            numPolys = regionPolyDict[species]
            outfile.write('\t' + str(numPolys))
        outfile.write('\n')
            
                
def polymorphismsSubstitutions(fofn):
    speciesDict = {'MH':['>HuRef'],'D':['>DenisovaPinky'], 'N':['>AltaiNea']}
    referenceDict = {'MH':'>HuRef','D':'>DenisovaPinky', 'N':'>AltaiNea'}
    speciesList = ['MH', 'D', 'N']
    geneticCode = {"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"} # standard code
    alignmentFiles = open(fofn,'r')
    featurePolyDict = {}
    featureSubDict = {}
    regionList = []
    featureLengthsDict = {}
    for alignment in alignmentFiles:
        seqDict = {}
        seqList = []
        currSeq = ''
        seqName = ''
        alignFile = open(alignment[0:-1],'r')
        for line in alignFile:
            if line[0] == '>':
                if seqName != '':
                    seqDict[seqName] = currSeq
                seqName = line
                while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                    seqName = seqName[0:-1]
       	        seqList.append(seqName)
       	        currSeq = ''
            else:	
                currSeq += line
                while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                    currSeq = currSeq[0:-1]
        seqDict[seqName] = currSeq
        alignFile.close()
        fileNameSplit = alignment[0:-1].split('_')
        feature = fileNameSplit[0]
        featureLengthsDict[feature] = len(seqDict[seqList[0]])
        regionList.append(feature)
        polyDict = {}
        for speciesA in speciesList:
            currSpeciesList = speciesDict[speciesA]
            codonDict,AADict = buildCodonDict(seqDict)
            invariantSites = range(len(seqDict[currSpeciesList[0] + '_' + feature]))
            polyDict[speciesA] = [invariantSites]
            speciesAPolyDict = polyDict[speciesA]
            invariantSitesA = speciesAPolyDict[0]
        featurePolyDict[feature] = polyDict
        speciesNum = 1
        subDict = {}
        for speciesA in speciesList:
            for speciesB in speciesList[speciesNum:]:
                logfile = open('logfile.txt','a')
                logfile.write('~~~~~~~~~~~' + feature + ' – Fixed Differences ~~~~~~~~~~~\n')
                logfile.close()
                comparison = speciesA + '-' + speciesB
                speciesBPolyDict = polyDict[speciesB]
                invariantSitesB = speciesBPolyDict[0]
                nonsynFixedDifferences,synFixedDifferences = fixedDifferencesCodons(seqDict,invariantSitesA,invariantSitesB,speciesA,speciesB,feature)
                totalFixedDifferences = synFixedDifferences + nonsynFixedDifferences
                subDict[comparison] = [totalFixedDifferences,nonsynFixedDifferences,synFixedDifferences]
            speciesNum += 1
        featureSubDict[feature] = subDict
    outfile = open('polymorphismsSubsititutions_N-mts.txt','a')
    outfile.write('\t\t\tMH-N\t\t\tMH-D\t\t\tN-D\nRegion\tLength\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\tNonsynonymous Fixed Differences\tSynonymous Fixed Differences\tTotal Fixed Differences\n')
    for feature in regionList:
        outfile.write(feature + '\t')
        outfile.write(str(featureLengthsDict[feature]))
        polyDict = featurePolyDict[feature]
        subDict = featureSubDict[feature]
        speciesNum = 1
        for speciesA in speciesList:
            for speciesB in speciesList[speciesNum:]:
                comparison = speciesA + '-' + speciesB
                compNumbers = subDict[comparison]
                outfile.write('\t' + str(compNumbers[1]) + '\t' + str(compNumbers[2]) + '\t' + str(compNumbers[0]))
            speciesNum += 1
        for species in speciesList:
            if len(speciesDict[species]) > 1:
                polNumbers = polyDict[species]
                outfile.write('\t' + str(polNumbers[1]) + '\t' + str(polNumbers[2]) + '\t' + str(polNumbers[3]))
        outfile.write('\n')

    outfile.close()

        
def polymorphicCodons(seqDict,seqList):
    speciesDict = {'>HuRef': 'MH', '>DenisovaPinky': 'D', '>AltaiNea': 'N'}
    geneticCode = {"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"} # standard code
    startCodons = {'TTG':'M','CTG':'M','ATG':'M'}
    codonDict,AADict = buildCodonDict(seqDict) 
    currSpecies = speciesDict[seqList[0]]
    seqLength = len(seqDict[seqList[0]])
    numCodons = len(codonDict[seqList[0]])
    i = 0
    numPolymorphisms = 0
    numNonsynPolymorphisms = 0
    numSynPolymorphisms = 0
    invariantSites = []
    polymorphicCodons = []
    while i/3 < numCodons:
        currSite = []
        for seq in seqList:
            currSeq = seqDict[seq]
            if currSeq[i] not in currSite and currSeq[i] != '-' and currSeq[i] != 'N':
                currSite.append(currSeq[i])
        if len(currSite) > 1:
            numPolymorphisms += len(currSite)
            codonNum = i/3
            if codonNum not in polymorphicCodons:
                polymorphicCodons.append(codonNum)
        elif len(currSite) == 1:
            invariantSites.append(i)
        i += 1 
    for codonNumber in polymorphicCodons:
        currCodon = []
        for seq in seqList:
            codons = codonDict[seq]
            if codons[codonNumber] not in currCodon and '-' not in codons[codonNumber] and 'N' not in codons[codonNumber]:
                currCodon.append(codons[codonNumber])
        if len(currCodon) > 1:
            AAs = []
            site1 = []
            site2 = []
            site3 = []
            for polyCodon in currCodon:
                if polyCodon[0] not in site1:
                    site1.append(polyCodon[0])
                if polyCodon[1] not in site2:
                    site2.append(polyCodon[1])
                if polyCodon[2] not in site3:
                    site3.append(polyCodon[2])
                if codonNumber == 0:
                    if polyCodon in startCodons:
                        currAA = startCodons[polyCodon]
                    else:
                        currAA = geneticCode[polyCodon]
                else:
                    currAA = geneticCode[polyCodon]
                if currAA not in AAs:
                    AAs.append(currAA)
            totalChanges = (len(site1) - 1) + (len(site2) - 1) + (len(site3) - 1)
            if totalChanges == 1:
                if len(AAs) == 1:
                    numSynPolymorphisms += 1
                    logfile = open('logfile.txt','a')
                    logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                    for allele in currCodon[1:]:
                        logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                    logfile.write(') -- ' + str(totalChanges) + ' Synonymous Polymorphism\n')
                    logfile.close()
                else:
                    numNonsynPolymorphisms += 1
                    logfile = open('logfile.txt','a')
                    logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                    for allele in currCodon[1:]:
                        logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                    logfile.write(') -- ' + str(totalChanges) + ' Nonsynonymous Polymorphism\n')
                    logfile.close()
            elif totalChanges == 2:
                numSynDiff = 0
                numNonSynDiff = 0
                if len(currCodon) == 3:
                    codon1 = currCodon[0]
                    codon2 = currCodon[1]
                    codon3 = currCodon[2]
                    if codonNumber == 0:
                        if codon1 in startCodons:
                            aa1 = startCodons[codon1]
                        else:
                            aa1 = geneticCode[codon1]
                        if codon2 in startCodons:
                            aa2 = startCodons[codon2]
                        else:
                            aa2 = geneticCode[codon2]
                        if codon3 in startCodons:
                            aa3 = startCodons[codon3]
                        else:
                            aa3 = geneticCode[codon3]
                    else:  
                        aa1 = geneticCode[codon1]
                        aa2 = geneticCode[codon2]
                        aa3 = geneticCode[codon3]
                    codon1_2 = 0
                    codon1_3 = 0
                    codon2_3 = 0
                    k = 0
                    while k < 3:
                        if codon1[k] != codon2[k]:
                            codon1_2 += 1
                        if codon1[k] != codon3[k]:
                            codon1_3 += 1
                        if codon2[k] != codon3[k]:
                            codon2_3 += 1
                        k += 1
                    if codon1_2 > codon1_3 and codon1_2 > codon2_3:
                        if aa1 != aa3:
                            numNonSynDiff += 1
                        else:
                            numSynDiff += 1 
                        if aa2 != aa3:
                            numNonSynDiff += 1
                        else:
                            numSynDiff += 1
                    elif codon1_3 > codon1_2 and codon1_3 > codon2_3:
                        if aa2 != aa1:
                            numNonSynDiff += 1
                        else:
                            numSynDiff += 1 
                        if aa2 != aa3:
                            numNonSynDiff += 1
                        else:
                            numSynDiff += 1
                    elif codon2_3 > codon1_3 and codon2_3 > codon1_2:
                        if aa1 != aa2:
                            numNonSynDiff += 1
                        else:
                            numSynDiff += 1 
                        if aa1 != aa3:
                            numNonSynDiff += 1
                        else:
                            numSynDiff += 1
                    else:
                        logfile = open('logfile.txt','a')
                        logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' has a complex evolutionary history --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                        for allele in currCodon[1:]:
                            logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                        logfile.write(') -- Polymorphism\n')
                        logfile.close()                            
                elif len(currCodon) == 2:
                    codon1 = currCodon[0]
                    codon2 = currCodon[1]
                    if codonNumber == 0:
                        if codon1 in startCodons:
                            aa1 = startCodons[codon1]
                        else:
                            aa1 = geneticCode[codon1]
                        if codon2 in startCodons:
                            aa2 = startCodons[codon2]
                        else:
                            aa2 = geneticCode[codon2]
                    else:  
                        aa1 = geneticCode[codon1]
                        aa2 = geneticCode[codon2]
                    sites = [site1,site2,site3]
                    int1 = ''
                    int2 = ''
                    changeNum = 0
                    siteNum = 0
                    for site in sites:
                        if len(site) == 1:
                            int1 += site[0]
                            int2 += site[0]
                        elif changeNum == 0:
                            int1 += codon1[siteNum]
                            int2 += codon2[siteNum]
                            changeNum +=1
                        elif changeNum == 1:
                            int1 += codon2[siteNum]
                            int2 += codon1[siteNum]
                        siteNum += 1
                    int1AA = geneticCode[int1]
                    int2AA = geneticCode[int2]
                    numSynDiff_a = 0
                    numSynDiff_b = 0
                    numNonSynDiff_a = 0
                    numNonSynDiff_b = 0
                    if aa1 != int1AA:
                        numNonSynDiff_a += 1
                        if aa2 != int1AA:
                            numNonSynDiff_a += 1
                        else:
                            numSynDiff_a += 1
                    else:
                        numSynDiff_a += 1
                        if aa2 != int1AA:
                            numNonSynDiff_a += 1
                        else:
                            numSynDiff_a += 1
                    if aa1 != int2AA:
                        numNonSynDiff_a += 1
                        if aa2 != int2AA:
                            numNonSynDiff_b += 1
                        else:
                            numSynDiff_b += 1
                    else:
                        numSynDiff_b += 1
                        if aa2 != int2AA:
                            numNonSynDiff_b += 1
                        else:
                            numSynDiff_b += 1
                    if int2AA != '*' and int1AA != '*':
                        if numSynDiff_a > numSynDiff_b:
                            numSynDiff = numSynDiff_a
                            numNonSynDiff = numNonSynDiff_a
                        elif numSynDiff_b > numSynDiff_a:
                            numSynDiff = numSynDiff_b
                            numNonSynDiff = numNonSynDiff_b
                        else:
                            numSynDiff = numSynDiff_a
                            numNonSynDiff = numNonSynDiff_a
                    elif int1AA == '*':
                        if int2AA == '*':
                            numSynDiff = numSynDiff_a
                            numNonSynDiff = numNonSynDiff_a
                        else:
                            numSynDiff = numSynDiff_b
                            numNonSynDiff = numNonSynDiff_b
                    else:
                        numSynDiff = numSynDiff_a
                        numNonSynDiff = numNonSynDiff_a
                    logfile = open('logfile.txt','a')
                    logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                    for allele in currCodon[1:]:
                        logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                    logfile.write(') -- ' + str(numSynDiff)  +' Synonymous Polymorphisms, ' + str(numNonSynDiff) + ' Nonsynyonmous Polymorphisms\n')
                    logfile.close()
                    numSynPolymorphisms += numSynDiff
                    numNonsynPolymorphisms += numNonSynDiff
                    numPolymorphisms += 2
            elif totalChanges == 3:
                logfile = open('logfile.txt','a')
                logfile.write(currSpecies + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + currCodon[0] + ' [' + geneticCode[currCodon[0]] + ']')
                for allele in currCodon[1:]:
                    logfile.write(', ' + allele + ' [' + geneticCode[allele] + ']')
                logfile.write(') -- has a complex evolutionary history\n')
                logfile.close()         
    return invariantSites, numNonsynPolymorphisms, numSynPolymorphisms
    
    
def fixedDifferencesCodons(seqDict,invariantSitesA,invariantSitesB,species,outgroup,feature):
    referenceDict = {'MH':('>HuRef_' + feature),'D':('>DenisovaPinky_' + feature), 'N':('>AltaiNea_' + feature)}
    geneticCode = {"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"} # standard code
    startCodons = {'TTG':'M','CTG':'M','ATG':'M'}
    nonsynFixedDifferences = 0
    synFixedDifferences = 0
    codonDict,AADict = buildCodonDict(seqDict)
    speciesA = referenceDict[species]
    speciesB = referenceDict[outgroup]
    speciesASeq = seqDict[speciesA]
    speciesBSeq = seqDict[speciesB] 
    speciesACodons = codonDict[speciesA]
    speciesBCodons = codonDict[speciesB]
    if invariantSitesA == False:
        invariantSitesA = range(len(speciesASeq))
    if invariantSitesB == False:
        invariantSitesB = range(len(speciesBSeq))
    invariantSites = []
    for siteA in invariantSitesA:
        for siteB in invariantSitesB:
            if siteA == siteB and siteA not in invariantSites:
                invariantSites.append(siteA)
    comparison = species + '-' + outgroup
    variableCodons = []
    for site in invariantSites:
        currSite = [speciesASeq[site],speciesBSeq[site]]
        if currSite[0] != currSite[1] and '-' not in currSite and 'N' not in currSite:
            codonValue = site/3
            if codonValue not in variableCodons:
                variableCodons.append(site/3)
    for codonNumber in variableCodons:
        if codonNumber < len(speciesACodons):
            currCodon = [speciesACodons[codonNumber],speciesBCodons[codonNumber]]
            codon1 = speciesACodons[codonNumber]
            codon2 = speciesBCodons[codonNumber]
            if '-' not in codon1 and 'N' not in codon1 and '-' not in codon2 and 'N' not in codon2:
                if codonNumber == 0:
                    if codon1 in startCodons:
                        aa1 = startCodons[codon1]
                    else:
                        aa1 = geneticCode[codon1]
                    if codon2 in startCodons:
                        aa2 = startCodons[codon2]
                    else:
                        aa2 = geneticCode[codon2]
                elif codon1 in geneticCode and codon2 in geneticCode:
                    aa1 = geneticCode[codon1]
                    aa2 = geneticCode[codon2]
                else:
                    if codon1 in geneticCode:
                        aa1 = geneticCode[codon1]
                    else:
                        aa1 = 'X'
                    if codon2 in geneticCode:
                        aa2 = geneticCode[codon2]
                    else:
                        aa2 = 'X'
                if aa1 != 'X' and aa2 != 'X':
                    AAs = [aa1,aa2]
                    site1 = [codon1[0],codon2[0]]
                    site2 = [codon1[1],codon2[1]]
                    site3 = [codon1[2],codon2[2]]
                    totalChanges = 0
                    sites = [site1,site2,site3]
                    for site in sites:
                        if site[0] != site[1]:
                            totalChanges += 1
                    if totalChanges == 1:
                        if aa1 != aa2:
                            nonsynFixedDifferences += 1
                            logfile = open('logfile.txt','a')
                            logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + codon1 + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- ' + str(totalChanges) + ' Nonsynonymous Fixed Difference\n')
                            logfile.close()
                        else:
                            synFixedDifferences += 1    
                            logfile = open('logfile.txt','a')
                            logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + codon1 + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- ' + str(totalChanges) + ' Synonymous Fixed Difference\n') 
                            logfile.close()
                    elif totalChanges > 1:
                        site1Num = codonNumber*3
                        site2Num = (codonNumber*3) + 1
                        site3Num = (codonNumber*3) + 2
                        if site1Num in invariantSites and site2Num in invariantSites and site3Num in invariantSites:
                            if totalChanges == 2:
                                if site1[0] == site1[1]:
                                    int1 = site1[0] + site2[1] + site2[0]
                                    int2 = site1[0] + site2[0] + site2[1]
                                elif site2[0] == site2[1]:
                                    int1 = site1[0] + site2[0] + site3[1]
                                    int2 = site1[1] + site2[0] + site3[0]
                                else:
                                    int1 = site1[0] + site2[1] + site3[0]
                                    int2 = site1[1] + site2[0] + site3[0]
                                if codonNumber == 0:
                                    if int1 in startCodons:
                                        int1AA = startCodons[int1]
                                    else:
                                        int1AA = geneticCode[int1]
                                    if int2 in startCodons:
                                        int2AA = startCodons[int2]
                                    else:
                                        int2AA = geneticCode[int2]
                                else:
                                    int1AA = geneticCode[int1]
                                    int2AA = geneticCode[int2]
                                numSynDiff_a = 0
                                numSynDiff_b = 0
                                numNonSynDiff_a = 0
                                numNonSynDiff_b = 0
                                if aa1 != int1AA:
                                    numNonSynDiff_a += 1
                                    if aa2 != int1AA:
                                        numNonSynDiff_a += 1
                                    else:
                                        numSynDiff_a += 1
                                else:
                                    numSynDiff_a += 1
                                    if aa2 != int1AA:
                                        numNonSynDiff_a += 1
                                    else:
                                        numSynDiff_a += 1
                                if aa1 != int2AA:
                                    numNonSynDiff_b += 1
                                    if aa2 != int2AA:
                                        numNonSynDiff_b += 1
                                    else:
                                        numSynDiff_b += 1
                                else:
                                    numSynDiff_b += 1
                                    if aa2 != int2AA:
                                        numNonSynDiff_b += 1
                                    else:
                                        numSynDiff_b += 1
                                if int2AA != '*' and int1AA != '*':
                                    if numSynDiff_a > numSynDiff_b:
                                        numSynDiff = numSynDiff_a
                                        numNonSynDiff = numNonSynDiff_a
                                    elif numSynDiff_b > numSynDiff_a:
                                        numSynDiff = numSynDiff_b
                                        numNonSynDiff = numNonSynDiff_b
                                    else:
                                        numSynDiff = numSynDiff_a
                                        numNonSynDiff = numNonSynDiff_a
                                elif int1AA == '*':
                                    if int2AA == '*':
                                        numSynDiff = numSynDiff_a
                                        numNonSynDiff = numNonSynDiff_a
                                    else:
                                        numSynDiff = numSynDiff_b
                                        numNonSynDiff = numNonSynDiff_b
                                else:
                                    numSynDiff = numSynDiff_a
                                    numNonSynDiff = numNonSynDiff_a
                                logfile = open('logfile.txt','a')
                                logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + codon1 + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- ' + str(numSynDiff) + ' Synonymous Fixed Differences, ' + str(numNonSynDiff) + ' Nonsynonymous Fixed Differences\n')
                                logfile.close()
                                synFixedDifferences += numSynDiff
                                nonsynFixedDifferences += numNonSynDiff
                            elif totalChanges == 3:
                                logfile = open('logfile.txt','a')
                                logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + str(codon1) + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- has a complex evolutionary history')
                                logfile.close()        
                        else:
                            logfile = open('logfile.txt','a')
                            logfile.write(comparison + ': Codon Num – ' + str(codonNumber + 1) + ' --> (' + str(codon1) + ' [' + aa1 + '] <--> ' + codon2 +  ' [' + aa2 + ']) -- some sites in codon are polymorphic in 1 or more of compared species\n')  
                            logfile.close()
    return nonsynFixedDifferences,synFixedDifferences
        

def buildCodonDict(seqDict):
    geneticCode = {"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"} # standard code
    startCodons = ['TTG','CTG','ATG']
    codonDict = {}
    AADict = {}
    for seq in seqDict:
        nucleotideSeq = seqDict[seq]
        codonList = []
        i = 2
        while i < len(nucleotideSeq):
            currCodon = nucleotideSeq[i-2] + nucleotideSeq[i-1] + nucleotideSeq[i]
            codonList.append(currCodon)
            i += 3
        codonDict[seq] = codonList
        AAseq = ''
        codonNum = 1
        for codon in codonList:
            if codonNum == 1:
                if codon in startCodons:
                    aa = 'M'
                elif codon in geneticCode:
                    aa = geneticCode[codon]
                else:
                    aa = 'X'
            elif codon in geneticCode:
                aa = geneticCode[codon]
            else:
                aa = 'X'
            AAseq += aa
            codonNum += 1
        if AAseq[-1] == '*':
            AAseq = AAseq[0:-1]
        AADict[seq] = AAseq
    return codonDict,AADict
    
def reverseComplement(seq):
    #seqDict = buildSeqDict(fasta)
    #for sequence in seqDict:
        #seq = seqDict[sequence]
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        elif nuc == 'M':
            seq_revc = 'K' + seq_revc
        elif nuc == 'R':
            seq_revc = 'Y' + seq_revc
        elif nuc == 'S':
            seq_revc = 'S' + seq_revc
        elif nuc == 'W':
            seq_revc = 'W' + seq_revc
        elif nuc == 'K':
            seq_revc = 'M' + seq_revc
        elif nuc == 'Y':
            seq_revc = 'R' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc
    
def seqDictGenerator(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict, scaffoldList

def countSites(fofn):
    geneticCode = {"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"} # standard code
    masterFile = open(fofn,'r')
    for line in masterFile:
        fasta = line[0:-1]
        infile = open(fasta,'r')
        geneSplit = fasta.split('_')
        numSeqs = 0
        totalSynSites = 0
        totalNonsynSites = 0
        for line in infile:
            if line[0] != '>':
                codonList = []
                i = 2
                nucleotideSeq = line[0:-1]
                while i < len(nucleotideSeq):
                    currCodon = nucleotideSeq[i-2] + nucleotideSeq[i-1] + nucleotideSeq[i]
                    codonList.append(currCodon)
                    i += 3        
                for codon in codonList:
                    if 'N' in codon or '-' in codon:
                        totalSynSites += 0.71875
                        totalNonsynSites += 2.28125
                    else:
                        site1 = codon[0]
                        site2 = codon[1]
                        site3 = codon[2]
                        if site1 == 'A':
                            mut1 = 'C' + site2 + site3
                            mut2 = 'G' + site2 + site3
                            mut3 = 'T' + site2 + site3
                        elif site1 == 'C':
                            mut1 = 'A' + site2 + site3
                            mut2 = 'G' + site2 + site3
                            mut3 = 'T' + site2 + site3
                        elif site1 == 'G':
                            mut1 = 'A' + site2 + site3
                            mut2 = 'C' + site2 + site3
                            mut3 = 'T' + site2 + site3
                        elif site1 == 'T':
                            mut1 = 'A' + site2 + site3
                            mut2 = 'C' + site2 + site3
                            mut3 = 'G' + site2 + site3
                        if site2 == 'A':
                            mut4 = site1 + 'C' + site3
                            mut5 = site1 + 'G' + site3
                            mut6 = site1 + 'T' + site3
                        elif site2 == 'C':
                            mut4 = site1 + 'A' + site3
                            mut5 = site1 + 'G' + site3
                            mut6 = site1 + 'T' + site3
                        elif site2 == 'G':
                            mut4 = site1 + 'A' + site3
                            mut5 = site1 + 'C' + site3
                            mut6 = site1 + 'T' + site3
                        elif site2 == 'T':
                            mut4 = site1 + 'A' + site3
                            mut5 = site1 + 'C' + site3
                            mut6 = site1 + 'G' + site3
                        if site3 == 'A':
                            mut7 = site1 + site2 + 'C'
                            mut8 = site1 + site2 + 'G'
                            mut9 = site1 + site2 + 'T'
                        elif site3 == 'C':
                            mut7 = site1 + site2 + 'A'
                            mut8 = site1 + site2 + 'G'
                            mut9 = site1 + site2 + 'T'
                        elif site3 == 'G':
                            mut7 = site1 + site2 + 'A'
                            mut8 = site1 + site2 + 'C'
                            mut9 = site1 + site2 + 'T'
                        elif site3 == 'T':
                            mut7 = site1 + site2 + 'A'
                            mut8 = site1 + site2 + 'C'
                            mut9 = site1 + site2 + 'G'
                        aaList = [geneticCode[mut1],geneticCode[mut2],geneticCode[mut3], geneticCode[mut4],geneticCode[mut5],geneticCode[mut6],geneticCode[mut7],geneticCode[mut8],geneticCode[mut9]]
                        currAA = geneticCode[codon]
                        synSites = 0.0
                        nonsynSites = 0.0
                        for aa in aaList:
                            if aa == currAA:
                                synSites += 1.0
                            else:
                                nonsynSites += 1.0
                        synSites = synSites/3.0
                        nonsynSites = nonsynSites/3.0
                        totalSynSites += synSites
                        totalNonsynSites += nonsynSites
                numSeqs += 1
        meanSynSites = totalSynSites/numSeqs
        meanNonsynSites = totalNonsynSites/numSeqs
        print geneSplit[0] + str((meanSynSites,meanNonsynSites))
    
def workingCode(fasta):
    infile = open(fasta,'r')
    seqDict = {}
    for line in infile:
        if line[0] == '>':
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            if seqName!= '>Potamopyrgus_estuarinus':
                seqDict[seqName] = 'antipodarum'
            else:
                seqDict[seqName] = 'estuarinus'
    return seqDict


def buildCommand(fofn,ploidy='diploid'):
    if ploidy == 'diploid':
        polymorphismsSubstitutions(fofn)
    else:
        polymorphismsSubstitutionsHaploid(fofn)

helpStatement = 'To run program in unix terminal, type the following:\n\n\tpython polymorphismsSubstitutions_v4.py <input> <ploidy>\n\n\t\tinput\tIf ploidy is diploid, input should be a file of file names for all the alignments \n\t\t\tfor which you plan to obtain polymorphism and substitution counts. If input is \n\t\t\thaploid (e.g., mtDNA sequence), input should be a single alignment file in fasta \n\t\t\tformat. Note that in haploid mode, you will need to adjust the positional \n\t\t\tinformation in the two variables "regions" and "regionList" in the \n\t\t\tpolymorphismsSubstitutionsHaploid command. \n\n\t\tploidy\tPloidy is set to diploid on default. The two options are diploid and haploid. Note \n\t\t\tthat diploid state assumes the standerd genetic code and the haploid state assumes \n\t\t\tthe vertebrate mitochondrial genetic code. For diploid data, the program is \n\t\t\tcurrently set up to only determine differences between two sequences in pairwise \n\t\t\tfashion. If you wish to investigate polymorphism, please contact the developer.\n\nAssumptions:\n\n\t*\tFor noncoding data, this program assumes that for sites with multiple hits, the number of \n\t\tchanges = # alleles -1. \n\n\t*\tFor codons with multiple hits, the minimal number of changes possible to \n\t\texplain the codon with # \n\t\tchanges/codon = (# site 1 alleles - 1) + (# site 2 alleles - 1) + (# site 3 alleles - 1). \n\n\t*\tIf the number of changes was > 2, the codon will be sent to the logfile for manual \n\t\tdetermination. \n\n\t*\tIf there were three codons present and were all equally distantly related (e.g., AAA, AAT, \n\t\tAAC), the codon will be sent to the logfile for manual determination. \n\n\t*\tIf there aret hree codons and one codon was equally similar to the other two, but the other \n\t\ttwo codons were more dissimilar, the codon sharing the highest degree of similarity is \n\t\tassumed to be the 	evolutionary intermediate (e.g., AAA(K) <--> ACA(T) <--> ACT(T), \n\t\tyielding 1 synonymous change and 1 nonsynonymous change). \n\n\t*\tIf only two codons were present, but the codon contained multiple hits, the intermediate \n\t\trequiring the fewest number of nonsynonymous changes and no stop codons was assumed (e.g., \n\t\tTTG(L) <--> CTG(L) <--> CTC(L) = 2 syn changes, TTG(L) <--> TTC(F) <--> CTC(L) = 1 syn, 1 \n\t\tnonsyn change, meaning the first scenario would be assumed). \n\n\t*\tIf there were 2 or more differences in the codons of two difference species, but not all of \n\t\tthe changes were at invariant sites (i.e., at least one of the changes was polymorphic \n\t\twithin the species), it will be sent to the logfile for manual determination and only counted the \n\t\tnumber of fixed differences in the final count.\n\nBE SURE TO DOUBLE CHECK THE LOGFILE. IT IS THE DEFINITIVE OUTPUT. ALL COMPLEX CODONS ARE SENT THERE, BUT NOT \nINCLUDED IN OUTPUT TABLES\n\n Questions should be directed to Joel Sharbrough at jsharbro[at]gmail.com\n\nCopyright (c) 2017 Joel Sharbrough\n'
if len(sys.argv) > 2:
    buildCommand(sys.argv[1],sys.argv[2])
elif sys.argv[1] == 'help':
    sys.stdout.write(helpStatement)
else:
    buildCommand(sys.argv[1])