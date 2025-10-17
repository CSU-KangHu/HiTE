CXX= g++

CXXFLAGS = -std=c++11 -g -O3 -fopenmp -fmessage-length=0 # Linux
#CXXFLAGS = -Xpreprocessor -fopenmp -lomp  -std=c++17 -g -O3 -fmessage-length=0 #MacOS



#
# Objects
#

ORepeatsDetector  = ../bin/RepeatsDetector.o
#
OInvalidInputException  = ../bin/exception/InvalidInputException.o
OInvalidStateException  = ../bin/exception/InvalidStateException.o
OFileDoesNotExistException  = ../bin/exception/FileDoesNotExistException.o
OInvalidOrderOfOperationsException  = ../bin/exception/InvalidOrderOfOperationsException.o
OInvalidScoreException  = ../bin/exception/InvalidScoreException.o
OInvalidOperationException = ../bin/exception/InvalidOperationException.o

# Utility
OUtil = ../bin/utility/Util.o
OLCS = ../bin/utility/LCS.o
OLCSLen = ../bin/utility/LCSLen.o
OLCSubStr = ../bin/utility/LCSubStr.o
OLocation = ../bin/utility/Location.o
OEmptyLocation = ../bin/utility/EmptyLocation.o
OTSD = ../bin/utility/TSD.o
OTailFinder = ../bin/utility/TailFinder.o
OTail = ../bin/utility/Tail.o
OEmptyTail = ../bin/utility/EmptyTail.o
OEmptyTSD = ../bin/utility/EmptyTSD.o
OGlobAlignE = ../bin/utility/GlobAlignE.o
OLocAlign = ../bin/utility/LocAlign.o
OALignUtility = ../../bin/utility/AlignUtility.o


# Non TR
OChromosome  = ../bin/nonltr/Chromosome.o
OChromosomeOneDigit = ../bin/nonltr/ChromosomeOneDigit.o
OChromosomeRandom = ../bin/nonltr/ChromosomeRandom.o
OChromListMaker = ../bin/nonltr/ChromListMaker.o
OTableBuilder = ../bin/nonltr/TableBuilder.o
# OScorer = ../bin/nonltr/Scorer.o
# OScorerBase = ../bin/nonltr/ScorerBase.o
# OScorerIndex = ../bin/nonltr/ScorerIndex.o
# OScorerIndexEnd = ../bin/nonltr/ScorerIndexEnd.o
ODetector = ../bin/nonltr/Detector.o
ODetectorMaxima = ../bin/nonltr/DetectorMaxima.o
OChromDetector = ../bin/nonltr/ChromDetector.o
OChromDetectorMaxima =  ../bin/nonltr/ChromDetectorMaxima.o
ODetection = ../bin/nonltr/Detection.o
ORefiner = ../bin/nonltr/Refiner.o
OHMM = ../bin/nonltr/HMM.o
OHMMAutoRegressive = ../bin/nonltr/HMMAutoRegressive.o
OSampler = ../bin/nonltr/Sampler.o
OKMeans = ../bin/nonltr/KMeans.o
OScanner = ../bin/nonltr/Scanner.o
OTrainer = ../bin/nonltr/Trainer.o
OCluster = ../bin/nonltr/Cluster.o
OMOG = ../bin/nonltr/MOG.o
OLocationList = ../bin/nonltr/LocationList.o
OLocationListCollection = ../bin/nonltr/LocationListCollection.o

#
# TR
#
OScorerTr = ../bin/tr/ScorerTr.o
ODetectorTr = ../bin/tr/DetectorTr.o
OMatchTr = ../bin/tr/MatchTr.o
OFilterTr = ../bin/tr/FilterTr.o
OTr = ../bin/tr/Tr.o
OForwardTr = ../bin/tr/ForwardTr.o
OBackwardTr = ../bin/tr/BackwardTr.o
OCandidate = ../bin/tr/Candidate.o
OPairContainer = ../bin/tr/PairContainer.o
OTrKVisitor = ../bin/tr/TrKVisitor.o
OTrCsVisitor = ../bin/tr/TrCSVisitor.o
OLtrTe = ../bin/tr/LtrTe.o
OTrCollector = ../bin/tr/TrCollector.o
OTrSineVisitor = ../bin/tr/TrSineVisitor.o


# Test
OTestDetector = ../bin/TestDetector.o
OTestAssembler = ../bin/TestAssembler.o
OTestTr = ../bin/TestTr.o

# Stand Alone
OSamplerForReCon = ../bin/SamplerForReCon.o
ONovelRepeatExtractor = ../bin/NovelRepeatExtractor.o
ONovelRepeatLenCalculator = ../bin/NovelRepeatLenCalculator.o
ORandomGenomeMaker = ../bin/RandomGenomeMaker.o
OMarkovSequenceProcessor = ../bin/MarkovSequenceProcessor.o

# MS
OCompositionMatrix = ../bin/ms/CompositionMatrix.o
OTransitionCompositionMatrix = ../bin/ms/TransitionCompositionMatrix.o
OScoresMaker = ../bin/ms/ScoresMaker.o
OScoresMatrixMaker = ../bin/ms/ScoresMatrixMaker.o
OScoresBiasMaker = ../bin/ms/ScoresBiasMaker.o
OCoordinator = ../bin/ms/Coordinator.o
OMsDetector = ../bin/MsDetector.o
OMsDetectorOptimized = ../bin/MsDetectorOptimized.o

# $(OScorerIndex) $(OScorerBase) $(OScorerIndexEnd) $(OScorer) $(OSampler) $(OKMeans) $(OMOG)
OBJS = $(ORepeatsDetector) $(OInvalidInputException) $(OInvalidStateException) $(OFileDoesNotExistException) $(OInvalidOrderOfOperationsException) $(OInvalidOperationException) $(OInvalidScoreException) $(OUtil) $(OChromosome) $(OChromosomeRC) $(OChromosomeOneDigit) $(OTableBuilder)  $(ODetector) $(OChromDetector) $(OChromDetectorMaxima) $(ODetection) $(ODetectorMaxima) $(OLCS) $(OLCSLen) $(OLocation) $(OHMM)  $(OScanner) $(OTrainer) $(OChromListMaker) $(OCluster) $(OLocationList) $(OLocationListCollection)  $(OHMMAutoRegressive)

#
# Targets
#

TRepeatsDetector = ../bin/RepeatsDetector
TTestDetector = ../bin/TestDetector
TTestAssembler = ../bin/TestAssembler
TTestTr = ../bin/LtrDetector
TSamplerForReCon = ../bin/SamplerForReCon
TRandomGenomeMaker = ../bin/RandomGenomeMaker
TNovelRepeatExtractor = ../bin/NovelRepeatExtractor
TNovelRepeatLenCalculator = ../bin/NovelRepeatLenCalculator
TMsDetector = ../bin/MsDetector
TMsDetectorOptimized = ../bin/MsDetectorOptimized
TMarkovSequenceProcessor = ../bin/MarkovSequenceProcessor
TAlignUtility = ../../bin/AlignUtility


#
# Make RepeatsDetector
#

$(TRepeatsDetector): $(OBJS)
	$(CXX) -o $(TRepeatsDetector) $(OBJS)

#
# RepeatsDetector
#

# nonltr/KmerHashTable.cpp
$(ORepeatsDetector): RepeatsDetector.cpp nonltr/KmerHashTable.h nonltr/TableBuilder.h nonltr/HMM.h nonltr/Scanner.h nonltr/Trainer.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c RepeatsDetector.cpp -o $(ORepeatsDetector)

#
# Exception
#
$(OInvalidInputException): exception/InvalidInputException.cpp exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidInputException.cpp -o $(OInvalidInputException)

$(OInvalidStateException): exception/InvalidStateException.cpp exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidStateException.cpp -o $(OInvalidStateException)

$(OFileDoesNotExistException): exception/FileDoesNotExistException.cpp exception/FileDoesNotExistException.h
	$(CXX) $(CXXFLAGS) -c exception/FileDoesNotExistException.cpp -o $(OFileDoesNotExistException)

$(OInvalidOrderOfOperationsException): exception/InvalidOrderOfOperationsException.cpp exception/InvalidOrderOfOperationsException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidOrderOfOperationsException.cpp -o $(OInvalidOrderOfOperationsException)

$(OInvalidScoreException): exception/InvalidScoreException.cpp exception/InvalidScoreException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidScoreException.cpp -o $(OInvalidScoreException)

$(OInvalidOperationException): exception/InvalidOperationException.cpp exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidOperationException.cpp -o $(OInvalidOperationException)

#
# Utility
#

$(OUtil): utility/Util.cpp utility/Util.h utility/Location.h exception/FileDoesNotExistException.h
	$(CXX) $(CXXFLAGS) -c utility/Util.cpp -o $(OUtil)

$(OLCS): utility/LCS.cpp utility/LCS.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/LCS.cpp -o $(OLCS)

$(OLCSLen): utility/LCSLen.cpp utility/LCSLen.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/LCSLen.cpp -o $(OLCSLen)

$(OLCSubStr): utility/LCSubStr.cpp utility/LCSubStr.h utility/ILocation.h utility/Location.h
	$(CXX) $(CXXFLAGS) -c utility/LCSubStr.cpp -o $(OLCSubStr)

$(OLocation): utility/Location.cpp utility/Location.h utility/ILocation.h exception/InvalidInputException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c utility/Location.cpp -o $(OLocation)

$(OEmptyLocation): utility/EmptyLocation.cpp utility/EmptyLocation.h utility/ILocation.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c utility/EmptyLocation.cpp -o $(OEmptyLocation)

$(OTSD): utility/TSD.cpp utility/TSD.h utility/ITSD.h utility/ILocation.h utility/Location.h utility/EmptyLocation.h utility/Util.h utility/LCSubStr.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/TSD.cpp -o $(OTSD)

$(OEmptyTSD): utility/EmptyTSD.cpp utility/EmptyTSD.h utility/ITSD.h utility/ILocation.h utility/Util.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c utility/EmptyTSD.cpp -o $(OEmptyTSD)

$(OTailFinder): utility/TailFinder.cpp utility/TailFinder.h utility/ILocation.h utility/Util.h exception/InvalidInputException.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/TailFinder.cpp -o $(OTailFinder)

$(OTail): utility/Tail.cpp utility/Tail.h utility/Util.h utility/ITail.h
	$(CXX) $(CXXFLAGS) -c utility/Tail.cpp -o $(OTail)

$(OEmptyTail): utility/EmptyTail.h utility/EmptyTail.cpp utility/ITail.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c utility/EmptyTail.cpp -o $(OEmptyTail)

$(OGlobAlignE): utility/GlobAlignE.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/GlobAlignE.cpp -o $(OGlobAlignE)

$(OLocAlign): utility/LocAlign.h
	$(CXX) $(CXXFLAGS) -c utility/LocAlign.cpp -o $(OLocAlign)


$(OAlignUtility): ../tr/GlobAlignE.h
	$(CXX) $(CXXFlAGS) -c ../utility/AlignUtility.cpp -o $(OAlignUtility)

#
# Non LTR
#

$(OChromosome): nonltr/Chromosome.cpp nonltr/Chromosome.h nonltr/IChromosome.h utility/Util.h exception/InvalidInputException.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Chromosome.cpp -o $(OChromosome)

$(OChromosomeOneDigit): nonltr/ChromosomeOneDigit.cpp nonltr/ChromosomeOneDigit.h nonltr/Chromosome.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromosomeOneDigit.cpp	-o $(OChromosomeOneDigit)

$(OChromosomeRandom): nonltr/ChromosomeRandom.cpp nonltr/ChromosomeRandom.h nonltr/IChromosome.h exception/InvalidInputException.h exception/InvalidStateException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromosomeRandom.cpp -o $(OChromosomeRandom)

# nonltr/KmerHashTable.cpp
$(OTableBuilder): nonltr/TableBuilder.cpp nonltr/TableBuilder.h utility/Util.h nonltr/ChromosomeOneDigit.h nonltr/ITableView.h nonltr/KmerHashTable.h  nonltr/EnrichmentMarkovView.h nonltr/EnrichmentMarkovView.cpp exception/InvalidStateException.h nonltr/ChromListMaker.h nonltr/IChromosome.h
	$(CXX) $(CXXFLAGS) -c nonltr/TableBuilder.cpp -o $(OTableBuilder)

# $(OScorerBase): nonltr/ScorerBase.cpp nonltr/ScorerBase.h nonltr/ChromosomeOneDigit.h
#	$(CXX) $(CXXFLAGS) -c nonltr/ScorerBase.cpp -o $(OScorerBase)

# $(OScorer): nonltr/Scorer.cpp nonltr/Scorer.h nonltr/ScorerBase.h nonltr/ChromosomeOneDigit.h utility/Util.h exception/InvalidStateException.h
#	$(CXX) $(CXXFLAGS) -c nonltr/Scorer.cpp -o $(OScorer)

# $(OScorerIndex): nonltr/ScorerIndex.cpp nonltr/ScorerIndex.h nonltr/ScorerBase.h nonltr/ChromosomeOneDigit.h
#	$(CXX) $(CXXFLAGS) -c nonltr/ScorerIndex.cpp -o $(OScorerIndex)

# $(OScorerIndexEnd): nonltr/ScorerIndexEnd.cpp nonltr/ScorerIndexEnd.h nonltr/ScorerBase.h nonltr/ChromosomeOneDigit.h
#	$(CXX) $(CXXFLAGS) -c nonltr/ScorerIndexEnd.cpp -o $(OScorerIndexEnd)

$(ODetector): nonltr/Detector.cpp nonltr/Detector.h exception/InvalidInputException.h utility/Util.h
	$(CXX) def PAIRCONTAINER_H_$(CXXFLAGS) -c nonltr/Detector.cpp -o $(ODetector)

$(ODetectorMaxima): nonltr/DetectorMaxima.cpp nonltr/DetectorMaxima.h utility/ILocation.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/DetectorMaxima.cpp -o $(ODetectorMaxima)

$(OChromDetector): nonltr/ChromDetector.cpp nonltr/ChromDetector.h nonltr/Detector.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromDetector.cpp -o $(OChromDetector)

$(OChromDetectorMaxima): nonltr/ChromDetectorMaxima.cpp nonltr/ChromDetectorMaxima.h nonltr/DetectorMaxima.h nonltr/ChromosomeOneDigit.h utility/Util.h utility/ILocation.h utility/Location.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromDetectorMaxima.cpp -o $(OChromDetectorMaxima)

$(ODetection): nonltr/Detection.cpp nonltr/Detection.h
	$(CXX) $(CXXFLAGS) -c nonltr/Detection.cpp -o $(ODetection)

$(ORefiner): nonltr/Refiner.cpp nonltr/Refiner.h utility/ILocation.h
	$(CXX) $(CXXFLAGS) -c nonltr/Refiner.cpp -o $(ORefiner)

$(OHMM): nonltr/HMM.cpp nonltr/HMM.h utility/ILocation.h exception/InvalidStateException.h exception/InvalidInputException.h exception/FileDoesNotExistException.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c nonltr/HMM.cpp -o $(OHMM)

$(OHMMAutoRegressive): nonltr/HMMAutoRegressive.cpp nonltr/HMMAutoRegressive.h nonltr/KmerHashTable.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/HMMAutoRegressive.cpp -o $(OHMMAutoRegressive)

$(OSampler): nonltr/Sampler.cpp nonltr/Sampler.h nonltr/Chromosome.h nonltr/Scorer.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Sampler.cpp -o $(OSampler)

$(OKMeans): nonltr/KMeans.cpp nonltr/KMeans.h utility/Util.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/KMeans.cpp -o $(OKMeans)

$(OScanner): nonltr/Scanner.cpp nonltr/Scanner.h nonltr/Chromosome.h nonltr/ChromosomeOneDigit.h nonltr/HMM.h nonltr/HMMAutoRegressive.h nonltr/ITableView.h nonltr/Scorer.h nonltr/ScorerIndex.h nonltr/ScorerBase.h nonltr/ScorerIndexEnd.h utility/Util.h utility/ILocation.h exception/InvalidInputException.h exception/InvalidStateException.h exception/FileDoesNotExistException.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Scanner.cpp -o $(OScanner)

# nonltr/KmerHashTable.cpp	nonltr/Sampler.h nonltr/KMeans.h nonltr/MOG.h
$(OTrainer): nonltr/Trainer.cpp nonltr/Trainer.h nonltr/TableBuilder.h nonltr/KmerHashTable.h  nonltr/HMM.h nonltr/HMMAutoRegressive.h nonltr/ChromDetectorMaxima.h nonltr/Scorer.h nonltr/ScorerIndex.h nonltr/ScorerIndexEnd.h nonltr/ChromListMaker.h utility/Util.h nonltr/LocationListCollection.h exception/InvalidOperationException.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Trainer.cpp -o $(OTrainer)

$(OChromListMaker): nonltr/ChromListMaker.cpp nonltr/ChromListMaker.h nonltr/Chromosome.h nonltr/ChromosomeOneDigit.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromListMaker.cpp -o $(OChromListMaker)

$(OCluster): nonltr/Cluster.cpp nonltr/Cluster.h utility/Util.h exception/InvalidStateException.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Cluster.cpp -o $(OCluster)

$(OMOG): nonltr/MOG.cpp nonltr/MOG.h nonltr/Cluster.h utility/Util.h exception/InvalidInputException.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/MOG.cpp -o $(OMOG)

$(OLocationList): nonltr/LocationList.cpp nonltr/LocationList.h utility/ILocation.h utility/Location.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/LocationList.cpp -o $(OLocationList)

$(OLocationListCollection): nonltr/LocationListCollection.cpp nonltr/LocationListCollection.h utility/Location.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/LocationListCollection.cpp -o $(OLocationListCollection)

#
# Tr
#
# nonltr/KmerHashTable.cpp
$(OScorerTr): tr/ScorerTr.cpp tr/ScorerTr.h nonltr/KmerHashTable.h  nonltr/ChromosomeOneDigit.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c tr/ScorerTr.cpp -o $(OScorerTr)

$(OMatchTr): tr/MatchTr.cpp tr/MatchTr.h tr/ForwardTr.h utility/Util.h exception/InvalidInputException.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c tr/MatchTr.cpp -o $(OMatchTr)

$(ODetectorTr): tr/DetectorTr.cpp tr/DetectorTr.h utility/Util.h exception/InvalidInputException.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c tr/DetectorTr.cpp -o $(ODetectorTr)

$(OTr): tr/Tr.cpp tr/Tr.h tr/ITrVisitor.h exception/InvalidStateException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c tr/Tr.cpp -o $(OTr)

$(OForwardTr): tr/ForwardTr.cpp tr/ForwardTr.h tr/Tr.h exception/InvalidStateException.h exception/InvalidOperationException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c tr/ForwardTr.cpp -o $(OForwardTr)

$(OBackwardTr): tr/BackwardTr.cpp tr/BackwardTr.h tr/Tr.h exception/InvalidStateException.h exception/InvalidOperationException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c tr/BackwardTr.cpp -o $(OBackwardTr)

$(OCandidate): tr/Candidate.cpp tr/Candidate.h
	$(CXX) $(CXXFLAGS) -c tr/Candidate.cpp -o $(OCandidate)

$(OPairContainer): tr/PairContainer.cpp tr/PairContainer.h tr/Candidate.h
	$(CXX) $(CXXFLAGS) -c tr/PairContainer.cpp -o $(OPairContainer)

$(OTrKVisitor): tr/TrKVisitor.cpp tr/TrKVisitor.h tr/Tr.h
	$(CXX) $(CXXFLAGS) -c tr/TrKVisitor.cpp -o $(OTrKVisitor)

$(OTrCsVisitor): tr/TrCsVisitor.cpp tr/TrCsVisitor.h tr/Tr.h utility/LCS.h
	$(CXX) $(CXXFLAGS) -c tr/TrCsVisitor.cpp -o $(OTrCsVisitor)

$(OLtrTe): tr/LtrTe.cpp tr/LtrTe.h tr/BackwardTr.h utility/ITSD.h utility/TSD.h utility/EmptyTSD.h utility/ITail.h utility/Tail.h utility/Util.h utility/EmptyTail.h  exception/InvalidOperationException.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c tr/LtrTe.cpp -o $(OLtrTe)

$(OFilterTr):tr/FilterTr.h tr/BackwardTr.h tr/LtrTe.h tr/TrKVisitor.h tr/TrCsVisitor.h utility/Util.h utility/TSD.h utility/EmptyTSD.h utility/TailFinder.h utility/Tail.h utility/EmptyTail.h utility/LocAlign.h exception/InvalidStateException.h exception/InvalidInputException.h tr/TrSineVisitor.h
	$(CXX) $(CXXFLAGS) -c tr/FilterTr.cpp -o $(OFilterTr)

$(OTrCollector): tr/TrCollector.cpp tr/TrCollector.h tr/ScorerTr.h tr/DetectorTr.h tr/FilterTr.h tr/LtrTe.h nonltr/ChromosomeOneDigit.h
	$(CXX) $(CXXFLAGS) -c tr/TrCollector.cpp -o $(OTrCollector)

$(OTrSineVisitor): tr/TrSineVisitor.cpp tr/TrSineVisitor.h utility/Location.h utility/TailFinder.h utility/TSD.h
	$(CXX) $(CXXFLAGS) -c tr/TrSineVisitor.cpp -o $(OTrSineVisitor)


#
# Make Test Detector
#

$(OTestDetector): test/TestDetector.cpp nonltr/Detector.h nonltr/DetectorMaxima.h nonltr/ChromListMaker.h nonltr/ChromosomeOneDigit.h nonltr/ChromosomeRandom.h utility/Util.h utility/ILocation.h utility/TSD.h nonltr/MOG.h nonltr/Sampler.h nonltr/LocationList.h nonltr/LocationListCollection.h
	$(CXX) $(CXXFLAGS) -c test/TestDetector.cpp -o $(OTestDetector)

# $(OScorer) $(OScorer)
$(TTestDetector): $(OTestDetector) $(ODetector) $(ODetectorMaxima) $(OInvalidInputException) $(OUtil) $(OLocation) $(OFileDoesNotExistException) $(OInvalidStateException) $(OChromListMaker) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom)  $(OChromDetectorMaxima) $(OHMM) $(OSampler) $(OKMeans) $(OTSD) $(OLCSubStr) $(OEmptyLocation) $(OInvalidOperationException) $(OCluster) $(OMOG) $(OChromosomeRC) $(OLocationList) $(OLocationListCollection)
	$(CXX) $(CXXFLAGS) -o $(TTestDetector) $(OTestDetector) $(ODetector) $(ODetectorMaxima) $(OInvalidInputException) $(OUtil) $(OLocation) $(OFileDoesNotExistException) $(OInvalidStateException) $(OChromListMaker) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom)  $(OChromDetectorMaxima) $(OHMM) $(OSampler) $(OKMeans) $(OTSD) $(OLCSubStr) $(OEmptyLocation) $(OInvalidOperationException) $(OCluster) $(OMOG) $(OChromosomeRC) $(OLocationList) $(OLocationListCollection)

#
# Make Test Assember
#
$(OTestAssembler): test/TestAssembler.cpp nonltr/ChromosomeOneDigit.h utility/TailFinder.h
	$(CXX) $(CXXFLAGS) -c test/TestAssembler.cpp -o $(OTestAssembler)

$(TTestAssembler): $(OTestAssembler) $(OTailFinder) $(OChromosomeOneDigit) $(OChromosome) $(OInvalidInputException) $(OUtil) $(OFileDoesNotExistException)
	$(CXX) $(CXXFLAGS) -o $(TTestAssembler) $(OTestAssembler) $(OTailFinder) $(OChromosomeOneDigit) $(OChromosome) $(OInvalidInputException) $(OUtil) $(OFileDoesNotExistException) $(OInvalidStateException) $(OLocation)

#
# Make Test TR
#

$(OTestTr): nonltr/Chromosome.h nonltr/ChromosomeRandom.h nonltr/ChromosomeOneDigit.h tr/ScorerTr.h tr/DetectorTr.h tr/FilterTr.h tr/Tr.h tr/ForwardTr.h tr/BackwardTr.h tr/TrCollector.h utility/Util.h utility/LCS.h utility/LCSLen.h utility/Location.h utility/ILocation.h utility/LCSubStr.h utility/TSD.h nonltr/ChromListMaker.h
	$(CXX) $(CXXFLAGS) -c test/TestTr.cpp -o $(OTestTr)

$(TTestTr): $(OTestTr) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException) $(OTr) $(OBackwardTr) $(OForwardTr) $(OCandidate) $(OPairContainer) $(OScorerTr) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom) $(OTrKVisitor) $(OTrCsVisitor) $(OTrSineVisitor) $(OFilterTr) $(OLtrTe) $(ODetectorTr) $(OMatchTr) $(OUtil) $(OLocation) $(OEmptyLocation) $(OLCS) $(OLCSLen) $(OLCSubStr) $(OTSD) $(OEmptyTSD) $(OTailFinder) $(OTail) $(OEmptyTail) $(OTrCollector) $(OGlobAlignE) $(OLocAlign) $(OChromListMaker)
	$(CXX) $(CXXFLAGS) -o $(TTestTr) $(OTestTr) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException) $(OTr) $(OBackwardTr) $(OForwardTr) $(OCandidate) $(OPairContainer) $(OScorerTr) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom) $(OTrKVisitor) $(OTrCsVisitor) $(OTrSineVisitor) $(OFilterTr) $(OLtrTe) $(ODetectorTr) $(OMatchTr) $(OUtil) $(OLocation) $(OEmptyLocation) $(OLCS) $(OLCSLen) $(OLCSubStr) $(OTSD) $(OEmptyTSD) $(OTailFinder) $(OTail) $(OEmptyTail) $(OTrCollector) $(OGlobAlignE) $(OLocAlign) $(OChromListMaker) -lstdc++fs

#
# Make Stand alone - SamplerForReCon
#
$(OSamplerForReCon): Utility/Util.h nonltr/Chromosome.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c standAlone/SamplerForReCon.cpp -o $(OSamplerForReCon)

$(TSamplerForReCon): $(OSamplerForReCon) $(OLocation) $(OUtil) $(OChromosome) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException)
	$(CXX) $(CXXFLAGS) -o $(TSamplerForReCon) $(OSamplerForReCon) $(OLocation) $(OUtil) $(OChromosome) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException)

$(ONovelRepeatExtractor): nonltr/Chromosome.h
	$(CXX) $(CXXFLAGS) -c standAlone/NovelRepeatExtractor.cpp -o $(ONovelRepeatExtractor)

$(TNovelRepeatExtractor): $(ONovelRepeatExtractor) $(OChromosome) $(OInvalidInputException) $(OInvalidOperationException) $(OFileDoesNotExistException) $(OUtil) $(OLocation)
	$(CXX) $(CXXFLAGS) -o $(TNovelRepeatExtractor) $(ONovelRepeatExtractor) $(OChromosome) $(OInvalidInputException) $(OInvalidOperationException) $(OFileDoesNotExistException) $(OUtil) $(OLocation)

$(ONovelRepeatLenCalculator): nonltr/Chromosome.h
	$(CXX) $(CXXFLAGS) -c standAlone/NovelRepeatLenCalculator.cpp -o $(ONovelRepeatLenCalculator)

$(TNovelRepeatLenCalculator): $(ONovelRepeatLenCalculator) $(OChromosome) $(OInvalidInputException) $(OInvalidOperationException) $(OFileDoesNotExistException) $(OUtil) $(OLocation)
	$(CXX) $(CXXFLAGS) -o $(TNovelRepeatLenCalculator) $(ONovelRepeatLenCalculator) $(OChromosome) cmake vs automake$(OInvalidInputException) $(OInvalidOperationException) $(OFileDoesNotExistException) $(OUtil) $(OLocation)

$(ORandomGenomeMaker): utility/Util.h nonltr/Chromosome.h nonltr/ChromosomeRandom.h
	$(CXX) $(CXXFLAGS) -c standAlone/RandomGenomeMaker.cpp -o $(ORandomGenomeMaker)

$(TRandomGenomeMaker): $(ORandomGenomeMaker) $(OChromosome) $(OChromosomeRandom) $(OLocation) $(OUtil) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException)
	$(CXX) $(CXXFLAGS) -o $(TRandomGenomeMaker) $(ORandomGenomeMaker) $(OChromosome) $(OChromosomeRandom) $(OLocation) $(OUtil) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException)

$(OMarkovSequenceProcessor): utility/Util.h nonltr/Chromosome.h
	$(CXX) $(CXXFLAGS) -c standAlone/MarkovSequenceProcessor.cpp -o $(OMarkovSequenceProcessor)

$(TMarkovSequenceProcessor): $(OMarkovSequenceProcessor) $(OChromosome) $(OLocation) $(OUtil) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException)
	$(CXX) $(CXXFLAGS) -o $(TMarkovSequenceProcessor) $(OMarkovSequenceProcessor) $(OChromosome) $(OLocation) $(OUtil) $(OInvalidInputException) $(OFileDoesNotExistException) $(OInvalidStateException) $(OInvalidOperationException)

#
# MS
#
$(OCompositionMatrix): ms/CompositionMatrix.cpp ms/CompositionMatrix.h utility/Util.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c ms/CompositionMatrix.cpp -o $(OCompositionMatrix)

$(OTransitionCompositionMatrix): ms/TransitionCompositionMatrix.cpp ms/TransitionCompositionMatrix.h ms/CompositionMatrix.h
	$(CXX) $(CXXFLAGS) -c ms/TransitionCompositionMatrix.cpp -o $(OTransitionCompositionMatrix)

$(OScoresMaker): ms/ScoresMaker.cpp ms/ScoresMaker.h utility/Util.h exception/InvalidScoreException.h
	$(CXX) $(CXXFLAGS) -c ms/ScoresMaker.cpp -o $(OScoresMaker)

$(OScoresMatrixMaker): ms/ScoresMatrixMaker.cpp ms/ScoresMatrixMaker.h ms/ScoresMaker.h
	$(CXX) $(CXXFLAGS) -c ms/ScoresMatrixMaker.cpp -o $(OScoresMatrixMaker)

$(OScoresBiasMaker): ms/ScoresBiasMaker.cpp ms/ScoresBiasMaker.h ms/ScoresMaker.h
	$(CXX) $(CXXFLAGS) -c ms/ScoresBiasMaker.cpp -o $(OScoresBiasMaker)

$(OCoordinator): ms/Coordinator.cpp ms/Coordinator.h ms/ScoresMaker.h ms/ScoresMatrixMaker.h ms/ScoresBiasMaker.h ms/CompositionMatrix.h exception/InvalidOrderOfOperationsException.h
	$(CXX) $(CXXFLAGS) -msse4 $(IDIR) -c ms/Coordinator.cpp -o $(OCoordinator)

$(OMsDetector): MsDetector.cpp ms/Coordinator.h ms/CompositionMatrix.h ms/TransitionCompositionMatrix.h utility/Util.h
	$(CXX) $(CXXFLAGS) -msse4 $(IDIR) -c MsDetector.cpp -o $(OMsDetector)

$(OMsDetectorOptimized): MsDetectorOptimized.cpp ms/Coordinator.h utility/Util.h
	$(CXX) $(CXXFLAGS) -msse4 $(IDIR) -c MsDetectorOptimized.cpp -o $(OMsDetectorOptimized)

$(TMsDetector): $(OScoresMaker) $(OScoresMatrixMaker) $(OScoresBiasMaker) $(OCoordinator) $(OInvalidOrderOfOperationsException) $(OInvalidScoreException) $(OUtil) $(OFileDoesNotExistException) $(OInvalidInputException) $(OCompositionMatrix) $(OTransitionCompositionMatrix) $(OMsDetector)
	$(CXX) $(CXXFLAGS) -msse4 -rdynamic $(LDIR) $(WL) -o $(TMsDetector) $(OScoresMaker) $(OScoresMatrixMaker) $(OScoresBiasMaker) $(OCoordinator) $(OUtil) $(OInvalidOrderOfOperationsException) $(OInvalidScoreException) $(OFileDoesNotExistException) $(OInvalidInputException) $(OCompositionMatrix) $(OTransitionCompositionMatrix) $(OMsDetector)

$(TMsDetectorOptimized): $(OScoresMaker) $(OScoresMatrixMaker) $(OScoresBiasMaker) $(OCoordinator)  $(OInvalidOrderOfOperationsException) $(OInvalidScoreException) $(OUtil) $(OFileDoesNotExistException) $(OInvalidInputException) $(OCompositionMatrix) $(OTransitionCompositionMatrix) $(OMsDetectorOptimized)
	$(CXX) $(CXXFLAGS) -msse4 -rdynamic $(LDIR) $(WL) -o $(TMsDetectorOptimized) $(OScoresMaker) $(OScoresMatrixMaker) $(OScoresBiasMaker) $(OCoordinator) $(OUtil) $(OInvalidOrderOfOperationsException) $(OInvalidScoreException) $(OFileDoesNotExistException) $(OInvalidInputException) $(OCompositionMatrix) $(OTransitionCompositionMatrix) $(OMsDetectorOptimized)

#
# Alignment Utility
#
$(TAlignUtility): $(OGlobAlignE)
	$(CXX) $(CXXFLAGS) -o $(TAlignUtility) $(OGlobAlignE)

#
# Make all
#

all: $(TRepeatsDetector) $(TTestDetector)  $(TTestAssembler) $(TMsDetector) $(TMsDetectorOptimized)
detector: $(TRepeatsDetector)
assembler: $(TTestAssembler)
samplerForReCon: $(TSamplerForReCon)
novel: $(TNovelRepeatExtractor)
novelLen: $(TNovelRepeatLenCalculator)
tr: $(TTestTr)
test: $(TTestDetector)
randomGenome: $(TRandomGenomeMaker)
processor: $(TMarkovSequenceProcessor)
util: $(TAlignUtility)

#
# Make initial directories
#
bin:
	mkdir ../bin
	mkdir ../bin/exception
	mkdir ../bin/utility
	mkdir ../bin/nonltr
	mkdir ../bin/tr
	mkdir ../bin/test

#
# Make clean
#

clean:
	rm -f ../bin/*.o ../bin/exception/*.o ../output/*.o ../bin/ms/*.o ../bin/nonltr/*.o ../bin/test/*.o ../bin/utility/*.o ../bin/tr/*.o *.o $(TRepeatsDetector) $(TTestDetector) $(TMsDetector) $(TMsDetectorOptimized) $(TTestAssembler) $(TTestTr) $(TSamplerForReCon) $(TRandomGenomeMaker)
