/*
 * TestTr.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: Hani Zakaria Girgis, PhD and Joseph Valencia
 */

#include "../nonltr/ChromosomeOneDigit.h"
#include "../nonltr/Chromosome.h"
#include "../nonltr/ChromosomeRandom.h"
#include "../nonltr/ChromListMaker.h" 

#include "../tr/ScorerTr.h"
#include "../tr/FilterTr.h"
#include "../tr/Tr.h"
#include "../tr/ForwardTr.h"
#include "../tr/BackwardTr.h"
#include "../tr/TrCollector.h"
#include "../utility/Util.h"
#include "../utility/LCS.h"
#include "../utility/LCSLen.h"
#include "../utility/Location.h"
#include "../utility/ILocation.h"
#include "../utility/LCSubStr.h"
#include "../utility/TSD.h"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <experimental/filesystem>


using namespace std;
using namespace nonltr;
using namespace utility;
using namespace tr;


namespace fs = std::experimental::filesystem;

int getOptionalArg(vector<string> * args, string option){
    string answer = "-1";
	auto it = std::find(args->begin(),args->end(),option);

	if(it!=args->end() & ++it !=args->end()){
		answer = * it;
	}
	return stoi(answer);
}

int main(int argc, char * argv[]) {
	int minDist = 400;
	int maxDist = 22000;
	int minLenLTR = 100;
	int maxLenLTR = 6000;
	int identity = 85;
	int k = 13;
	int minPlateauLength = 10;
	int gapTol = 200;
	int threads = 1;
	string chromDir ="";
	string outputDir = "";
	bool seqLevel = false;
	bool printRawScores = false;
	bool printCleanScores = false;
	bool displayHelp = false;
	bool bedOutput = false;
	bool displayNested = false;


	string helpMessage = 
	"| -arg     | Description | Default |\n"
	"| ---------------- | ----------- | ------- |\n"
	"| -fasta | Directory containing files to be scanned. Files must have .fa extension. | required |\n"
	"| -destDir | Output directory where the results are stored. | required |\n"
	"( IMPORTANT: Files under the output directory are deleted at the start of the program.)\n"
	"| -minLen | Minimum length of complete LTR-RTs. Constrains scoring system and filters. | 400 |\n"
	"| -maxLen |  Maximum length of complete LTR-RTs. Constrains scoring system and filters. | 22000 |\n"
	"| -minLenLTR | Minimum length of LTR direct repeat. Constrains filters. | 100 |\n"
	"| -maxLenLTR | Maximum length of LTR direct repeat. Constrains filters. | 6000 |\n"
	"( Note run time is highly dependent on this parameter, as it provides an upper bound for alignment length in the boundary adjustment step.)\n"               
	"| -id | Minimum identity [0-100] between 5' and 3' LTRs. | 85 |\n"
	"| -k  | Length of k-mers to adjust scoring system. Tradeoff between noise and resistance to mutation. | 13 |\n"
	"| -plateauSeed | Minimum length of plateaus to be initially considered 'Keep' in merging step. | 10 |\n"
	"| -nThreads | Number of cores to be used. | 1 |\n"
	"| -gapTol | Number of base pairs that two plateaus can differ by in height/distance. Affects both plateau merging and pairing steps. | 200 |\n"
	"|-seqLevel| Forces parallel execution on sequences within multi-FASTA file. Loads all sequences into memory | disabled |\n"
	"| -rawScores | prints the raw scores to a file called xxxxRawScores.txt under the output directory. | disabled |\n"
	"| -cleanedScores | prints the scores after merging to a file called xxxxCleanedScores.txt under the output directory. | disabled |\n"
	"| -nested | searches for nested elements. Results are stored in seperate files (marked as xxxxNestedDetector.bed) under the output directory | disabled |\n"
	"| -bedFormat | prints BED format without additional annotations (PPT and TSD). | disabled |\n"
	"| -help | prints this help message. | disabled |\n"
	;

	std::vector<string> * argList = new vector<string>();

	for(int i = 1;i<argc;i++){

		argList->push_back(string(argv[i]));
	}

	auto help = std::find(argList->begin(),argList->end(),"-help");

	if(help!=argList->end()){
		cout <<helpMessage<<endl;
		exit(1);
	}

	auto printRaw = std::find(argList->begin(),argList->end(),"-rawScores");

	if(printRaw!=argList->end()){

		printRawScores = true;
	}

	auto printClean = std::find(argList->begin(),argList->end(),"-cleanedScores");

		if(printClean!=argList->end()){

		printCleanScores = true;
	}


	auto bed = std::find(argList->begin(),argList->end(),"-bedFormat");

			if(bed!=argList->end()){

			bedOutput= true;
			cout <<"Bed format activated"<<endl;
		}

	auto nest = std::find(argList->begin(),argList->end(),"-nested");

			if(nest!=argList->end()){

			displayNested= true;
			
		}

	auto src = std::find(argList->begin(),argList->end(),"-fasta");

	if(src!=argList->end() & ++src !=argList->end()){

		 chromDir = *src;
	}

	else{
		cout<<" \'-fasta\' is a required argument"<<endl;
		exit(1);
	}

	auto dest = std::find(argList->begin(),argList->end(),"-destDir");

	if(dest!=argList->end() & ++dest != argList->end()){

		outputDir = *dest;
	}
	else{
		cout<<" \'-destDir\' is a required argument"<<endl;
		exit(1);
	}

	int test = getOptionalArg(argList, "-minLen");

	if(test!=-1){
		minDist = test;
	}

	test = getOptionalArg(argList, "-maxLen");

	if(test!=-1){
		maxDist = test;
	}

	test = getOptionalArg(argList,"-minLenLTR");

	if(test!=-1){

		minLenLTR = test;
	}

	test = getOptionalArg(argList, "-maxLenLTR");

	if(test!=-1){
		maxLenLTR = test;
	}

	test = getOptionalArg(argList,"-id");

	if(test!= -1){

		identity = test;
	}

	test = getOptionalArg(argList, "-k");

	if(test!=-1){
		k = test;
	}

	test = getOptionalArg(argList, "-plateauSeed");

	if(test!=-1){
		minPlateauLength = test;
	}

	test = getOptionalArg(argList, "-gapTol");

	if(test!=-1){
		gapTol = test;
	}

	test = getOptionalArg(argList,"-nThreads");

	if(test!=-1){
		threads = test;
	}

	auto level = std::find(argList->begin(),argList->end(),"-seqLevel");

			if(level!=argList->end() && threads >1){

			seqLevel= true;
			
		}


	vector <string> * chromList = new vector<string>();

	if(!fs::exists(outputDir)){

		bool status = fs::create_directory(outputDir);

		if(!status){
			cerr<<"Output Directory could not be made"<<endl;
			throw std::exception();
		}
	}

	//Util::deleteFilesUnderDirectory(outputDir);


	//delete previous output files

	for (const auto & entry : fs::directory_iterator(outputDir)){

		string fileName = entry.path();

        cout << fileName <<endl;

		const string detectorEnding = "Detector.bed";

		if(fileName.length()>= detectorEnding.length() && !fileName.compare(fileName.length()-detectorEnding.length(),detectorEnding.length(),detectorEnding)){

			std::remove(fileName.c_str());
			cout<<"Deleted old output file: "<<fileName<<endl;

		}
	}
	

	cout << "Reading chromosome directory ..." << endl;
	Util::readChromList(chromDir, chromList,"fa");

	ChromosomeOneDigit *chrom;

	cout<<"Num threads: "<<threads<<endl;
	cout<<minDist<<endl;
	cout<<maxDist<<endl;
	cout<<minLenLTR<<endl;
	cout<<maxLenLTR<<endl;
	cout<<identity<<endl;
	cout<<k<<endl;
	cout<<minPlateauLength<<endl;
	cout<<gapTol<<endl;
	
	if(!seqLevel){
		#pragma omp parallel for schedule(dynamic) num_threads(threads)
		for (int i = 0; i < chromList->size(); i++)
		{

			string chromFile = chromList->at(i);

			int nameBegin = chromFile.find_last_of("/")+1;
			int nameEnd = chromFile.find_last_of(".") ;
			int len = nameEnd - nameBegin ;

			string name = chromFile.substr(nameBegin,len);

			string outFile = outputDir+"/"+name;
			ChromListMaker chromListMaker(chromFile);
			const vector<ChromosomeOneDigit *> *  chromList = chromListMaker.makeChromOneDigitList();
			for(auto chrom : *chromList){
				TrCollector * collector = new TrCollector(chrom, outFile,name, minDist, maxDist, minLenLTR, maxLenLTR,identity,k, minPlateauLength, gapTol,printCleanScores,printRawScores,bedOutput,displayNested);
				
				delete collector;
			}
		}
	}else{
		for (int i = 0; i < chromList->size(); i++)
		{

			string chromFile = chromList->at(i);

			int nameBegin = chromFile.find_last_of("/")+1;
			int nameEnd = chromFile.find_last_of(".") ;
			int len = nameEnd - nameBegin ;

			string name = chromFile.substr(nameBegin,len);

			string outFile = outputDir+"/"+name;
			ChromListMaker chromListMaker(chromFile);
			const vector<ChromosomeOneDigit *> *  chromList = chromListMaker.makeChromOneDigitList();
			#pragma omp parallel for schedule(dynamic) num_threads(threads)
			for(unsigned int j = 0; j < chromList->size(); j++){
				auto chrom = chromList->at(j);
				TrCollector * collector = new TrCollector(chrom, outFile,name, minDist, maxDist, minLenLTR, maxLenLTR,identity,k, minPlateauLength, gapTol,printCleanScores,printRawScores,bedOutput,displayNested);
				
				delete collector;
			}
	}

	chromList->clear();
	delete chromList;

	return 0;
	}
}

