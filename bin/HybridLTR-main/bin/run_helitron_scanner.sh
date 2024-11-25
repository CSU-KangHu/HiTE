#!/bin/bash

GENOME=$1
GENOMEFA=$2
HSDIR=$3
HSJAR=$4

#### the base path to find all the lcv alignments
#HSDIR=/public/home/hpc194701009/repeat_detect_tools
### where to find HelitronScanner.jar
#HSJAR=/public/home/hpc194701009/repeat_detect_tools/HelitronScanner/HelitronScanner.jar


### helitron scanner needs some memory to load each chromosome in, so remember that when picking a queue
MEMGB=2

###########################
##   DIRECT ORIENTATION  ##
###########################

##find helitron heads
### will load each chromosome into memory, without splitting into 1Mb batches (-buffer_size option ==0)
java -Xmx${MEMGB}g -jar ${HSJAR} scanHead -lcv_filepath ${HSDIR}/head.lcvs -g $GENOMEFA -buffer_size 0 -output ${GENOME}.HelitronScanner.head
## helitron tails
java -Xmx${MEMGB}g -jar ${HSJAR} scanTail -lcv_filepath ${HSDIR}/tail.lcvs -g $GENOMEFA -buffer_size 0 -output ${GENOME}.HelitronScanner.tail

## pair the ends to generate possible helitrons
java -Xmx${MEMGB}g -jar ${HSJAR} pairends -head_score ${GENOME}.HelitronScanner.head -tail_score ${GENOME}.HelitronScanner.tail -output ${GENOME}.HelitronScanner.pairends

## draw the helitrons into fastas
java -Xmx${MEMGB}g -jar ${HSJAR} draw -pscore ${GENOME}.HelitronScanner.pairends -g $GENOMEFA -output ${GENOME}.HelitronScanner.draw -pure_helitron

############################
##    REVERSE COMPLEMENT  ##
############################

##find helitron heads
### will load each chromosome into memory, without splitting into 1Mb batches (-buffer_size option ==0)
java -Xmx${MEMGB}g -jar ${HSJAR} scanHead -lcv_filepath ${HSDIR}/head.lcvs -g $GENOMEFA -buffer_size 0 --rc -output ${GENOME}.HelitronScanner.rc.head
## helitron tails
java -Xmx${MEMGB}g -jar ${HSJAR} scanTail -lcv_filepath ${HSDIR}/tail.lcvs -g $GENOMEFA -buffer_size 0 --rc -output ${GENOME}.HelitronScanner.rc.tail

## pair the ends to generate possible helitrons
java -Xmx${MEMGB}g -jar ${HSJAR} pairends -head_score ${GENOME}.HelitronScanner.rc.head -tail_score ${GENOME}.HelitronScanner.rc.tail --rc -output ${GENOME}.HelitronScanner.rc.pairends

## draw the helitrons
java -Xmx${MEMGB}g -jar ${HSJAR} draw -pscore ${GENOME}.HelitronScanner.rc.pairends -g $GENOMEFA -output ${GENOME}.HelitronScanner.draw.rc -pure_helitron