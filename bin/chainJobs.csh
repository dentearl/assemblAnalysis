#!/bin/csh -fe
zcat psl/ce9.2bit:chrV:0-10010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrV:0-10010000.chain
zcat psl/ce9.2bit:chrV:10000000-20010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrV:10000000-20010000.chain
zcat psl/ce9.2bit:chrV:20000000-20924143.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrV:20000000-20924143.chain
zcat psl/ce9.2bit:chrX:0-10010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrX:0-10010000.chain
zcat psl/ce9.2bit:chrX:10000000-17718854.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrX:10000000-17718854.chain
zcat psl/ce9.2bit:chrIV:0-10010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrIV:0-10010000.chain
zcat psl/ce9.2bit:chrIV:10000000-17493784.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrIV:10000000-17493784.chain
zcat psl/ce9.2bit:chrII:0-10010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrII:0-10010000.chain
zcat psl/ce9.2bit:chrII:10000000-15279323.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrII:10000000-15279323.chain
zcat psl/ce9.2bit:chrI:0-10010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrI:0-10010000.chain
zcat psl/ce9.2bit:chrI:10000000-15072421.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrI:10000000-15072421.chain
zcat psl/ce9.2bit:chrIII:0-10010000.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrIII:0-10010000.chain
zcat psl/ce9.2bit:chrIII:10000000-13783685.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrIII:10000000-13783685.chain
zcat psl/ce9.2bit:chrM:0-13794.*.psl.gz \
    | axtChain -psl -verbose=0 -minScore=5000 -linearGap=medium \
	stdin /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdout \
   | chainAntiRepeat /hive/data/genomes/ce9/bed/testLastz/ce9.2bit /hive/data/genomes/ce9/bed/testLastz/cb3.2bit stdin chain/ce9.2bit:chrM:0-13794.chain
find ./chain -name "*.chain" | chainMergeSort -inputList=stdin | gzip -c > ce9.cb3.all.chain.gz
