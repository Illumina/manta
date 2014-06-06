// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Xiaoyu Chen
///


#include "assembly/IterativeAssembler.hh"
#include "blt_util/set_util.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"

#include <cassert>

#include <vector>
#include <algorithm>


// compile with this macro to get verbose output:
#define DEBUG_ASBL


// stream used by DEBUG_ASBL:
#ifdef DEBUG_ASBL
#include "blt_util/log.hh"
#include <iostream>

static
void print_unsignSet(const std::set<unsigned>& unsignSet)
{
    log_os << "[";
    BOOST_FOREACH(const unsigned us, unsignSet)
    {
        log_os << us << ",";
    }
    log_os << "]\n";
}

static
void print_stringSet(const std::set<std::string>& strSet)
{
    log_os << "[";
    BOOST_FOREACH(const std::string& str, strSet)
    {
        log_os << str << ",";
    }
    log_os << "]\n";
}
#endif


// maps kmers to positions in read
typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;
// maps kmers to support reads
typedef boost::unordered_map<std::string,std::set<unsigned> > str_set_uint_map_t;
typedef boost::unordered_map<std::string, std::pair<unsigned,unsigned> > str_pair_uint_map_t;



/**
 * Adds base @p base to the end (isEnd is true) or start (otherwise) of the contig.
 *
 *	@return The extended contig.
 */
static
std::string
addBase(
    const std::string& contig,
    const char base,
    const bool isEnd)
{
    if (isEnd) return contig + base;
    else       return base + contig;
}



/**
 * Returns a suffix (isEnd is true) or prefix (otherwise) of @p contig with length @p length.
 *
 *	@return The suffix or prefix.
 */
static
std::string
getEnd(
    const std::string& contig,
    const unsigned length,
    const bool isEnd)
{

    const unsigned csize(contig.size());
    assert(length <= csize);

    if (isEnd) return contig.substr((csize-length),length);
    else       return contig.substr(0,length);
}


#if 0
/// adapt Ole's function to represent word hash as a graph:
///
void
wordHashToDot(
    const str_uint_map_t& wordCount,
    std::ostream& os)
{
    static const unsigned MIN_KMER_FREQ(1);

    // kmer nodes with coverage higher than this get a different color
    //static const int lowCovGraphVisThreshold(3);
    static const unsigned lowCovGraphVisThreshold(MIN_KMER_FREQ);
    static const std::string lowCovNodeColor("red");
    static const std::string highCovNodeColor("green");

    os << "graph {\n";
    os << "node [ style = filled ];\n";
    str_uint_map_t aliasH;
    unsigned n(0);
    BOOST_FOREACH(const str_uint_map_t::value_type& val, wordCount)
    {
        const std::string& word(val.first);
        const unsigned cov(val.second);
        aliasH[word] = n++;
        const std::string& color(cov>lowCovGraphVisThreshold ? highCovNodeColor : lowCovNodeColor);
        os << word << "[label=\"cov" << cov << "\" color=" << color << "]\n";
    }

    // need to add edges here
    static const bool isEnd(true);
    BOOST_FOREACH(const str_uint_map_t::value_type& val, wordCount)
    {
        const std::string& word(val.first);
        const std::string tmp(getEnd(word,word.size()-1,isEnd));
        BOOST_FOREACH(const char symbol, alphabet)
        {
            const std::string newKey(addBase(tmp,symbol,isEnd));
            if (wordCount.find(newKey) != wordCount.end())
            {
                os << aliasH[word] << " -- " << aliasH[newKey] << ";\n";
            }
        }
    }
    os << "}\n";
}
#endif



/**
 * Extends the seed contig (aka most frequent k-mer)
 *
 */
static
bool
walk(const IterativeAssemblerOptions& opt,
     const std::string& seed,
     const unsigned wordLength,
     const str_uint_map_t& wordCount,
     const str_set_uint_map_t& wordReads,
     const std::set<std::string>& repeatWords,
     std::set<std::string>& unusedWords,
     AssembledContig& contig)
{
	const str_uint_map_t::const_iterator wordCountEnd(wordCount.end());
    const str_set_uint_map_t::const_iterator wordReadsEnd(wordReads.end());

    // we start with the seed
    str_set_uint_map_t::const_iterator wordReadsIter(wordReads.find(seed));
    assert(wordReadsIter != wordReadsEnd);
    contig.supportReads = wordReadsIter->second;
    contig.seq = seed;

    unusedWords.erase(seed);

    if (repeatWords.find(seed) != repeatWords.end())
    {
    #ifdef DEBUG_ASBL
    	log_os << "The seed is a repeat word " << seed << ". Stop walk.\n";
    #endif
    	return true;
    }

    bool isRepeatFound(false);

    // 0 => walk to the right, 1 => walk to the left
    for (unsigned mode(0); mode<2; ++mode)
    {
        const bool isEnd(mode==0);

        while (true)
        {
            const std::string tmp(getEnd(contig.seq, wordLength-1, isEnd));

#ifdef DEBUG_ASBL
            log_os << "# current contig : " << contig.seq << " size : " << contig.seq.size() << "\n"
                   << " getEnd : " << tmp << "\n";
            log_os << "contig rejecting reads : ";
            print_unsignSet(contig.rejectReads);
            log_os << "contig supporting reads : ";
            print_unsignSet(contig.supportReads);
#endif

            unsigned maxBaseCount(0);
            unsigned maxSharedReadCount(0);
            char maxBase(opt.alphabet[0]);
            std::set<unsigned> maxSharedReads;
            std::set<unsigned> supportReads2Remove;
            std::set<unsigned> supportReads2Add;
            std::set<unsigned> rejectReads2Add;

            BOOST_FOREACH(const char symbol, opt.alphabet)
            {
                const std::string newKey(addBase(tmp, symbol, isEnd));
#ifdef DEBUG_ASBL
                log_os << "Extending end : base " << symbol << " " << newKey << "\n";
#endif
                const str_uint_map_t::const_iterator wordCountIter(wordCount.find(newKey));
                if (wordCountIter == wordCountEnd) continue;
                const unsigned currWordCount(wordCountIter->second);

                wordReadsIter= wordReads.find(newKey);
                if (wordReadsIter == wordReadsEnd) continue;
                const std::set<unsigned>& currWordReads(wordReadsIter->second);

                // get the shared supporting reads between the contig and the current word
                std::set<unsigned> sharedReads;
                std::set_intersection(contig.supportReads.begin(), contig.supportReads.end(),
                                      currWordReads.begin(), currWordReads.end(),
                                      std::inserter(sharedReads, sharedReads.begin()));
#ifdef DEBUG_ASBL
                log_os << "Word supporting reads : ";
                print_unsignSet(currWordReads);
                log_os << "Contig-word shared reads : ";
                print_unsignSet(sharedReads);
#endif

                if (sharedReads.empty()) continue;

                const unsigned sharedReadCount(sharedReads.size());
                if (sharedReadCount > maxSharedReadCount)
                {
                    // the old shared reads support an unselected allele
                    // remove them from the contig's supporting reads
                    if (!maxSharedReads.empty())
                        supportReads2Remove.insert(maxSharedReads.begin(), maxSharedReads.end());
                    // the old supporting reads is for an unselected allele
                    // they become rejecting reads for the currently selected allele
                    if (!supportReads2Add.empty())
                        rejectReads2Add.insert(supportReads2Add.begin(), supportReads2Add.end());
                    // new supporting reads for the currently selected allele
                    supportReads2Add = currWordReads;

                    maxSharedReadCount = sharedReadCount;
                    maxSharedReads = sharedReads;
                    maxBaseCount = currWordCount;
                    maxBase = symbol;
                }
                else
                {
                    supportReads2Remove.insert(sharedReads.begin(), sharedReads.end());
                    rejectReads2Add.insert(currWordReads.begin(), currWordReads.end());
                }
            }
#ifdef DEBUG_ASBL
            log_os << "Winner is : " << maxBase << " with " << maxBaseCount << " occurrences." << "\n";
#endif

            if (maxBaseCount < opt.minCoverage)
            {
#ifdef DEBUG_ASBL
                log_os << "Coverage or error rate below threshold.\n"
                       << "maxBaseCount : " << maxBaseCount << " minCverage: " << opt.minCoverage << "\n";
#endif
                break;
            }


#ifdef DEBUG_ASBL
            log_os << "Adding base " << contig.seq << " " << maxBase << " " << mode << "\n";
#endif
            contig.seq = addBase(contig.seq, maxBase, isEnd);
#ifdef DEBUG_ASBL
            log_os << "New contig : " << contig.seq << "\n";
#endif

            // TODO: can add threshold for the count or percentage of shared reads
            {
#ifdef DEBUG_ASBL
                log_os << "Adding rejecting reads " << "\n"
                       << " Old : ";
                print_unsignSet(contig.rejectReads);
                log_os << " To be added : ";
                print_unsignSet(rejectReads2Add);
#endif
                // update rejecting reads
                // add reads that support the unselected allele
                BOOST_FOREACH(const unsigned rd, rejectReads2Add)
                {
                    contig.rejectReads.insert(rd);
                }
#ifdef DEBUG_ASBL
                log_os << " New : ";
                print_unsignSet(contig.rejectReads);
#endif

#ifdef DEBUG_ASBL
                log_os << "Updating supporting reads " << "\n"
                       << " Old : ";
                print_unsignSet(contig.supportReads);
                log_os << " To be added : ";
                print_unsignSet(supportReads2Add);
#endif
                // update supporting reads
                // add reads that support the selected allel
                BOOST_FOREACH(const unsigned rd, supportReads2Add)
                {
                    if (contig.rejectReads.find(rd) == contig.rejectReads.end())
                        contig.supportReads.insert(rd);
#ifdef DEBUG_ASBL
                    if (contig.rejectReads.find(rd) != contig.rejectReads.end())
                        log_os << "  Excluding rejected " << rd << "\n";
#endif
                }

#ifdef DEBUG_ASBL
                log_os << " To be removed : ";
                print_unsignSet(supportReads2Remove);
#endif
                // remove reads that do NOT support the selected allel anymore
                BOOST_FOREACH(const unsigned rd, supportReads2Remove)
                {
                    contig.supportReads.erase(rd);
                }
#ifdef DEBUG_ASBL
                log_os << " New : ";
                print_unsignSet(contig.supportReads);
#endif
            }

            // remove the last word from the unused list, so it cannot be used as the seed in finding the next contig
            const std::string lastWord(addBase(tmp, maxBase, isEnd));
            unusedWords.erase(lastWord);
            // stop walk in the current mode after seeing one repeat word
            if (repeatWords.find(lastWord) != repeatWords.end())
            {
#ifdef DEBUG_ASBL
                log_os << "Seen a repeat word " << lastWord << ". Stop walk in the current mode " << mode << "\n";
#endif
                isRepeatFound = true;
                break;
             }
        }

#ifdef DEBUG_ASBL
        log_os << "mode change. Current mode " << mode << "\n";
#endif
    }

    return isRepeatFound;
}



/// \params isFindRepeatReads if true record all reads with repeated words
///
static
void
getKmerCounts(
    const AssemblyReadInput& reads,
    const unsigned wordLength,
    str_uint_map_t& wordCount,
    str_set_uint_map_t& wordSupportReads)
{
    const unsigned readCount(reads.size());

    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        // stores the index of a kmer in a read sequence
        const std::string& seq(reads[readIndex]);
        const unsigned readLen(seq.size());

        // this read is unusable for assembly:
        if (readLen < wordLength) continue;

        // track all words from the read, including repetitive words
        std::set<std::string> readWords;
        for (unsigned j(0); j<=(readLen-wordLength); ++j)
        {
            const std::string word(seq.substr(j,wordLength));
            readWords.insert(word);
        }

        // total occurrences from this read
        BOOST_FOREACH(const std::string& word, readWords)
        {
            wordCount[word]++;
            // record the supporting read
            wordSupportReads[word].insert(readIndex);
        }
    }
}


static
unsigned
searchRepeats(
		const IterativeAssemblerOptions& opt,
		const unsigned index,
		const std::string& word,
		str_pair_uint_map_t& wordIndices,
		std::vector<std::string>& wordStack,
		std::set<std::string>& repeatWords)
{
	// set the depth index for the current word to the smallest unused index
	wordIndices[word] = std::pair<unsigned,unsigned>(index, index);
	unsigned nextIndex = index + 1;
	wordStack.push_back(word);

	const std::string tmp(getEnd(word, word.size()-1, true));
	BOOST_FOREACH(const char symbol, opt.alphabet)
	{
		// candidate successor of the current word
		const std::string nextWord(addBase(tmp, symbol, true));

		// homopolymer
		if (word == nextWord)
		{
			repeatWords.insert(word);
			continue;
		}

		// the successor word does not exist in the reads
		if (wordIndices.count(nextWord) == 0) continue;

		const unsigned nextWordIdx = wordIndices[nextWord].first;
		if (nextWordIdx == 0)
		{
			// the successor word has not been visited
			// recurse on it
			nextIndex = searchRepeats(opt, nextIndex, nextWord, wordIndices, wordStack, repeatWords);
			// update the current word's lowlink
			const unsigned wordLowLink = wordIndices[word].second;
			const unsigned nextWordLowLink = wordIndices[nextWord].second;
			wordIndices[word].second = std::min(wordLowLink, nextWordLowLink);
		}
		else
		{
			const bool isContained(std::find(wordStack.begin(), wordStack.end(), nextWord) != wordStack.end());
			if (isContained)
			{
				// the successor word is in stack and therefore in the current circle of words
				// only update the current word's lowlink
				const unsigned wordLowLink = wordIndices[word].second;
				wordIndices[word].second = std::min(wordLowLink, nextWordIdx);
			}
		}
	}

	// if the current word is a root node,
	if (wordIndices[word].second == index)
	{
		// exclude singletons
		bool isSingleton(wordStack.back() == word);
		if (isSingleton)
		{
			wordStack.pop_back();
		}
		// record identified repeat words (i.e. words in the current circle)
		else
		{
			while (true)
			{
				const std::string repeatWd = wordStack.back();
				repeatWords.insert(repeatWd);
				wordStack.pop_back();

				if (repeatWd == word) break;
			}
		}
	}

	return nextIndex;
}


static
void
getRepeatKmers(
		const IterativeAssemblerOptions& opt,
		const str_uint_map_t& wordCount,
		std::set<std::string>& repeatWords)
{
	str_pair_uint_map_t wordIndices;
	BOOST_FOREACH(const str_uint_map_t::value_type& wdct, wordCount)
	{
		wordIndices[wdct.first] = std::pair<unsigned,unsigned>(0, 0);
	}

	unsigned index = 1;
	std::vector<std::string> wordStack;
	BOOST_FOREACH(const str_pair_uint_map_t::value_type& wdidx, wordIndices)
	{
		const std::string word = wdidx.first;
		const unsigned wordIdx = wdidx.second.first;
		if (wordIdx == 0)
			index = searchRepeats(opt, index, word, wordIndices, wordStack, repeatWords);
	}
}


static
bool
buildContigs(
    const IterativeAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    const unsigned wordLength,
    Assembly& contigs)
{
#ifdef DEBUG_ASBL
    static const std::string logtag("buildContigs: ");
    log_os << logtag << "Building contigs with " << reads.size() << " reads.\n";
#endif

    contigs.clear();
    bool isAssemblySuccess(true);

    // counts the number of occurrences for each kmer in all reads
    str_uint_map_t wordCount;
    // records the supporting reads for each kmer
    str_set_uint_map_t wordSupportReads;
    // get counts and supporting reads for each kmer
    getKmerCounts(reads, wordLength, wordCount, wordSupportReads);

    // identify repeat kmers (i.e. circles from the de bruijn graph)
    std::set<std::string> repeatWords;
    getRepeatKmers(opt, wordCount, repeatWords);
#ifdef DEBUG_ASBL
    log_os << logtag << "Identified " << repeatWords.size() << " repeat words.\n";
    print_stringSet(repeatWords);
#endif

    // track kmers can be used as seeds for searching for the next contig
    std::set<std::string> unusedWords;
    BOOST_FOREACH(const str_uint_map_t::value_type& wdct, wordCount)
    {
    	// filter out kmers with too few coverage
    	if (wdct.second >= opt.minCoverage)
    		unusedWords.insert(wdct.first);
    }

    while (!unusedWords.empty())
    {
    	std::string maxWord;
    	unsigned maxWordCount(0);
    	// get the kmers corresponding the highest count
    	BOOST_FOREACH(const std::string& word, unusedWords)
    	{
    		assert (wordCount.count(word) > 0);
    		const unsigned currWordCount = wordCount.at(word);
    		if (currWordCount > maxWordCount)
    		{
    			maxWord = word;
    			maxWordCount = currWordCount;
    		}
    	}

    	// solve for a best contig in the graph by a heuristic greedy maxflow-ish criteria
    	AssembledContig contig;
    	bool isRepeatFound = walk(opt, maxWord, wordLength, wordCount, wordSupportReads, repeatWords, unusedWords, contig);
    	if (isRepeatFound) isAssemblySuccess = false;

#ifdef DEBUG_ASBL
    	log_os << logtag << "Found one contig of length " << contig.seq.size()
    		   << ", with supporting reads: ";
    	print_unsignSet(contig.supportReads);
    	log_os << ", with rejecting reads: ";
    	print_unsignSet(contig.rejectReads);
    	log_os << ". Contig seq: \n" << contig.seq << "\n";
#endif

    	contigs.push_back(contig);
    }

    // done with this now
    wordCount.clear();
    wordSupportReads.clear();

    return isAssemblySuccess;
}


static
void
selectContigs(
		const IterativeAssemblerOptions& opt,
		AssemblyReadOutput& readInfo,
		const unsigned normalReadCount,
		Assembly candidateContigs,
		Assembly& finalContigs)
{
#ifdef DEBUG_ASBL
	static const std::string logtag("selectContigs: ");
    log_os << logtag << "Start selecting contigs to be returned.\n";
#endif

    finalContigs.clear();
    unsigned finalContigCount(0);
    // a set of reads that has been used to construct contigs, including pseudo ones.
    std::set<unsigned> usedReads;
    // a set of pseudo reads that has been used to construct contigs
    std::set<unsigned> usedPseudoReads;

    while ((candidateContigs.size() > 0) && (finalContigCount < opt.maxAssemblyCount))
    {
    	// count unused reads that are not pseudo reads
    	const unsigned usedNormalReads = usedReads.size() - usedPseudoReads.size();
    	const unsigned unusedNormalReads = normalReadCount - usedNormalReads;
#ifdef DEBUG_ASBL
    	log_os << logtag << "# of candidateContigs: " << candidateContigs.size() << "\n";
    	log_os << logtag << "# of unused normal reads: " << unusedNormalReads << "\n";
#endif
        if (unusedNormalReads < opt.minUnusedReads) return;

        unsigned contigIndex(0);
        std::set<unsigned> contigs2Remove;

        AssembledContig selectedContig;
        unsigned selectedContigIndex;
        unsigned maxSupport(0);
        unsigned maxLength(0);
        BOOST_FOREACH(const AssembledContig& contig, candidateContigs)
        {
        	// identify new support reads that were not used for the previously identified contigs
        	std::set<unsigned> newSupportReads;
        	std::set_difference(contig.supportReads.begin(), contig.supportReads.end(),
        			            usedReads.begin(), usedReads.end(),
        			            std::inserter(newSupportReads, newSupportReads.end()));
#ifdef DEBUG_ASBL
        	log_os << logtag << "Contig #" << contigIndex << " newSupportReads=";
        	print_unsignSet(newSupportReads);
#endif

        	// count the number of new support reads that are not pseudo reads
        	unsigned newNormalSupport(0);
        	BOOST_FOREACH(const unsigned rd, newSupportReads)
        	{
        		const AssemblyReadInfo& rinfo(readInfo[rd]);
        		if (!rinfo.isPseudo) newNormalSupport++;
        	}
        	if (newNormalSupport < opt.minSupportReads)
        	{
#ifdef DEBUG_ASBL
        		log_os << logtag << "Contig #" << contigIndex << " to be skipped: too few non-pseudo support reads that has not been used for previously identified contigs.\n";
#endif
        		contigs2Remove.insert(contigIndex);
        		contigIndex++;
        		continue;
        	}

        	// either more support reads that were not used
        	// or the same number of supports but longer contig
        	const unsigned currNewSupport = newSupportReads.size();
        	const unsigned currContigLen = contig.seq.size();
        	bool isBetterContig((currNewSupport > maxSupport) ||
        			            ((currNewSupport == maxSupport) && (currContigLen > maxLength)));
        	if (isBetterContig)
        	{
        		selectedContig = contig;
        		selectedContigIndex = contigIndex;
        		maxSupport = currNewSupport;
        		maxLength = currContigLen;
        	}

        	contigIndex++;
        }


#ifdef DEBUG_ASBL
        log_os << logtag << "Contig #" << selectedContigIndex << " selected.\n";
#endif
    	// select one contig
    	finalContigs.push_back(selectedContig);

        contigs2Remove.insert(selectedContigIndex);
        // remove selected & failed contigs
        BOOST_REVERSE_FOREACH(const unsigned cix, contigs2Remove)
        	candidateContigs.erase(candidateContigs.begin()+cix);
#ifdef DEBUG_ASBL
    	log_os << logtag << "Removed " << contigs2Remove.size() << " contigs.\n";
#endif

        // update the info about used reads
        BOOST_FOREACH(const unsigned rd, selectedContig.supportReads)
        {
        	usedReads.insert(rd);
        	AssemblyReadInfo& rinfo(readInfo[rd]);
        	// read info record the ID of the very first contig that the read supports
        	// TODO: may need to record the IDs of all contigs that the read supports?
        	if (!rinfo.isUsed)
        	{
        		rinfo.isUsed = true;
        		rinfo.contigId = finalContigCount;
        	}
        	if (rinfo.isPseudo) usedPseudoReads.insert(rd);
        }
#ifdef DEBUG_ASBL
    	log_os << logtag << "Updated used reads: \n";
    	print_unsignSet(usedReads);
#endif

    	finalContigCount++;
    }
}


void
runIterativeAssembler(
    const IterativeAssemblerOptions& opt,
    AssemblyReadInput& reads,
    AssemblyReadOutput& readInfo,
    Assembly& contigs)
{
	const unsigned normalReadCount(reads.size());
#ifdef DEBUG_ASBL
    static const std::string logtag("runIterativeAssembler: ");
    log_os << logtag << "Starting assembly with " << normalReadCount << " reads.\n";
#endif
    assert(opt.alphabet.size()>1);

    readInfo.clear();
    readInfo.resize(reads.size());
    Assembly iterativeContigs;

    for (unsigned wordLength(opt.minWordLength); wordLength<=opt.maxWordLength; wordLength+=opt.wordStepSize)
    {
#ifdef DEBUG_ASBL
    	log_os << logtag << "Try " << wordLength << "-mer.\n";
#endif
    	const bool isAssemblySuccess = buildContigs(opt, reads, wordLength, iterativeContigs);
    	if (isAssemblySuccess)
    	{
#ifdef DEBUG_ASBL
    		log_os << logtag << "Assembly succeeded with " << wordLength << "-mer.\n";
#endif
    		break;
    	}

#ifdef DEBUG_ASBL
    	log_os << logtag << "Repeats encountered with " << wordLength << "-mer.\n";
#endif
    	// remove pseudo reads from the previous iteration
    	const unsigned readCount(reads.size());
    	for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    	{
    		AssemblyReadInfo& rinfo(readInfo[readIndex]);
    		if (rinfo.isPseudo)
    		{
    			reads.erase(reads.begin()+readIndex, reads.end());
    			readInfo.erase(readInfo.begin()+readIndex, readInfo.end());
#ifdef DEBUG_ASBL
    			log_os << logtag << "Removed " << (readCount - readIndex) << " pseudo reads (from the previous iteration).\n";
#endif
    			break;
    		}
    	}

    	unsigned addedCount(0);
    	//  Add contigs from the current iteration as pseudo reads
    	BOOST_FOREACH(const AssembledContig& contig, iterativeContigs)
    	{
    		if (contig.seq.size() > (wordLength+opt.wordStepSize))
    		{
#ifdef DEBUG_ASBL
    			log_os << logtag << "Adding a contig as pseudo read: " << contig.seq << ".\n";
#endif
    			reads.push_back(contig.seq);

    			AssemblyReadInfo rinfo;
    			rinfo.isPseudo = true;
    			readInfo.push_back(rinfo);

    			addedCount++;
    		}
    	}
#ifdef DEBUG_ASBL
    	log_os << logtag << "Added " << addedCount << " pseudo reads.\n";
#endif
    }

    // greedy selection of contigs to be returned
    selectContigs(opt, readInfo, normalReadCount, iterativeContigs, contigs);

#ifdef DEBUG_ASBL
    log_os << logtag << "Selected " << contigs.size() << "contigs.\n";
    unsigned index(1);
    BOOST_REVERSE_FOREACH(const AssembledContig& ctg, contigs)
    {
    	log_os << logtag <<"Selected contig # " << index << ": " << ctg.seq << "\n";
        log_os << logtag << "Contig supporting reads: ";
        print_unsignSet(ctg.supportReads);
        log_os << logtag << "Contig rejecting reads: ";
        print_unsignSet(ctg.rejectReads);
    	index++;
    }
#endif
}

