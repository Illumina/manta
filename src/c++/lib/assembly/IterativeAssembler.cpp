//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Xiaoyu Chen
/// \author Ole Schulz-Trieglaff
///

#include "assembly/IterativeAssembler.hpp"

#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <vector>

#include "boost/foreach.hpp"

#include "blt_util/set_util.hpp"

// compile with this macro to get verbose output:
//#define DEBUG_ASBL
//#define DEBUG_WALK

// stream used by DEBUG_ASBL:
#ifdef DEBUG_ASBL
#include <iostream>
#include "blt_util/log.hpp"

static void print_unsignSet(const std::set<unsigned>& unsignSet)
{
  log_os << "[";
  for (const unsigned us : unsignSet) {
    log_os << us << ",";
  }
  log_os << "]\n";
}

static void print_stringSet(const std::set<std::string>& strSet)
{
  log_os << "[";
  for (const std::string& str : strSet) {
    log_os << str << ",";
  }
  log_os << "]\n";
}
#endif

// maps kmers to positions in read
typedef std::unordered_map<std::string, unsigned> str_uint_map_t;
// maps kmers to supporting reads
typedef std::unordered_map<std::string, std::set<unsigned>>            str_set_uint_map_t;
typedef std::unordered_map<std::string, std::pair<unsigned, unsigned>> str_pair_uint_map_t;

/// Adds base @p base to the end (isEnd is true) or start (otherwise) of the contig.
///
/// \return the extended contig.
static std::string addBase(const std::string& contig, const char base, const bool isEnd)
{
  if (isEnd)
    return contig + base;
  else
    return base + contig;
}

/// Returns a suffix (isEnd is true) or prefix (otherwise) of the input contig with the specified length.
static std::string getEnd(const std::string& contig, const unsigned length, const bool isEnd)
{
  const unsigned csize(contig.size());
  assert(length <= csize);

  if (isEnd)
    return contig.substr((csize - length), length);
  else
    return contig.substr(0, length);
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
    for (const str_uint_map_t::value_type& val : wordCount)
    {
        const std::string& word(val.first);
        const unsigned cov(val.second);
        aliasH[word] = n++;
        const std::string& color(cov>lowCovGraphVisThreshold ? highCovNodeColor : lowCovNodeColor);
        os << word << "[label=\"cov" << cov << "\" color=" << color << "]\n";
    }

    // need to add edges here
    static const bool isEnd(true);
    for (const str_uint_map_t::value_type& val : wordCount)
    {
        const std::string& word(val.first);
        const std::string tmp(getEnd(word,word.size()-1,isEnd));
        for (const char symbol : alphabet)
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

/// Construct a contig from the seed provided
///
/// The contig is extened in both directions until
/// either the there is no sufficient support evidence
/// or the last k-mer is repeatitive (i.e. part of a bubble in the graph)
///
/// \return True if the contig runs into a repeatitive k-mer when extending in either mode
static bool walk(
    const IterativeAssemblerOptions& opt,
    const std::string&               seed,
    const unsigned                   wordLength,
    const str_uint_map_t&            wordCount,
    const str_set_uint_map_t&        wordReads,
    const std::set<std::string>&     repeatWords,
    std::set<std::string>&           unusedWords,
    AssembledContig&                 contig)
{
  const str_uint_map_t::const_iterator     wordCountEnd(wordCount.cend());
  const str_set_uint_map_t::const_iterator wordReadsEnd(wordReads.cend());

#ifdef DEBUG_WALK
  log_os << "\nSeed: " << seed << "\n";
#endif
  // we start with the seed
  str_set_uint_map_t::const_iterator wordReadsIter(wordReads.find(seed));
  assert(wordReadsIter != wordReadsEnd);
  contig.supportReads = wordReadsIter->second;
  contig.seq          = seed;
  unusedWords.erase(seed);

  if (repeatWords.find(seed) != repeatWords.end()) {
#ifdef DEBUG_WALK
    log_os << "The seed is a repeat word " << seed << ". Stop walk.\n";
#endif
    contig.conservativeRange.set_begin_pos(0);
    contig.conservativeRange.set_end_pos(wordLength);
    return true;
  }

  // collecting words used to build the contig
  std::set<std::string> wordsInContig;
  wordsInContig.insert(seed);

  const std::string tmpTrunk = getEnd(seed, wordLength - 1, false);
  // collecting rejecting reads for the seed from the unselected branches
  for (const char symbol : opt.alphabet) {
    // the seed itself
    if (symbol == seed[wordLength - 1]) continue;

    // add rejecting reads from an unselected word/branch
    const std::string newKey(addBase(tmpTrunk, symbol, true));
#ifdef DEBUG_WALK
    log_os << "Extending the seed trunk: base " << symbol << " " << newKey << "\n";
#endif

    wordReadsIter = wordReads.find(newKey);
    if (wordReadsIter == wordReadsEnd) continue;
    const std::set<unsigned>& unselectedReads(wordReadsIter->second);
#ifdef DEBUG_WALK
    log_os << "Supporting reads for the non-seed word : ";
    print_unsignSet(unselectedReads);
#endif

    contig.rejectReads.insert(unselectedReads.begin(), unselectedReads.end());
#ifdef DEBUG_WALK
    log_os << "seed's rejecting reads : ";
    print_unsignSet(contig.rejectReads);
#endif
  }

  bool isRepeatFound(false);

  // 0 => walk to the right, 1 => walk to the left
  for (unsigned mode(0); mode < 2; ++mode) {
    const bool isEnd(mode == 0);
    unsigned   conservativeEndOffset(0);

    while (true) {
      const std::string previousWord = getEnd(contig.seq, wordLength, isEnd);
      const std::string trunk(getEnd(contig.seq, wordLength - 1, isEnd));
#ifdef DEBUG_WALK
      log_os << "# current contig : " << contig.seq << " size : " << contig.seq.size() << "\n"
             << " getEnd : " << trunk << "\n";
      log_os << "contig rejecting reads : ";
      print_unsignSet(contig.rejectReads);
      log_os << "contig supporting reads : ";
      print_unsignSet(contig.supportReads);
#endif

      unsigned           maxBaseCount(0);
      unsigned           maxContigWordReadCount(0);
      char               maxBase(opt.alphabet[0]);
      std::string        maxWord;
      std::set<unsigned> maxWordReads;
      std::set<unsigned> maxContigWordReads;
      std::set<unsigned> previousWordReads;
      std::set<unsigned> supportReads2Remove;
      std::set<unsigned> rejectReads2Add;

      for (const char symbol : opt.alphabet) {
        const std::string newKey(addBase(trunk, symbol, isEnd));
#ifdef DEBUG_WALK
        log_os << "Extending end : base " << symbol << " " << newKey << "\n";
#endif
        const str_uint_map_t::const_iterator wordCountIter(wordCount.find(newKey));
        if (wordCountIter == wordCountEnd) continue;
        const unsigned currWordCount(wordCountIter->second);

        wordReadsIter = wordReads.find(newKey);
        if (wordReadsIter == wordReadsEnd) continue;
        const std::set<unsigned>& currWordReads(wordReadsIter->second);

        // get the shared supporting reads between the contig and the current word
        std::set<unsigned> contigWordReads;
        std::set_intersection(
            contig.supportReads.begin(),
            contig.supportReads.end(),
            currWordReads.begin(),
            currWordReads.end(),
            std::inserter(contigWordReads, contigWordReads.begin()));

        // get the shared supporting reads across two alleles
        std::set<unsigned> sharedReads;
        std::set_intersection(
            maxContigWordReads.begin(),
            maxContigWordReads.end(),
            currWordReads.begin(),
            currWordReads.end(),
            std::inserter(sharedReads, sharedReads.begin()));
#ifdef DEBUG_WALK
        log_os << "Word supporting reads : ";
        print_unsignSet(currWordReads);
        log_os << "Contig-word shared reads : ";
        print_unsignSet(contigWordReads);
        log_os << "Allele shared reads : ";
        print_unsignSet(sharedReads);
#endif

        if (contigWordReads.empty()) continue;

        const unsigned contigWordReadCount(contigWordReads.size());
        if (contigWordReadCount > maxContigWordReadCount) {
          // the old shared reads support an unselected allele if they don't support the new word
          // remove them from the contig's supporting reads
          if (!maxContigWordReads.empty()) {
            std::set<unsigned> toRemove;
            std::set_difference(
                maxContigWordReads.begin(),
                maxContigWordReads.end(),
                sharedReads.begin(),
                sharedReads.end(),
                std::inserter(toRemove, toRemove.begin()));
            if (!toRemove.empty()) supportReads2Remove.insert(toRemove.begin(), toRemove.end());
          }

          // the old supporting reads is for an unselected allele if they don't support the new word
          // they become rejecting reads for the currently selected allele
          if (!maxWordReads.empty()) {
            std::set<unsigned> toAdd;
            std::set_difference(
                maxWordReads.begin(),
                maxWordReads.end(),
                sharedReads.begin(),
                sharedReads.end(),
                std::inserter(toAdd, toAdd.begin()));
            if (!toAdd.empty()) rejectReads2Add.insert(toAdd.begin(), toAdd.end());
          }
          // new supporting reads for the currently selected allele
          maxWordReads = currWordReads;

          maxContigWordReadCount = contigWordReadCount;
          maxContigWordReads     = contigWordReads;
          maxBaseCount           = currWordCount;
          maxBase                = symbol;
          maxWord                = newKey;
        } else {
          std::set<unsigned> toRemove;
          std::set_difference(
              contigWordReads.begin(),
              contigWordReads.end(),
              sharedReads.begin(),
              sharedReads.end(),
              std::inserter(toRemove, toRemove.begin()));
          if (!toRemove.empty()) supportReads2Remove.insert(toRemove.begin(), toRemove.end());

          std::set<unsigned> toAdd;
          std::set_difference(
              currWordReads.begin(),
              currWordReads.end(),
              sharedReads.begin(),
              sharedReads.end(),
              std::inserter(toAdd, toAdd.begin()));
          if (!toAdd.empty()) rejectReads2Add.insert(toAdd.begin(), toAdd.end());
        }
      }

#ifdef DEBUG_WALK
      log_os << "Winner is : " << maxBase << " with " << maxBaseCount << " occurrences."
             << "\n";
#endif

      if (maxBaseCount < opt.minCoverage) {
#ifdef DEBUG_WALK
        log_os << "Coverage or error rate below threshold.\n"
               << "maxBaseCount : " << maxBaseCount << " minCoverage: " << opt.minCoverage << "\n";
#endif
        break;
      }

      // stop walk in the current mode after seeing one repeat word
      if (wordsInContig.find(maxWord) != wordsInContig.end()) {
#ifdef DEBUG_WALK
        log_os << "Seen a repeat word " << maxWord << ".\n Stop walk in the current mode " << mode << "\n";
#endif
        isRepeatFound = true;
        break;
      }

#ifdef DEBUG_WALK
      log_os << "Adding base " << contig.seq << " " << maxBase << " " << mode << "\n";
#endif
      contig.seq = addBase(contig.seq, maxBase, isEnd);
#ifdef DEBUG_WALK
      log_os << "New contig : " << contig.seq << "\n";
#endif

      if ((conservativeEndOffset != 0) || (maxBaseCount < opt.minConservativeCoverage))
        conservativeEndOffset += 1;
#ifdef DEBUG_WALK
      log_os << "conservative end offset : " << conservativeEndOffset << "\n";
#endif

      // TODO: can add threshold for the count or percentage of shared reads
      {
        // walk backwards for one step at a branching point
        if (maxWordReads != previousWordReads) {
          const char tmpSymbol = (isEnd ? previousWord[0] : previousWord[wordLength - 1]);
          for (const char symbol : opt.alphabet) {
            // the selected branch: skip the backward word itself
            if (symbol == tmpSymbol) continue;

            // add rejecting reads from an unselected branch
            const std::string newKey(addBase(trunk, symbol, !isEnd));
#ifdef DEBUG_WALK
            log_os << "Extending end backwards: base " << symbol << " " << newKey << "\n";
#endif
            // the selected branch: skip the word just extended
            if (newKey == maxWord) continue;

            wordReadsIter = wordReads.find(newKey);
            if (wordReadsIter == wordReadsEnd) continue;

            const std::set<unsigned>& backWordReads(wordReadsIter->second);
#ifdef DEBUG_WALK
            log_os << "Supporting reads for the backwards word : ";
            print_unsignSet(backWordReads);
#endif
            // get the shared supporting reads across two alleles
            std::set<unsigned> sharedReads;
            std::set_intersection(
                maxContigWordReads.begin(),
                maxContigWordReads.end(),
                backWordReads.begin(),
                backWordReads.end(),
                std::inserter(sharedReads, sharedReads.begin()));

            std::set<unsigned> toUpdate;
            std::set_difference(
                backWordReads.begin(),
                backWordReads.end(),
                sharedReads.begin(),
                sharedReads.end(),
                std::inserter(toUpdate, toUpdate.begin()));
            if (!toUpdate.empty()) {
              rejectReads2Add.insert(toUpdate.begin(), toUpdate.end());
              supportReads2Remove.insert(toUpdate.begin(), toUpdate.end());
            }

#ifdef DEBUG_WALK
            log_os << "rejectReads2Add upated : ";
            print_unsignSet(rejectReads2Add);
            log_os << "supportReads2Remove updated : ";
            print_unsignSet(supportReads2Remove);
#endif
          }
        }
        previousWordReads = maxWordReads;

#ifdef DEBUG_WALK
        log_os << "Adding rejecting reads "
               << "\n"
               << " Old : ";
        print_unsignSet(contig.rejectReads);
        log_os << " To be added : ";
        print_unsignSet(rejectReads2Add);
#endif
        // update rejecting reads
        // add reads that support the unselected allele
        for (const unsigned rd : rejectReads2Add) {
          contig.rejectReads.insert(rd);
        }
#ifdef DEBUG_WALK
        log_os << " New : ";
        print_unsignSet(contig.rejectReads);
#endif

#ifdef DEBUG_WALK
        log_os << "Updating supporting reads "
               << "\n"
               << " Old : ";
        print_unsignSet(contig.supportReads);
        log_os << " To be added : ";
        print_unsignSet(maxWordReads);
#endif
        // update supporting reads
        // add reads that support the selected allel
        for (const unsigned rd : maxWordReads) {
          if (contig.rejectReads.find(rd) == contig.rejectReads.end()) contig.supportReads.insert(rd);
#ifdef DEBUG_WALK
          if (contig.rejectReads.find(rd) != contig.rejectReads.end())
            log_os << "  Excluding rejected read " << rd << " for the contig\n";
#endif
        }

#ifdef DEBUG_WALK
        log_os << " To be removed : ";
        print_unsignSet(supportReads2Remove);
#endif
        // remove reads that do NOT support the selected allel anymore
        for (const unsigned rd : supportReads2Remove) {
          contig.supportReads.erase(rd);
        }
#ifdef DEBUG_WALK
        log_os << " New : ";
        print_unsignSet(contig.supportReads);
#endif
      }

      // remove the last word from the unused list, so it cannot be used as the seed in finding the next
      // contig
      unusedWords.erase(maxWord);
      // collect the words used to build the contig
      wordsInContig.insert(maxWord);
    }

    // set conservative coverage range for the contig
    if (mode == 0)
      contig.conservativeRange.set_end_pos(conservativeEndOffset);
    else
      contig.conservativeRange.set_begin_pos(conservativeEndOffset);

#ifdef DEBUG_WALK
    log_os << "mode change. Current mode " << mode << "\n";
#endif
  }

  contig.conservativeRange.set_end_pos(contig.seq.size() - contig.conservativeRange.end_pos());

  return isRepeatFound;
}

/// Construct k-mer maps
/// k-mer ==> number of reads containing the k-mer
/// k-mer ==> a list of read IDs containg the k-mer
static void getKmerCounts(
    const IterativeAssemblerOptions& opt,
    const AssemblyReadInput&         reads,
    AssemblyReadOutput&              readInfo,
    const unsigned                   wordLength,
    str_uint_map_t&                  wordCount,
    str_set_uint_map_t&              wordSupportReads)
{
  const unsigned readCount(reads.size());

  for (unsigned readIndex(0); readIndex < readCount; ++readIndex) {
    // stores the index of a kmer in a read sequence
    const std::string& seq(reads[readIndex]);
    const unsigned     readLen(seq.size());

    // this read is unusable for assembly:
    if (readLen < wordLength) continue;

    // track all words from the read, including repetitive words
    std::set<std::string> readWords;
    for (unsigned j(0); j <= (readLen - wordLength); ++j) {
      const std::string word(seq.substr(j, wordLength));

      // filter words with "N" (either directly from input alignment
      // or marked due to low basecall quality:
      if (word.find('N') != std::string::npos) continue;

      readWords.insert(word);
    }

    AssemblyReadInfo& rinfo(readInfo[readIndex]);
    unsigned          wordCountAdd = 1;
    // pseudo reads must have passed coverage check with smaller kmers
    // Assigning minCoverage (instead of 1) to a pseudo read allows the pseudo read to rescue the regions
    // where coverage is as low as minCoverage and reads overlap is small
    if (rinfo.isPseudo) wordCountAdd = opt.minCoverage;

    // total occurrences from this read
    for (const std::string& word : readWords) {
      wordCount[word] += wordCountAdd;
      // record the supporting read
      wordSupportReads[word].insert(readIndex);
    }
  }
}

/// Identify repetitive k-mers
/// i.e. k-mers that form a circular subgraph
///
static unsigned searchRepeats(
    const IterativeAssemblerOptions& opt,
    const unsigned                   index,
    const std::string&               word,
    str_pair_uint_map_t&             wordIndices,
    std::vector<std::string>&        wordStack,
    std::set<std::string>&           repeatWords)
{
  // set the depth index for the current word to the smallest unused index
  wordIndices[word]  = std::pair<unsigned, unsigned>(index, index);
  unsigned nextIndex = index + 1;
  wordStack.push_back(word);

  const std::string tmp(getEnd(word, word.size() - 1, true));
  for (const char symbol : opt.alphabet) {
    // candidate successor of the current word
    const std::string nextWord(addBase(tmp, symbol, true));

    // homopolymer
    if (word == nextWord) {
      repeatWords.insert(word);
      continue;
    }

    // the successor word does not exist in the reads
    if (wordIndices.count(nextWord) == 0) continue;

    const unsigned nextWordIdx = wordIndices[nextWord].first;
    if (nextWordIdx == 0) {
      // the successor word has not been visited
      // recurse on it
      nextIndex = searchRepeats(opt, nextIndex, nextWord, wordIndices, wordStack, repeatWords);
      // update the current word's lowlink
      const unsigned wordLowLink     = wordIndices[word].second;
      const unsigned nextWordLowLink = wordIndices[nextWord].second;
      wordIndices[word].second       = std::min(wordLowLink, nextWordLowLink);
    } else {
      const bool isContained(std::find(wordStack.begin(), wordStack.end(), nextWord) != wordStack.end());
      if (isContained) {
        // the successor word is in stack and therefore in the current circle of words
        // only update the current word's lowlink
        const unsigned wordLowLink = wordIndices[word].second;
        wordIndices[word].second   = std::min(wordLowLink, nextWordIdx);
      }
    }
  }

  // if the current word is a root node,
  if (wordIndices[word].second == index) {
    const std::string& lastWord(wordStack.back());
    // exclude singletons
    bool isSingleton(lastWord == word);
    if (isSingleton) {
      wordStack.pop_back();
    } else {
      // record identified repeat words (i.e. words in the current circle) if the circle is small

      const unsigned lastWordIndex(wordIndices[lastWord].first);
      const bool     isSmallCircle((lastWordIndex - index) <= 50);
      while (true) {
        const std::string repeatWd = wordStack.back();
        if (isSmallCircle) repeatWords.insert(repeatWd);
        wordStack.pop_back();

        if (repeatWd == word) break;
      }
    }
  }

  return nextIndex;
}

static void getRepeatKmers(
    const IterativeAssemblerOptions& opt, const str_uint_map_t& wordCount, std::set<std::string>& repeatWords)
{
  str_pair_uint_map_t wordIndices;
  for (const str_uint_map_t::value_type& wdct : wordCount) {
    wordIndices[wdct.first] = std::pair<unsigned, unsigned>(0, 0);
  }

  unsigned                 index = 1;
  std::vector<std::string> wordStack;
  for (const str_pair_uint_map_t::value_type& wdidx : wordIndices) {
    const std::string word    = wdidx.first;
    const unsigned    wordIdx = wdidx.second.first;
    if (wordIdx == 0) index = searchRepeats(opt, index, word, wordIndices, wordStack, repeatWords);
  }
}

static bool buildContigs(
    const IterativeAssemblerOptions& opt,
    const AssemblyReadInput&         reads,
    AssemblyReadOutput&              readInfo,
    const unsigned                   wordLength,
    Assembly&                        contigs)
{
#ifdef DEBUG_ASBL
  static const std::string logtag("buildContigs: ");
  log_os << logtag << "Building contigs with " << reads.size() << " reads.\n";
  for (unsigned readIndex(0); readIndex < reads.size(); readIndex++) {
    log_os << "read #" << readIndex << ": " << reads[readIndex] << "\n";
  }
#endif

  contigs.clear();
  bool isAssemblySuccess(true);

  // counts the number of occurrences for each kmer in all reads
  str_uint_map_t wordCount;
  // records the supporting reads for each kmer
  str_set_uint_map_t wordSupportReads;
  // get counts and supporting reads for each kmer
  getKmerCounts(opt, reads, readInfo, wordLength, wordCount, wordSupportReads);

  // identify repeat kmers (i.e. circles from the de bruijn graph)
  std::set<std::string> repeatWords;
  getRepeatKmers(opt, wordCount, repeatWords);
#ifdef DEBUG_ASBL
  log_os << logtag << "Identified " << repeatWords.size() << " repeat words.\n";
  print_stringSet(repeatWords);
#endif

  // track kmers can be used as seeds for searching for the next contig
  std::set<std::string> unusedWords;
  for (const str_uint_map_t::value_type& wdct : wordCount) {
    // filter out kmers with too few coverage
    if (wdct.second >= opt.minCoverage) unusedWords.insert(wdct.first);
  }

  // limit the number of contigs generated for the seek of speed
  while ((!unusedWords.empty()) && (contigs.size() < 2 * opt.maxAssemblyCount)) {
    std::string maxWord;
    unsigned    maxWordCount(0);
    // get the kmers corresponding the highest count
    for (const std::string& word : unusedWords) {
      assert(wordCount.count(word) > 0);
      const unsigned currWordCount = wordCount.at(word);
      if (currWordCount > maxWordCount) {
        maxWord      = word;
        maxWordCount = currWordCount;
      }
    }

    // solve for a best contig in the graph by a heuristic greedy maxflow-ish criteria
    AssembledContig contig;
    bool            isRepeatFound =
        walk(opt, maxWord, wordLength, wordCount, wordSupportReads, repeatWords, unusedWords, contig);
    if (isRepeatFound) isAssemblySuccess = false;

#ifdef DEBUG_ASBL
    log_os << logtag << "Found one contig of length " << contig.seq.size() << ", with supporting reads: ";
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

static void selectContigs(
    const IterativeAssemblerOptions& opt,
    AssemblyReadOutput&              readInfo,
    const unsigned                   normalReadCount,
    Assembly                         candidateContigs,
    Assembly&                        finalContigs)
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

  // contig are selected based on the number of supporting reads that are not pseudo
  while ((!candidateContigs.empty()) && (finalContigCount < opt.maxAssemblyCount)) {
    // count unused reads that are not pseudo reads
    const unsigned usedNormalReads   = usedReads.size() - usedPseudoReads.size();
    const unsigned unusedNormalReads = normalReadCount - usedNormalReads;
#ifdef DEBUG_ASBL
    log_os << logtag << "# of candidateContigs: " << candidateContigs.size() << "\n";
    log_os << logtag << "# of unused normal reads: " << unusedNormalReads << "\n";
#endif
    if (unusedNormalReads < opt.minUnusedReads) return;

    unsigned           contigIndex(0);
    std::set<unsigned> contigs2Remove;

    AssembledContig selectedContig;
    unsigned        selectedContigIndex;
    unsigned        maxSupport(0);
    unsigned        maxLength(0);
    for (const AssembledContig& contig : candidateContigs) {
      // identify new supporting reads that were not used for the previously identified contigs
      std::set<unsigned> newSupportReads;
      std::set_difference(
          contig.supportReads.begin(),
          contig.supportReads.end(),
          usedReads.begin(),
          usedReads.end(),
          std::inserter(newSupportReads, newSupportReads.end()));
#ifdef DEBUG_ASBL
      log_os << logtag << "Contig #" << contigIndex << " newSupportReads=";
      print_unsignSet(newSupportReads);
#endif

      // count the number of new supporting reads that are not pseudo reads
      unsigned newNormalSupport(0);
      for (const unsigned rd : newSupportReads) {
        const AssemblyReadInfo& rinfo(readInfo[rd]);
        if (!rinfo.isPseudo) newNormalSupport++;
      }
      if (newNormalSupport < opt.minSupportReads) {
#ifdef DEBUG_ASBL
        log_os
            << logtag << "Contig #" << contigIndex
            << " to be skipped: too few non-pseudo supporting reads that has not been used for previously identified contigs.\n";
#endif
        contigs2Remove.insert(contigIndex);
        contigIndex++;
        continue;
      }

      // either more supporting reads that were not used
      // or the same number of supports but longer contig
      const unsigned currNewSupport = newSupportReads.size();
      const unsigned currContigLen  = contig.seq.size();
      const bool     isBetterContig(
          (currNewSupport > maxSupport) || ((currNewSupport == maxSupport) && (currContigLen > maxLength)));
      if (isBetterContig) {
        selectedContig      = contig;
        selectedContigIndex = contigIndex;
        maxSupport          = currNewSupport;
        maxLength           = currContigLen;
      }

      contigIndex++;
    }

    // no more contigs selected, selection is done.
    if (maxSupport == 0) break;

#ifdef DEBUG_ASBL
    log_os << logtag << "Contig #" << selectedContigIndex << " selected.\n";
#endif
    // select one contig
    finalContigs.push_back(selectedContig);

    contigs2Remove.insert(selectedContigIndex);
    // remove selected & failed contigs
    BOOST_REVERSE_FOREACH(const unsigned cix, contigs2Remove)
    {
      candidateContigs.erase(candidateContigs.begin() + cix);
    }
#ifdef DEBUG_ASBL
    log_os << logtag << "Removed " << contigs2Remove.size() << " contigs.\n";
#endif

    // update the info about used reads
    for (const unsigned rd : selectedContig.supportReads) {
      usedReads.insert(rd);
      AssemblyReadInfo& rinfo(readInfo[rd]);
      // read info record the ID of all contigs that the read supports
      if (!rinfo.isUsed) rinfo.isUsed = true;
      rinfo.contigIds.push_back(finalContigCount);

      if (rinfo.isPseudo) usedPseudoReads.insert(rd);
    }
#ifdef DEBUG_ASBL
    log_os << logtag << "Updated used reads: \n";
    print_unsignSet(usedReads);
#endif

    finalContigCount++;
  }
}

void runIterativeAssembler(
    const IterativeAssemblerOptions& opt,
    AssemblyReadInput&               reads,
    AssemblyReadOutput&              readInfo,
    Assembly&                        contigs)
{
  const unsigned normalReadCount(reads.size());
#ifdef DEBUG_ASBL
  static const std::string logtag("runIterativeAssembler: ");
  log_os << logtag << "Starting assembly with " << normalReadCount << " reads.\n";
#endif

  assert(opt.alphabet.size() > 1);

  readInfo.clear();
  readInfo.resize(reads.size());
  Assembly iterativeContigs;

  for (unsigned wordLength(opt.minWordLength); wordLength <= opt.maxWordLength;
       wordLength += opt.wordStepSize) {
#ifdef DEBUG_ASBL
    log_os << logtag << "Try " << wordLength << "-mer.\n";
#endif
    const bool isAssemblySuccess = buildContigs(opt, reads, readInfo, wordLength, iterativeContigs);

    // remove pseudo reads from the previous iteration
    const unsigned readCount(reads.size());

    if (isAssemblySuccess) {
#ifdef DEBUG_ASBL
      log_os << logtag << "Assembly succeeded with " << wordLength << "-mer.\n";
#endif
      break;
    }

#ifdef DEBUG_ASBL
    log_os << logtag << "Repeats encountered with " << wordLength << "-mer.\n";
#endif
    for (unsigned readIndex(0); readIndex < readCount; ++readIndex) {
      AssemblyReadInfo& rinfo(readInfo[readIndex]);
      if (rinfo.isPseudo) {
        reads.erase(reads.begin() + readIndex, reads.end());
        readInfo.erase(readInfo.begin() + readIndex, readInfo.end());
#ifdef DEBUG_ASBL
        log_os << logtag << "Removed " << (readCount - readIndex)
               << " pseudo reads (from the previous iteration).\n";
#endif
        break;
      }
    }

    unsigned addedCount(0);
    //  Add contigs from the current iteration as pseudo reads
    for (const AssembledContig& contig : iterativeContigs) {
      if (contig.seq.size() > (wordLength + opt.wordStepSize)) {
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
  log_os << logtag << "Selected " << contigs.size() << " contigs.\n";
  unsigned index(1);
  for (const AssembledContig& ctg : contigs) {
    log_os << logtag << "Selected contig # " << index << ": " << ctg.seq << "\n";
    log_os << logtag << "Contig supporting reads: ";
    print_unsignSet(ctg.supportReads);
    log_os << logtag << "Contig rejecting reads: ";
    print_unsignSet(ctg.rejectReads);
    index++;
  }
#endif
}
