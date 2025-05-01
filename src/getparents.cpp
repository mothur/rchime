/*
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "getparents.h"

// default - minchuck 64 chunk 4
/******************************************************************************/
GetParents::GetParents(Options o) {
	minChunk = o.getMinchunk();
	numChunks = o.getChunks();
	abskew = o.getAbskew();
}
/******************************************************************************/
vector<unsigned> GetParents::getChunkInfo(unsigned L, unsigned &Length) {
	vector<unsigned> Los;

	if (L <= minChunk) {
		Length = L;
		Los.push_back(0);
		return Los;
	}

	Length = (L - 1)/numChunks + 1;
	if (Length < minChunk) {
		Length = minChunk;
	}

	unsigned Lo = 0;
	for (;;) {
		if (Lo + Length >= L) {
			Lo = L - Length - 1;
			Los.push_back(Lo);
			return Los;
		}
		Los.push_back(Lo);
		Lo += Length;
	}
}
/******************************************************************************/
vector<unsigned> GetParents::getCandidateParents(SeqDB* data, const SeqData & query) {

	vector<unsigned> parents;
	set<unsigned> ParentIndexes;

	unsigned queryLength = query.getSeqLength();
	unsigned ChunkLength;
	vector<unsigned> ChunkLos = getChunkInfo(queryLength, ChunkLength);

	// after breaking the query seqs into chunks look for parents for each chunk
	unsigned ChunkCount = SIZE(ChunkLos);
	for (unsigned ChunkIndex = 0; ChunkIndex < ChunkCount; ++ChunkIndex) {
		unsigned Lo = ChunkLos[ChunkIndex];

		string chunk = query.getSeq().substr(Lo, ChunkLength);
		Lo += ChunkLength;

		data->addPotentialParents(chunk, ParentIndexes);
	}

	// add all potential parents that meet abundance requirements
	float cutoff = abskew*query.getAbund();
	for (auto p = ParentIndexes.begin(); p != ParentIndexes.end(); ++p) {
		unsigned ParentIndex = *p;

		string targetName = data->getName(ParentIndex);
		float targetAbund = data->getAbundance(ParentIndex);

		if (data->isDenovo()) {
			if (!(targetAbund < cutoff)) {
				parents.push_back(ParentIndex);
			}
		}else {
			parents.push_back(ParentIndex);
		}
	}
	return parents;
}
/******************************************************************************/


