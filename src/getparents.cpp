/*
 *  The GetParents class is based on getparents.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 2/27/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "getparents.h"
#include "addtargets2.h"

// default - minchuck 64 chunk 4
/******************************************************************************/
vector<unsigned> GetParents::getChunkInfo(unsigned L, unsigned &Length)
	{
	vector<unsigned> Los;

	if (L <= opt->getMinchunk()) {
		Length = L;
		Los.push_back(0);
		return Los;
	}

	Length = (L - 1)/opt->getChunks() + 1;
	if (Length < opt->getMinchunk()) {
		Length = opt->getMinchunk();
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
	set<unsigned> TargetIndexes;

	unsigned queryLength = query.getSeqLength();

	SeqData QuerySD = query;

	unsigned ChunkLength;
	vector<unsigned> ChunkLos = getChunkInfo(queryLength, ChunkLength);

	// after breaking the query seqs into chunks look for parents for each chunk

	unsigned ChunkCount = SIZE(ChunkLos);
	PotentialParents parentCollector;
	for (unsigned ChunkIndex = 0; ChunkIndex < ChunkCount; ++ChunkIndex) {
		unsigned Lo = ChunkLos[ChunkIndex];

		string chunk = query.getSeq().substr(Lo);

		//QuerySD.Label = ChunkLabel.c_str();
		QuerySD.setSeq(chunk);

		parentCollector.addTargets(data, QuerySD, TargetIndexes);

		Lo += ChunkLength;
	}

	for (auto p = TargetIndexes.begin(); p != TargetIndexes.end(); ++p) {
		unsigned TargetIndex = *p;
		bool Accept = true;
		if (query.getAbund() > 0.0f) {
			string targetName = data->getName(TargetIndex);
			float targetAbund = data->getAbundance(TargetIndex);
			if (targetAbund > 0.0f && targetAbund < opt->getAbskew()*query.getAbund()) {
				Accept = false;
			}
		}

		if (Accept) {
			parents.push_back(TargetIndex);
		}
	}

	return parents;
}
/******************************************************************************/
