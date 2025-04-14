
#include "searchchime.h"
#include "alignchime.h"
#include "globalalign2.h"

/******************************************************************************/
SearchChime::SearchChime(Options* o) { 
	opt = o; 
	chimeAlign = new AlignChimes(opt->getXn(), opt->getDn(), opt->getXa(),
								 opt->getSkipgaps(), opt->getSkipgaps2());
}
/******************************************************************************/
SearchChime::~SearchChime() { delete chimeAlign; }
/******************************************************************************/
vector<unsigned> SearchChime::getSmoothedIdVec(const SeqData &queryData, const SeqData &parentData,
 									const string &Path, unsigned d) {
	vector<unsigned> IdVec;
	const unsigned ColCount = SIZE(Path);

	string query = queryData.getSeq();
	string parent = parentData.getSeq();

	const unsigned QL = query.length();

	if (QL <= d) {
		IdVec.resize(QL, 0);
		return IdVec;
	}

	unsigned QPos = 0;
	unsigned PPos = 0;

	vector<bool> SameVec;
	SameVec.reserve(QL);
	for (unsigned Col = 0; Col < ColCount; ++Col) {
		char c = Path[Col];

		bool Same = false;
		if (c == 'M') {
		    Byte q = query[QPos];
		    Byte p = parent[PPos];
			Same = (toupper(q) == toupper(p));
		}

		if (c == 'M' || c == 'D') {
			++QPos;
			SameVec.push_back(Same);
		}

		if (c == 'M' || c == 'I') {
			++PPos;
		}
	}
	unsigned n = 0;
	for (unsigned QPos = 0; QPos < d; ++QPos) {
		if (SameVec[QPos]) {
			++n;
		}
		IdVec.push_back(n);
	}

	for (unsigned QPos = d; QPos < QL; ++QPos) {
		if (SameVec[QPos]) {
			++n;
		}
		IdVec.push_back(n);
		if (SameVec[QPos-d]) {
			--n;
		}
	}

	return IdVec;
}
/******************************************************************************/
bool SearchChime::searchChime(SeqDB* database, const SeqData & queryData, 
  								ChimeHit2 &Hit) {

	//float MinFractId = 0.95f;	

	Hit.QLabel = queryData.getName();
	Hit.AbQ = queryData.getAbund();

	GetParents parentFinder(opt);
	vector<unsigned> Parents = parentFinder.getCandidateParents(database, queryData);

	unsigned ParentCount = SIZE(Parents);
	if (ParentCount <= 1) {
		return false;
	}

    // set in getChunkInfo
	unsigned ChunkLength;
	vector<unsigned> ChunkLos = parentFinder.getChunkInfo(queryData.getSeqLength(), ChunkLength);
	const unsigned ChunkCount = SIZE(ChunkLos);

	vector<unsigned> ChunkIndexToBestId(ChunkCount, 0);
	vector<unsigned> ChunkIndexToBestParentIndex(ChunkCount, UINT_MAX);

	vector<SeqData> PSDs;
	vector<string> Paths;
	double TopPctId = 0.0;
	unsigned QL = queryData.getSeqLength();
	vector<unsigned> MaxIdVec(QL, 0);
	for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex) {
		unsigned ParentSeqIndex = Parents[ParentIndex];

		SeqData parentData = database->getSeqData(ParentSeqIndex);

		PSDs.push_back(parentData);

		// align query to each parent
		PathData PD;
		bool Found = aligner.globalAlign(queryData, parentData, PD);
		if (!Found) {
			Paths.push_back("");
			continue;
		}

		double PctId = 100.0*getFractIdGivenPath(queryData.getSeq(), parentData.getSeq(), PD.Start);
		
		if (PctId > TopPctId) { TopPctId = PctId; }

		string Path = PD.Start;
		Paths.push_back(Path);

		vector<unsigned> IdVec = getSmoothedIdVec(queryData, parentData, Path, opt->getIdsmoothwindow());

		for (unsigned QPos = 0; QPos < QL; ++QPos) {
			if (IdVec[QPos] > MaxIdVec[QPos]) {
				MaxIdVec[QPos] = IdVec[QPos];
			}
		}
	}

	vector<unsigned> BestParents;
	for (unsigned k = 0; k < opt->getMaxp(); ++k)
		{
		unsigned BestParent = UINT_MAX;
		unsigned BestCov = 0;
		for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
			{
			const SeqData &parentData = PSDs[ParentIndex];
			const string &Path = Paths[ParentIndex];
			if (Path == "")
				continue;

			vector<unsigned> IdVec = getSmoothedIdVec(queryData, parentData, Path, 
											opt->getIdsmoothwindow());
			unsigned Cov = 0;
			for (unsigned QPos = 0; QPos < QL; ++QPos)
				if (IdVec[QPos] == MaxIdVec[QPos])
					++Cov;

			if (Cov > BestCov)
				{
				BestParent = ParentIndex;
				BestCov = Cov;
				}
			}

		if (BestParent == UINT_MAX)
			break;

		BestParents.push_back(BestParent);

		const SeqData &parentData = PSDs[BestParent];
		const string &Path = Paths[BestParent];
		vector<unsigned> IdVec = getSmoothedIdVec(queryData, parentData, Path,
								  opt->getIdsmoothwindow());
		for (unsigned QPos = 0; QPos < QL; ++QPos)
			if (IdVec[QPos] == MaxIdVec[QPos])
				MaxIdVec[QPos] = UINT_MAX;
		}

	unsigned BestParentCount = SIZE(BestParents);

	bool chimeric = false;
	for (unsigned k1 = 0; k1 < BestParentCount; ++k1) {
		unsigned i1 = BestParents[k1];

		const SeqData &PSD1 = PSDs[i1];
		const string &Path1 = Paths[i1];

		for (unsigned k2 = k1 + 1; k2 < BestParentCount; ++k2) {
			unsigned i2 = BestParents[k2];

			const SeqData &PSD2 = PSDs[i2];
			const string &Path2 = Paths[i2];

			ChimeHit2 Hit2 = chimeAlign->alignChime(queryData, PSD1, PSD2, Path1, Path2);
			Hit2.PctIdQT = TopPctId;

			if (Hit2.Accept(opt->getMinh(), opt->getMindiv(), opt->getMindiffs())) {
				chimeric = true;
			}
			if (Hit2.Score > Hit.Score) {
				Hit = Hit2;
			}
		}
	}

	return chimeric;
}
/******************************************************************************/
double SearchChime::getFractIdGivenPath(string A,  string B, const char *Path) {
	
	unsigned PosA = 0;
	unsigned PosB = 0;
	unsigned Ids = 0;
	unsigned Wilds = 0;

	unsigned ColCount = (unsigned) strlen(Path);
	if (ColCount == 0) {
		return 0.0;
	}

	for (const char *p = Path; *p; ++p) {

		char c = *p;
		if (c == 'M') {
			Byte a = toupper(A[PosA]);
			Byte b = toupper(B[PosB]);
			if (g_IsACGTU[a] && g_IsACGTU[b]) {
				if (a == b) {
					++Ids;
				}
			}else {
				++Wilds;
			}
		}
		if (c == 'M' || c == 'D')
			++PosA;
		if (c == 'M' || c == 'I')
			++PosB;
	}

	unsigned MinLen = min(PosA, PosB) - Wilds;
	double FractId = (MinLen == 0 ? 0.0 : double(Ids)/double(MinLen));
	
	return FractId;
}
/******************************************************************************/
