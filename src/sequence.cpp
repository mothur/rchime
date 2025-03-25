/*
 *  sequence.cpp
 *
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "sequence.h"

/***********************************************************************/
Sequence::Sequence(){
    initialize();
}
/***********************************************************************/
Sequence::Sequence(string newName, string sequence) {
    initialize();
    name = newName;

    //setUnaligned removes any gap characters for us
    setUnaligned(sequence);
    setAligned(sequence);
}
/***********************************************************************/
void Sequence::initialize(){
	name = "";
	unaligned = "";
	aligned = "";
	pairwise = "";
    comment = "";

	numBases = 0;
	alignmentLength = 0;
	startPos = -1;
	endPos = -1;
	longHomoPolymer = -1;
	ambigBases = -1;
}
/***********************************************************************/
void Sequence::setName(string seqName) {
	if(seqName[0] == '>')	{	seqName = seqName.substr(1);	}
	name = seqName;
}
/***********************************************************************/
void Sequence::setUnaligned(string sequence){

	if(sequence.find_first_of('.') != string::npos ||
        sequence.find_first_of('-') != string::npos) {
		string temp = "";
		for(int j=0;j<sequence.length();j++) {
			if(isalpha(sequence[j]))	{	temp += sequence[j];	}
		}
		unaligned = temp;
	} else {
		unaligned = sequence;
	}
	numBases = unaligned.length();
    startPos = -1;
    endPos = -1;
    longHomoPolymer = -1;
    ambigBases = -1;
}
/***********************************************************************/
void Sequence::setAligned(string sequence){

    for (auto i = 0; i < sequence.length(); i++) {
        sequence[i] = ::toupper(sequence[i]);
        if(sequence[i] == 'U'){ sequence[i] = 'T'; }
        if(sequence[i] != '.' && sequence[i] != '-' && sequence[i] != 'A' &&
           sequence[i] != 'T' && sequence[i] != 'G'  && sequence[i] != 'C' &&
           sequence[i] != 'N'){
            sequence[i] = 'N';
        }
    }

	//if the alignment starts or ends with a gap, replace it with a period to indicate missing data
	aligned = sequence;
	alignmentLength = aligned.length();
	setUnaligned(sequence);

	if(aligned[0] == '-'){
		for(int i=0;i<alignmentLength;i++){
			if(aligned[i] == '-'){
				aligned[i] = '.';
			}
			else{
				break;
			}
		}
		for(int i=alignmentLength-1;i>=0;i--){
			if(aligned[i] == '-'){
				aligned[i] = '.';
			}
			else{
				break;
			}
		}
	}
}
/***********************************************************************/
bool Sequence::isAligned(){
    for (int i = 0; i < aligned.length(); i++) {
        if ((aligned[i] == '.') || (aligned[i] == '-')) { return true; }
    }
    return false;
}
/***********************************************************************/
string Sequence::convert2ints() {
	if(unaligned == "")	{	/* need to throw an error */	}

	string processed = unaligned;

	for (auto i = 0; i < processed.length(); i++) {
	    if(processed[i] == 'A')                {    processed[i] = '0';    }
	    else if(processed[i] == 'C')           {    processed[i] = '1';    }
	    else if(processed[i] == 'G')           {    processed[i] = '2';    }
	    else if(processed[i] == 'T')           {    processed[i] = '3';    }
	    else if(processed[i] == 'U')           {    processed[i] = '3';    }
	    else                                   {    processed[i] = '4';    }
	}
    //iterate through string - replace bases with ints

	return processed;
}
/***********************************************************************/
string Sequence::getInlineSeq(){
	return name + '\t' + aligned;
}
/***********************************************************************/
int Sequence::getNumNs(){
    int numNs = 0;
	for (int i = 0; i < unaligned.length(); i++) {
        if(unaligned[i] == 'N') { numNs++; }
    }
    return numNs;
}
/***********************************************************************/
int Sequence::getAmbigBases(){
	if(ambigBases == -1){
		ambigBases = 0;
		for(int j=0;j<numBases;j++){
			if(unaligned[j] != 'A' && unaligned[j] != 'T' &&
                unaligned[j] != 'G' && unaligned[j] != 'C'){
				ambigBases++;
			}
		}
	}

	return ambigBases;
}
/***********************************************************************/
void Sequence::removeAmbigBases(){
	for(int j=getStartPos()-1;j<getEndPos();j++){
		if(aligned[j] != 'A' && aligned[j] != 'T' && aligned[j] != 'G' &&
            aligned[j] != 'C'){
			aligned[j] = '-';
		}
	}
	setUnaligned(aligned);
}
/***********************************************************************/
int Sequence::getLongHomoPolymer(){
    if(longHomoPolymer == -1){
        longHomoPolymer = 1;
        int homoPolymer = 1;
        for(int j=1;j<numBases;j++){
            if(unaligned[j] == unaligned[j-1]){
                homoPolymer++;
            }
            else{
                if(homoPolymer > longHomoPolymer) {
                    longHomoPolymer = homoPolymer;
                }
                homoPolymer = 1;
            }
        }
        if(homoPolymer > longHomoPolymer){	longHomoPolymer = homoPolymer;	}
    }
    return longHomoPolymer;
}
/***********************************************************************/
int Sequence::getStartPos(){
    bool seqAligned = isAligned();

	if(startPos == -1){
		for(int j = 0; j < alignmentLength; j++) {
			if((aligned[j] != '.')&&(aligned[j] != '-')){
				startPos = j + 1;
				break;
            }
		}
	}

	if(!seqAligned){	startPos = 1;	}

	return startPos;
}
/***********************************************************************/
void Sequence::filterToPos(int start){

    if (start > aligned.length()) { start = aligned.length(); }

	for(int j = 0; j < start; j++) { aligned[j] = '.'; }

    //things like ......----------AT become ................AT
    for(int j = start; j < aligned.length(); j++) {
        if (isalpha(aligned[j])) { break; }
        else { aligned[j] = '.'; }
    }
    setUnaligned(aligned);
}
/***********************************************************************/
void Sequence::filterFromPos(int end){

    if (end > aligned.length()) { end = aligned.length(); }

	for(int j = end; j < aligned.length(); j++) {
		aligned[j] = '.';
	}

    for(int j = aligned.length()-1; j >= 0; j--) {
        if (isalpha(aligned[j])) { break; }
        else { aligned[j] = '.'; }
    }

    setUnaligned(aligned);
}
/***********************************************************************/
int Sequence::getEndPos(){
    bool seqAligned = isAligned();

	if(endPos == -1){
		for(int j=alignmentLength-1;j>=0;j--){
			if((aligned[j] != '.')&&(aligned[j] != '-')){
				endPos = j + 1;
				break;
            }
		}
	}
	if(!seqAligned){	endPos = numBases;	}

	return endPos;
}
/***********************************************************************/
void Sequence::padToPos(int start){

    for(int j = getStartPos()-1; j < start-1; j++) {
        aligned[j] = '.';
    }
    startPos = start;

}
/***********************************************************************/
void Sequence::padFromPos(int end){

	for(int j = end; j < getEndPos(); j++) {
		aligned[j] = '.';
	}
	endPos = end;

}
/***********************************************************************/
void Sequence::reverseComplement(){

	string temp;
	for(int i=numBases-1;i>=0;i--){
		if(unaligned[i] == 'A')		{	temp += 'T';	}
		else if(unaligned[i] == 'T'){	temp += 'A';	}
		else if(unaligned[i] == 'G'){	temp += 'C';	}
		else if(unaligned[i] == 'C'){	temp += 'G';	}
		else						{	temp += 'N';	}
	}

	setAligned(temp);
}
/***********************************************************************/
void Sequence::trim(int length){

	if(numBases > length){
		unaligned = unaligned.substr(0,length);
		numBases = length;
        setAligned(unaligned);
	}
}
/***********************************************************************/
