/*********************************************************************************
*  Copyright (c) 2010-2011, Elliott Cooper-Balis
*                             Paul Rosenfeld
*                             Bruce Jacob
*                             University of Maryland 
*                             dramninjas [at] gmail [dot] com
*  All rights reserved.
*  
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*  
*     * Redistributions of source code must retain the above copyright notice,
*        this list of conditions and the following disclaimer.
*  
*     * Redistributions in binary form must reproduce the above copyright notice,
*        this list of conditions and the following disclaimer in the documentation
*        and/or other materials provided with the distribution.
*  
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
*  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
*  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*********************************************************************************/




//MemorySystem.cpp
//
//Class file for JEDEC memory system wrapper
//

#include "MemorySystem.h"
#include "IniReader.h"
#include <unistd.h>
#include <math.h> 

using namespace std;


ofstream cmd_verify_out; //used in Rank.cpp and MemoryController.cpp if VERIFICATION_OUTPUT is set

unsigned NUM_DEVICES;
unsigned NUM_RANKS;
unsigned NUM_RANKS_LOG;

namespace DRAMSim {

powerCallBack_t MemorySystem::ReportPower = NULL;

MemorySystem::MemorySystem(const string &refreshMapFilename_,unsigned id, unsigned int megsOfMemory, CSVWriter &csvOut_, ostream &dramsim_log_) :

		dramsim_log(dramsim_log_),
		ReturnReadData(NULL),
		WriteDataDone(NULL),
		systemID(id),
		csvOut(csvOut_),
		refreshMapFilename(refreshMapFilename_)
{
	currentClockCycle = 0;
	DEBUG("===== MemorySystem "<<systemID<<" =====");
	//	NUM_RANKS = 1;

	//calculate the total storage based on the devices the user selected and the number of

	//calculate number of devices
	/************************
	  This code has always been problematic even though it's pretty simple. I'll try to explain it 
	  for my own sanity. 

	  There are two main variables here that we could let the user choose:
	  NUM_RANKS or TOTAL_STORAGE.  Since the density and width of the part is
	  fixed by the device ini file, the only variable that is really
	  controllable is the number of ranks. Users care more about choosing the
	  total amount of storage, but with a fixed device they might choose a total
	  storage that isn't possible. In that sense it's not as good to allow them
	  to choose TOTAL_STORAGE (because any NUM_RANKS value >1 will be valid).

	  However, users don't care (or know) about ranks, they care about total
	  storage, so maybe it's better to let them choose and just throw an error
	  if they choose something invalid. 

	  A bit of background: 

	  Each column contains DEVICE_WIDTH bits. A row contains NUM_COLS columns.
	  Each bank contains NUM_ROWS rows. Therefore, the total storage per DRAM device is: 
	  		PER_DEVICE_STORAGE = NUM_ROWS*NUM_COLS*DEVICE_WIDTH*NUM_BANKS (in bits)

	 A rank *must* have a 64 bit output bus (JEDEC standard), so each rank must have:
	  		NUM_DEVICES_PER_RANK = 64/DEVICE_WIDTH  
			(note: if you have multiple channels ganged together, the bus width is 
			effectively NUM_CHANS * 64/DEVICE_WIDTH)
	 
	If we multiply these two numbers to get the storage per rank (in bits), we get:
			PER_RANK_STORAGE = PER_DEVICE_STORAGE*NUM_DEVICES_PER_RANK = NUM_ROWS*NUM_COLS*NUM_BANKS*64 

	Finally, to get TOTAL_STORAGE, we need to multiply by NUM_RANKS
			TOTAL_STORAGE = PER_RANK_STORAGE*NUM_RANKS (total storage in bits)

	So one could compute this in reverse -- compute NUM_DEVICES,
	PER_DEVICE_STORAGE, and PER_RANK_STORAGE first since all these parameters
	are set by the device ini. Then, TOTAL_STORAGE/PER_RANK_STORAGE = NUM_RANKS 

	The only way this could run into problems is if TOTAL_STORAGE < PER_RANK_STORAGE,
	which could happen for very dense parts.
	*********************/
	std::array<std::vector<std::tuple<unsigned , unsigned>> * , 262144> t1;
	std::vector<std::tuple<unsigned , unsigned>> * t2 ;
	// number of bytes per rank
	unsigned long megsOfStoragePerRank = ((((long long)NUM_ROWS * (NUM_COLS * DEVICE_WIDTH) * NUM_BANKS) * ((long long)JEDEC_DATA_BUS_BITS / DEVICE_WIDTH)) / 8) >> 20;

	// If this is set, effectively override the number of ranks
	if (megsOfMemory != 0)
	{
		NUM_RANKS = megsOfMemory / megsOfStoragePerRank;
		NUM_RANKS_LOG = dramsim_log2(NUM_RANKS);
		if (NUM_RANKS == 0)
		{
			PRINT("WARNING: Cannot create memory system with "<<megsOfMemory<<"MB, defaulting to minimum size of "<<megsOfStoragePerRank<<"MB");
			NUM_RANKS=1;
		}
	}

//////////////////////*******************//////////////////////////////////
	rowsRefreshMap = new uint64_t*[NUM_RANKS];
	for (unsigned i=0; i < NUM_RANKS; i++)
		rowsRefreshMap[i] = new uint64_t[8192];
//////////////////////////////////////////////////////////////////////////
	
	for (unsigned i =0 ; i < NUM_RANKS ; i++){
		//std::array<std::vector<std::tuple<unsigned , unsigned>> * , 8192> t3;
		for (unsigned j=0; j < 262144 ; j++)
			t1.at(j) = new std::vector<std::tuple<unsigned,unsigned>>();
		VRTTable.insert(std::pair<unsigned , std::array<std::vector<std::tuple<unsigned , unsigned>> * , 262144>>(i ,t1 ));
	}

	//////////////////////////////////////////////////////////////////////
	RowsStates = new unsigned* [NUM_RANKS];
	for (unsigned i=0; i < NUM_RANKS;i++){
		RowsStates[i] = new unsigned[262144];
	}

	for (unsigned i =0 ; i < NUM_RANKS ; i++)
		for (unsigned j=0; j < 262144 ; j++){
			RowsStates[i][j] = HIGH_STATE;	}
	//////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////
	std::srand(time(NULL));
	ifstream deviceRefreshMap;
	size_t lineNumber = 0;
	string line;
	//float ranNum;
	
	float sigma = 1;
	float mu = 7;
	float refTime ;
	float temp = 10000000.0;
	std::string::size_type sz;
	double time_;
	unsigned state_;
	unsigned cellNumber = 0;
	unsigned cell_id = 0;
	uint64_t clock_;
	uint64_t tick = (REFRESH_PERIOD/tCK);
	//unsigned** cellInEachRow;

	


	deviceRefreshMap.open(refreshMapFilename.c_str());

	if (deviceRefreshMap.is_open()){

		while (!deviceRefreshMap.eof()){

			lineNumber++;
			getline(deviceRefreshMap , line);

			if (deviceRefreshMap.bad()){
				ERROR("Cannot read device refresh map '" << deviceRefreshMap << "'");
				exit(-1);
			}
			
			if (line.size() == 0)
				continue;

			if (RefreshMap){
				if (line == "##"){
					RefreshMap = false;
					continue;
				}
				refTime = sigma * (std::stof(line)) + mu;
				refTime = refTime * pow(0.535 , (temperature-45)/10);

				if (refTime < temp)
					temp = refTime;

				if (lineNumber % 64 == 0){
					for (unsigned i =0 ; i < NUM_RANKS ; i++){
						for (unsigned j =0 ; j < 32 ; j++){
						rowsRefreshMap[i][((lineNumber/64)-1) + (256*j)] = (uint64_t)((temp / 0.064)); 
						//printf("%lu , %f\n",(uint64_t)((temp / 0.064)) , temp);
						}
					}
					if ((unsigned)((temp / 0.064) >(2*maxRefTest) ))
						maxRefTest *= 2;
					temp = 10000000.0;
				}

			}
			else{
				if (line == "##"){
					cell_id++;
					unsigned modu = 262144;
					cellNumber = std::rand() % modu;
					//modu = NUM_BANKS * NUM_ROWS / 8192;
					//cellNumber = cellNumber / modu;
					//printf("cellNumber:%u cell_id:%u\n", cellNumber , cell_id );
					for(unsigned i = 0; i < NUM_RANKS; i++){
						t1 = VRTTable.at(i);
						t2 = t1.at(cellNumber);
						t2->push_back(std::make_tuple(cell_id , HIGH_STATE));
					}
					continue;
					}
				time_ = std::stod(line , &sz);
				state_= std::stof(line.substr(sz));
				//printf("TIME: %20.18f STATE: %u\n" ,time_ ,  state_);
				clock_ = (time_ * 1e9 ) / tCK ;
				//printf("TIME: %20.18f CLOCK:%lu\n", time_ , clock_ );
				//printf("cellNumber:%u cell_id:%u\n", cellNumber , cell_id );
				VRTCells.insert(std::pair<uint64_t , std::tuple<unsigned,unsigned,unsigned>>( clock_ , std::make_tuple(cellNumber , state_ ,cell_id))) ;
			}
				// read VRT Cells
			}
			
		}
	

	else {
		ERROR("Cannot load Device refresh map '" << refreshMapFilename << "'");
		abort();
	}

//	for (MapIterator ite = VRTCells.begin() ; ite != VRTCells.end() ; ite++)
//		printf("VRTCell: %lu\n", ite->first );

	MapIterator iter = VRTCells.begin();
	nextVRTTableUpdate = iter->first;
	//printf("Next Update at: %lu\n", nextVRTTableUpdate );


//////////////////////////***********************/////////////////////////////////


	NUM_DEVICES = JEDEC_DATA_BUS_BITS/DEVICE_WIDTH;
	TOTAL_STORAGE = (NUM_RANKS * megsOfStoragePerRank); 

	DEBUG("CH. " <<systemID<<" TOTAL_STORAGE : "<< TOTAL_STORAGE << "MB | "<<NUM_RANKS<<" Ranks | "<< NUM_DEVICES <<" Devices per rank");


	memoryController = new MemoryController(this, csvOut, dramsim_log);

	// TODO: change to other vector constructor?
	ranks = new vector<Rank *>();

	for (size_t i=0; i<NUM_RANKS; i++)
	{
		Rank *r = new Rank(dramsim_log);
		r->setId(i);
		r->attachMemoryController(memoryController);
		ranks->push_back(r);
	}

	memoryController->attachRanks(ranks);
/*
	for (unsigned i= 0; i < NUM_RANKS; i++ ){
		t1 = VRTTable.at(i);
		for (unsigned j=0; j < 262144 ; j++){
			t2 = t1.at(j);
			for (std::vector<std::tuple<unsigned,unsigned>>::iterator iter = t2->begin(); iter != t2->end() ; iter++){
				printf("VRTTable: ROW:%u CELLID:%u STATE:%u\n", j, std::get<0>(*iter) , std::get<1>(*iter));
			}
		}
	}
*/	
///////////////////////////////////Initial value of nextTestCycle and testSteps////////////////////////////
	// 8192 * tick equals to around 64ms; depending on what step we want to start n could be changed
	testSteps = 2 * 8192 * tick ;
	nextTestCycle = testSteps;
	//printf("TestSteps: %lu  nextTestCycle:%lu \n", testSteps , nextTestCycle);
}



MemorySystem::~MemorySystem()
{
	/* the MemorySystem should exist for all time, nothing should be destroying it */  
//	ERROR("MEMORY SYSTEM DESTRUCTOR with ID "<<systemID);
//	abort();

	delete(memoryController);

	for (size_t i=0; i<NUM_RANKS; i++)
	{
		delete (*ranks)[i];
	}
	ranks->clear();
	delete(ranks);

	if (VERIFICATION_OUTPUT)
	{
		cmd_verify_out.flush();
		cmd_verify_out.close();
	}
}

bool MemorySystem::WillAcceptTransaction()
{
	return memoryController->WillAcceptTransaction();
}

bool MemorySystem::addTransaction(bool isWrite, uint64_t addr)
{
	TransactionType type = isWrite ? DATA_WRITE : DATA_READ;
	Transaction *trans = new Transaction(type,addr,NULL);
	// push_back in memoryController will make a copy of this during
	// addTransaction so it's kosher for the reference to be local 

	if (memoryController->WillAcceptTransaction()) 
	{
		return memoryController->addTransaction(trans);
	}
	else
	{
		pendingTransactions.push_back(trans);
		return true;
	}
}

bool MemorySystem::addTransaction(Transaction *trans)
{
	return memoryController->addTransaction(trans);
}

//prints statistics
void MemorySystem::printStats(bool finalStats)
{
	memoryController->printStats(finalStats);

	
}


//update the memory systems state
void MemorySystem::update()
{

	//PRINT(" ----------------- Memory System Update ------------------");

	//updates the state of each of the objects
	// NOTE - do not change order
	uint64_t tick = (REFRESH_PERIOD/tCK);
	uint64_t tick_test = 8192 * tick * 16; // around 1 sec
	uint64_t scrubing_tick = tick_test * 40;
	std::array<std::vector<std::tuple<unsigned , unsigned>> * , 262144> t1;
	std::vector<std::tuple<unsigned , unsigned>> * t2 ;
	bool SECDEC = true;
	bool DECTEC = false;
	uint64_t testEnergy_tick = tick_test * 3600;

	for (size_t i=0;i<NUM_RANKS;i++)
	{
		(*ranks)[i]->update();
	}

	//pendingTransactions will only have stuff in it if MARSS is adding stuff
	if (pendingTransactions.size() > 0 && memoryController->WillAcceptTransaction())
	{
		memoryController->addTransaction(pendingTransactions.front());
		pendingTransactions.pop_front();
	}
	memoryController->update();

//////////////////////////TEST SESSION=> shuold be changed to consider testing and updating new refresh values///////////////////////////////////////////////////////////
	if (this->getClock() == nextTestCycle){

		printf("TEST SESSION=%lu\n", this->getClock());
		nextTestCycle += testSteps;

		unsigned clockInSixty = testSteps / 8192;
		unsigned numInSixty = clockInSixty / tick;
		unsigned BUFF_SIZE = 4;
		printf("numInSixty: %u\n", numInSixty );
		for (unsigned k =0; k < BUFF_SIZE ; k++){
			for (unsigned i=0; i < NUM_RANKS; i++){
				//uint64_t temp = memoryController->getRefreshTable(i , nextRowupdate);
				unsigned MaxRefTime = 10000000;
				unsigned MaxRefTime_= 0;
				unsigned offset = 262144/8192;
				for (unsigned m =0; m < offset ; m++){
					memoryController->testEnergy[i]++;
					MaxRefTime_ =  getMaxRefTime(rowsRefreshMap[i][nextRowupdate] , RowsStates[i][(nextRowupdate*offset)+m]);
					if (MaxRefTime_ < MaxRefTime){
						MaxRefTime = MaxRefTime_;
					}
				}
				if (numInSixty < MaxRefTime){
					memoryController->updateRefreshTable(i , nextRowupdate,numInSixty);
					printf("Row:%u Updated at:%lu newRefresh:%u \n", nextRowupdate, this->getClock() , numInSixty );
				}
				else {
					//memoryController->updateRefreshTable(i , nextRowupdate, rowsRefreshMap[i][nextRowupdate]/4);
					memoryController->updateRefreshTable(i , nextRowupdate, 1);
					printf("Row:%u Could'nt update at:%lu newRefresh:%d\n",nextRowupdate, this->getClock() , 1 );
				}
				
			}

			nextRowupdate++;
			if (nextRowupdate == 8192){
				nextRowupdate = 0;
				nextTestCycle -= testSteps;
				unsigned clockInSixty_ = testSteps / 8192;
				unsigned numInSixty_ = clockInSixty_ / tick;
				if (numInSixty_ == (maxRefTest)/2)
					testSteps = 2 * 8192 * tick;
				else 
					testSteps = testSteps * 2;
				nextTestCycle += testSteps;
			} 
		}
	}

//////////////////VRT Cells Session*///////////////////////////////////
	

	if (this->getClock() >= nextVRTTableUpdate){
		//double sec =  this->getClock()* tCK * 1e-9;
		//printf("SECOND:%lu\n" ,this->getClock() );
		for (unsigned i= 0; i < NUM_RANKS; i++ ){
			t1 = VRTTable.at(i);
			unsigned row_n = std::get<0>(VRTCells.at(nextVRTTableUpdate));
			unsigned state_n = std::get<1>(VRTCells.at(nextVRTTableUpdate));
			unsigned cellid_n = std::get<2>(VRTCells.at(nextVRTTableUpdate));
			//printf("NEW: r:%u s:%u c:%u\n", row_n , state_n , cellid_n );
			t2 = t1.at(row_n);
			for (std::vector<std::tuple<unsigned,unsigned>>::iterator iter = t2->begin(); iter != t2->end() ; iter++)
				if (std::get<0>(*iter) == cellid_n){
				//	printf("Before Update: CELLID:%u STATE:%u\n", std::get<0>(*iter) , std::get<1>(*iter));
					std::get<1>(*iter) = state_n;
					break;
				}
				/*
			for (std::vector<std::tuple<unsigned,unsigned>>::iterator iter = t2->begin(); iter != t2->end() ; iter++)
				if (std::get<0>(*iter) == cellid_n){
					printf("After Update: CELLID:%u STATE:%u\n", std::get<0>(*iter) , std::get<1>(*iter));
					break;
				}
				*/
			//printf("Updated at: %lu\n", nextVRTTableUpdate );
			
			unsigned newState = HIGH_STATE;
			for (std::vector<std::tuple<unsigned,unsigned>>::iterator iter = t2->begin(); iter != t2->end() ; iter++){
				if (std::get<1>(*iter) < HIGH_STATE)
					newState = std::get<1>(*iter);
			}
			RowsStates[i][row_n] = newState;
		}


		VRTCells.erase(nextVRTTableUpdate);

		MapIterator iter = VRTCells.begin();
		nextVRTTableUpdate = iter->first;
		//printf("Next Update at: %lu\n", nextVRTTableUpdate );

	}
	

	//////////////////////////////////scrubing//////////////////////////////////////////////////
	if (this->getClock() % scrubing_tick == 0){

		printf("SCRUBING SESSION:%lu\n" ,this->getClock() );
		unsigned SCRUBS_FOUND=0;
		unsigned offset = 262144/8192;
		for (unsigned i= 0; i < NUM_RANKS; i++ ){
			memoryController->scrubingEnergy[i]++;
			t1 = VRTTable.at(i);
			//printf("NEW: r:%u s:%u c:%u\n", row_n , state_n , cellid_n );
			for (unsigned r = 0 ; r < 8192 ; r++){
				uint64_t refNow = memoryController->getRefreshTable(i , r);
				for (unsigned m = 0; m < offset ; m++){
					SCRUBS_FOUND = 0;
					t2 = t1.at(r * offset + m);
					// use RowsStates to decice if a crrectable or uncorectabele error is happened
					for (std::vector<std::tuple<unsigned,unsigned>>::const_iterator iter = t2->begin(); iter != t2->end() ; iter++){
						unsigned MaxRefTime = getMaxRefTime(rowsRefreshMap[i][r] , std::get<1>(*iter));
						if (refNow > MaxRefTime){
							SCRUBS_FOUND++;
						}
					}
					// number of correctable and uncorrectable errors 
					if (SECDEC){	
						if (SCRUBS_FOUND == 1){
							memoryController->updateRefreshTable(i , r , 1);
						}
							
					
						else if (SCRUBS_FOUND > 1){
							printf("scrubing couldn't correct errors at row:%u\n" , r);
						}
					}
					else if (DECTEC){
						if (SCRUBS_FOUND == 2){
							memoryController->updateRefreshTable(i , r , 1);
						}
						else if (SCRUBS_FOUND > 2){
							printf("scrubing couldn't correct errors at row:%u\n" , r);
						}
						
					}
				}
			}
		}
	}
	
/////////////////////////////////Error correction section////////////////////////////////////////////////////
	
if (this->getClock() % tick_test == 0){
		printf("ERROR SESSION:%lu\n" ,this->getClock() );
		unsigned ERRORS_FOUND=0;
		unsigned offset = 262144/8192;
		for (unsigned i= 0; i < NUM_RANKS; i++ ){
			t1 = VRTTable.at(i);
			//printf("NEW: r:%u s:%u c:%u\n", row_n , state_n , cellid_n );
			for (unsigned r = 0 ; r < 8192 ; r++){
				uint64_t refNow = memoryController->getRefreshTable(i , r);
				for (unsigned m = 0; m < offset ; m++){
					ERRORS_FOUND = 0;
					t2 = t1.at(r * offset + m);
					// use RowsStates to decice if a crrectable or uncorectabele error is happened
					for (std::vector<std::tuple<unsigned,unsigned>>::const_iterator iter = t2->begin(); iter != t2->end() ; iter++){
						unsigned MaxRefTime = getMaxRefTime(rowsRefreshMap[i][r] , std::get<1>(*iter));
						if (refNow > MaxRefTime){
							ERRORS_FOUND++;
						}
					}
					// number of correctable and uncorrectable errors 
					if (SECDEC){	
						if (ERRORS_FOUND == 1){
							memoryController->correctableErrors[i][r * offset + m]++;
							//printf("An correctable error happened at row:%u\n" , r);
						}
						else if (ERRORS_FOUND > 1){
							memoryController->uncorrectableErrors[i][r * offset + m]++;
							printf("An uncorrectable error happened at row:%u\n" , r);
						}
					}
					else if (DECTEC){
						if (ERRORS_FOUND == 2){
							memoryController->correctableErrors[i][r * offset + m]++;
							//printf("An correctable error happened at row:%u\n" , r);
						}
						else if (ERRORS_FOUND > 2){
							memoryController->uncorrectableErrors[i][r * offset + m]++;
							printf("An uncorrectable error happened at row:%u\n" , r);
							for (std::vector<std::tuple<unsigned,unsigned>>::const_iterator iter = t2->begin(); iter != t2->end() ; iter++)
								printf("Cells: CELLID:%u STATE:%u\n", std::get<0>(*iter) , std::get<1>(*iter));
						
						}
						
					}
				}
			}
		}

	}
	
	/////////////////////////////////////////////////////////////////////////////////
	if (this->getClock() % testEnergy_tick == 0){
		for (unsigned i=0; i < NUM_RANKS; i++)
				printf("Rank: %u scrubingEnergy: %lu testEnergy: %lu\n" , i , memoryController->scrubingEnergy[i] , memoryController->testEnergy[i] );

		for (unsigned i=0; i < NUM_RANKS; i++)
			for (unsigned j=0 ; j < 262144; j++)
				if (memoryController->correctableErrors[i][j] > 0 || memoryController->uncorrectableErrors[i][j] > 0)
					printf("Rank: %u Row: %u correctableErrors: %u uncorrectableErrors: %u\n" , i ,j, memoryController->correctableErrors[i][j], memoryController->uncorrectableErrors[i][j] );
	}
	
	//simply increments the currentClockCycle field for each object
	for (size_t i=0;i<NUM_RANKS;i++)
	{
		(*ranks)[i]->step();
	}
	memoryController->step();
	this->step();

	//PRINT("\n"); // two new lines
}

void MemorySystem::RegisterCallbacks( Callback_t* readCB, Callback_t* writeCB,
                                      void (*reportPower)(double bgpower, double burstpower,
                                                          double refreshpower, double actprepower))
{
	ReturnReadData = readCB;
	WriteDataDone = writeCB;
	ReportPower = reportPower;
}

} /*namespace DRAMSim */



unsigned MemorySystem::getMaxRefTime( uint64_t rowsRefreshMap , unsigned RowsStates){
	if (RowsStates == HIGH_STATE)
		return (unsigned)rowsRefreshMap;
	else if (RowsStates == MED_STATE){
		unsigned t = (unsigned)((rowsRefreshMap/4)*3);
		return t;
	}

	else if (RowsStates == LOW_STATE){
		unsigned t = (unsigned)((rowsRefreshMap/4)*2);
		return t;
	}
	else if (RowsStates == CRIT_STATE){
		unsigned t = (unsigned)((rowsRefreshMap/4));
		return t;
	}

	else {
		ERROR("WRONG STATE IN RowsStates!");
		abort();
	}

}

// This function can be used by autoconf AC_CHECK_LIB since
// apparently it can't detect C++ functions.
// Basically just an entry in the symbol table
extern "C"
{
	void libdramsim_is_present(void)
	{
		;
	}
}

