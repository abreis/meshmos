

#ifndef __NS2_WIMSH_MOS_SCHEDULER_H
#define __NS2_WIMSH_MOS_SCHEDULER_H

#include <wimax_common.h>
#include <timer-handler.h>

#include <wimsh_mac.h>

#include <videodata.h>
#include <ns-process.h>

#include <vector>

class WimshMac;
class WimshMOSScheduler;
class WimshSchedulerFairRR;

enum MOStraffic { M_VOD, M_VOIP, M_FTP, M_NTRAFFIC };

struct MOSFlowInfo {
	//! Flow ID
	int fid_;
	//! Flow type
	MOStraffic traffic_;
	//! Last packet UID received
	int lastuid_;
	//! Packet count
	unsigned int count_;
	//! Lost packet count
	unsigned int lostcount_;
	//! Packet loss estimate
	float loss_;

	//! Delay information
	double delay_;

	//! MOS of the flow
	float mos_;

	//! distortion information
	// list of lost video packet MSEs
	vector<float> mse_;
	// list of lost video packet IDs
	vector<int> vod_id_;


	//! Constructor
	MOSFlowInfo (int fid = 0, int uid = 0) {
		fid_ = fid; count_ = 0;
		lostcount_ = 0; lastuid_ = uid; loss_ = 0;
		delay_ = 0; mos_ = 0;
	}
};

class MOStimer : public TimerHandler {
public:
	MOStimer(WimshMOSScheduler *a) : TimerHandler() { a_ = a; }
protected:
	virtual void expire(Event *e);
	WimshMOSScheduler *a_;
};

class WimshMOSScheduler : public TclObject {
public:
	//! Do nothing.
	WimshMOSScheduler ();
	//! Do nothing.
	virtual ~WimshMOSScheduler () { }

	//! Called each time the timer fires
	void trigger(void);
	//! Return the timer object
	MOStimer& gettimer() { return timer_; }
	//! Handle the timer event
	void handle ();

	//! Apply the scheduler algoritm to the buffers
	void bufferMOS(void);

	//! Process an SDU for statistics
	void statSDU(WimaxSdu* sdu);
	//! Process a dropped PDU for statistics
	void dropPDU(WimaxPdu* pdu);

	//! Obtain MOS for an audio flow
	float audioMOS (double delay, float loss);
	//! Obtain MOS for a video flow
	float videoMOS (vector<float>* mse, float loss);
	//! Obtain deltaMOS for a video flow
	float deltaVideoMOS (vector<float>* mse, float distincrease, float loss, unsigned packetdrops);
	//! Obtain MOS for a data flow
	float dataMOS (float loss, float rate);

	//! Vector of MOSFlowInfo structs to keep track of data
	std::vector <MOSFlowInfo> stats_;
private:

protected:
	//! Tcl interface.
	virtual int command(int argc, const char*const* argv);

	//! Parameters for data MOS calculation (set via TCL)
//	float data_a=0, data_b=0;

	//! Trigger timer
	MOStimer timer_;

	//! Pointer to the MAC layer.
	WimshMac* mac_;
	//! Pointer to the Scheduler
	WimshSchedulerFairRR* sched_;

	//! convert a long int to a binary value, stored in a vector<bool>
	void dec2bin(long decimal, vector<bool>* binary);
};



#endif __NS2_WIMSH_MOS_SCHEDULER_H
