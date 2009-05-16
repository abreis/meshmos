

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

struct MOSFlowInfo {
	//! Flow ID
	int fid_;
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
//		//! Timestamp (for statistics collection).
//		double timestamp_;

	//! TODO: distortion information


	//! Constructor
	MOSFlowInfo (int fid = 0, int uid = 0) {
		fid_ = fid; count_ = 0;
		lostcount_ = 0; lastuid_ = uid; loss_ = 0;
		delay_ = 0;
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

	//! Vector of MOSFlowInfo structs to keep track of data
	std::vector <MOSFlowInfo> stats_;
private:

protected:
	//! Tcl interface.
	virtual int command(int argc, const char*const* argv);

	//! Trigger timer
	MOStimer timer_;

	//! Pointer to the MAC layer.
	WimshMac* mac_;
	//! Pointer to the Scheduler
	WimshSchedulerFairRR* sched_;
};



#endif __NS2_WIMSH_MOS_SCHEDULER_H
