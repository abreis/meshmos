

#ifndef __NS2_WIMSH_MOS_SCHEDULER_H
#define __NS2_WIMSH_MOS_SCHEDULER_H

#include <wimax_common.h>
#include <timer-handler.h>

#include <wimsh_mac.h>

#include <videodata.h>
#include <ns-process.h>

class WimshMac;
class WimshMOSScheduler;
class WimshSchedulerFairRR;

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

	//! Process an SDU for statistics
	void statSDU(const WimaxSdu* sdu);

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
