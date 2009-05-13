

#ifndef __NS2_WIMSH_MOS_SCHEDULER_H
#define __NS2_WIMSH_MOS_SCHEDULER_H

#include <wimax_common.h>
#include <timer-handler.h>

#include <scheduler.h>
#include <rng.h>

#include <math.h>


class WimshMac;
class WimshMshDsch;
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

	MOStimer& gettimer() { return timer_; }
	//! Handle the propagation timer: dispatch burst.
	void handle ();

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
