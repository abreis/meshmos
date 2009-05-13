

#ifndef __NS2_WIMSH_MOS_SCHEDULER_H
#define __NS2_WIMSH_MOS_SCHEDULER_H

#include <wimax_common.h>
#include <t_timers.h>

#include <scheduler.h>
#include <rng.h>

#include <math.h>


class WimshMac;
class WimshMshDsch;

class WimshMOSScheduler : public TclObject {
public:
	//! Do nothing.
	WimshMOSScheduler ();
	//! Do nothing.
	virtual ~WimshMOSScheduler () { }

	//! Handle the propagation timer: dispatch burst.
	void handle ();

private:

protected:
	//! Tcl interface.
	virtual int command(int argc, const char*const* argv);

	//! Pointer to the MAC layer.
	WimshMac* mac_;
};



#endif __NS2_WIMSH_MOS_SCHEDULER_H
