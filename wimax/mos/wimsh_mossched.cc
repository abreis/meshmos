#include <wimsh_mossched.h>
#include <wimsh_topology.h>
#include <wimsh_phy.h>
#include <wimsh_packet.h>

#include <stat.h>

static class WimshMOSSchedulerClass : public TclClass {
public:
	WimshMOSSchedulerClass() : TclClass("WimshMOSScheduler") {}
	TclObject* create(int, const char*const*) {
		return (new WimshMOSScheduler);
	}
} class_wimsh_mosscheduler;

WimshMOSScheduler::WimshMOSScheduler () // : timer_ (this)
{
	fprintf(stderr, "Initialized a MOS Scheduler\n");
}

int
WimshMOSScheduler::command(int argc, const char*const* argv)
{
	if ( argc == 3 && strcmp (argv[1], "mac") == 0 ) {
		mac_ = (WimshMac*) TclObject::lookup(argv[2]);
		return TCL_OK;
//	} else if ( argc == 3 && strcmp (argv[1], "propagation") == 0 ) {
//		propagation_ = 1.0e-6 * atof (argv[2]);   // in us
//		return TCL_OK;
//	} else if ( argc == 3 && strcmp (argv[1], "id") == 0 ) {
//		uid_ = (unsigned int) atoi (argv[2]);
//		return TCL_OK;
	}
	return TCL_ERROR;
}

void
WimshMOSScheduler::handle ()
{

}
