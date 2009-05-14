#include <wimsh_mossched.h>
#include <wimsh_scheduler_frr.h>
#include <stat.h>

static class WimshMOSSchedulerClass : public TclClass {
public:
	WimshMOSSchedulerClass() : TclClass("WimshMOSScheduler") {}
	TclObject* create(int, const char*const*) {
		return (new WimshMOSScheduler);
	}
} class_wimsh_mosscheduler;

WimshMOSScheduler::WimshMOSScheduler () : timer_(this)
{
	fprintf(stderr, "Initialized a MOS Scheduler\n");
	timer_.resched(0.010);
}

int
WimshMOSScheduler::command(int argc, const char*const* argv)
{
	if ( argc == 3 && strcmp (argv[1], "mac") == 0 ) {
		mac_ = (WimshMac*) TclObject::lookup(argv[2]);
		return TCL_OK;
	} else if ( argc == 3 && strcmp (argv[1], "scheduler") == 0 ) {
		sched_ = (WimshSchedulerFairRR*) TclObject::lookup(argv[2]);
		return TCL_OK;
	}
	return TCL_ERROR;

//	} else if ( argc == 3 && strcmp (argv[1], "propagation") == 0 ) {
//		propagation_ = 1.0e-6 * atof (argv[2]);   // in us
//		return TCL_OK;
//	} else if ( argc == 3 && strcmp (argv[1], "id") == 0 ) {
//		uid_ = (unsigned int) atoi (argv[2]);
//		return TCL_OK;
}

void
WimshMOSScheduler::trigger(void)
{
	fprintf (stderr,"%.9f Timer PING on node %d\n", NOW, mac_->nodeId());

	// pointer to this node's scheduler
	WimshSchedulerFairRR* sched_ = (WimshSchedulerFairRR*)mac_->scheduler();


	// this scheduler won't work with per-flow or per-link buffer sharing
	if (sched_->BufferMode() == WimshSchedulerFairRR::PER_FLOW ||
			sched_->BufferMode() == WimshSchedulerFairRR::PER_LINK)
		{
			fprintf (stderr, "Warning: no buffer sharing, MOS scheduler will not function.\n");
			return;
		}

	// TODO: scheduler code enters here
	// Probably we need to do SHARED or PER_LINK buffer instead of PER_FLOW
	// WimaxPdu* pdu
	// desc.size_ -> buffer size
	// pdu->size()
	// link_[ndx][s] -> link queue
	// 	.rr_ -> Round robin list of flow descriptors
	// 	.queue_ -> Packet queue to go to the destination dst_ (std::queue<WimaxPdu*>)

	// a std::queue cannot be manipulated, so we either:
	//	- change all queues to vectors or lists
	//	- dump the queues, apply the algorithms, and then rebuild them

	// RDscheduler(WimaxPdu* pdu)

//	if(pdu->sdu()->ip()->datalen())
//		if(pdu->sdu()->ip()->userdata()->type() == VOD_DATA) {
//			VideoData* vodinfo_ = (VideoData*)pdu->sdu()->ip()->userdata();
//			fprintf (stderr, "\tGot a VOD_DATA packet, distortion %f\n", vodinfo_->distortion());
//		}

	fprintf (stderr, "\tPackets in the buffers:\n");

	// vector array to store all PDUs in the buffers
	// [ndx][serv][queueindex][pdu]
//	std::vector< std::vector< std::vector< std::vector< WimaxPdu* > > > > pdulist_;
	std::vector< std::vector< std::vector< WimaxPdu* > > > pdulist_;

	// resize the array
	pdulist_.resize(mac_->nneighs());
//	for(unsigned i=0; i < mac_->nneighs(); i++)
//		pdulist_[i].resize(wimax::N_SERV_CLASS);


	// pop all PDUs from the queues, to process
	// run all neighbors
	for(unsigned i=0; i < mac_->nneighs(); i++) {
		// run all services
//		for(unsigned j=0; j < wimax::N_SERV_CLASS; j++) {
			// get the list of packet queues
			std::list<WimshSchedulerFairRR::FlowDesc> list_ = sched_->Link()[i].rr_.list();
			list<WimshSchedulerFairRR::FlowDesc>::iterator iter1 = list_.begin();

			// queue index, for indexing
			unsigned int qindex = 0;

			// resize the array for n packet queues
			unsigned lsize = list_.size();
			pdulist_[i].resize(lsize);

			// for each packet queue
			while( iter1 != list_.end() ) {
				// for each PDU
				while(!iter1->queue_.empty()) {
					// pop a PDU from the list into the vector
					WimaxPdu* temppdu = iter1->queue_.front();
					pdulist_[i][qindex].push_back(temppdu);
					iter1->queue_.pop();
				}
				// increment the queue index
				qindex++;
				// get the next packet queue
				iter1++;
			}
//		}
	}


	// show a list of all VOD_DATA packets in this node's buffers
	for(unsigned i=0; i < pdulist_.size(); i++)
//		for(unsigned j=0; j < pdulist_[i].size(); j++)
			for(unsigned k=0; k < pdulist_[i].size(); k++)
				for(unsigned l=0; l < pdulist_[i][k].size(); l++) {
					if(pdulist_[i][k][l]->sdu()->ip()->datalen())
						if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOD_DATA) {
							VideoData* vodinfo_ = (VideoData*)pdulist_[i][k][l]->sdu()->ip()->userdata();
							fprintf (stderr, "\t\tVOD_DATA fid %d ndx %d distortion %f\n", pdulist_[i][k][l]->sdu()->flowId(), i, vodinfo_->distortion());
						}
				}

	// fill a vector with all flow ids of the packets in the buffers
	std::vector <int> flowids_;

	for(unsigned i=0; i < pdulist_.size(); i++)
//		for(unsigned j=0; j < pdulist_[i].size(); j++)
			for(unsigned k=0; k < pdulist_[i].size(); k++)
				for(unsigned l=0; l < pdulist_[i][k].size(); l++) {
					// vector empty, push
					if(flowids_.size() == 0) {
						flowids_.push_back(pdulist_[i][k][l]->sdu()->flowId());
					} else {
						// see if the fid is already in the list
						bool inlist_ = false;
						for(unsigned int m=0; m < flowids_.size(); m++) {
							if(flowids_[m] == pdulist_[i][k][l]->sdu()->flowId()) {
								inlist_ = true;
							}
						}
						// not on the list, push
						if(!inlist_)
							flowids_.push_back(pdulist_[i][k][l]->sdu()->flowId());
					}
				}
//
	// print all flowIDs
	fprintf (stderr, "\tFlow IDs in the buffer: ");
	for(unsigned int i=0; i < flowids_.size(); i++)
		fprintf (stderr, "%d ", flowids_[i]);
	fprintf (stderr, " \n");

	// see how full the buffer is
	unsigned int buffusage = 0;
	for(unsigned i=0; i < pdulist_.size(); i++)
//		for(unsigned j=0; j < pdulist_[i].size(); j++)
			for(unsigned k=0; k < pdulist_[i].size(); k++)
				for(unsigned l=0; l < pdulist_[i][k].size(); l++) {
					buffusage += pdulist_[i][k][l]->size();
				}
	fprintf (stderr, "\tBuffer usage %d/%d, %.2f%%\n",
			buffusage, sched_->bufSize(), ((float)buffusage/(float)sched_->maxBufSize())*100 );

	// reconstruct queues
	for(unsigned i=0; i < mac_->nneighs(); i++) {
		// run all services
//		for(unsigned j=0; j < wimax::N_SERV_CLASS; j++) {
			std::list<WimshSchedulerFairRR::FlowDesc> list_ = sched_->Link()[i].rr_.list();
			list<WimshSchedulerFairRR::FlowDesc>::iterator iter1 = list_.begin();

			// queue index, for indexing
			unsigned int qindex = 0;

			// for each packet queue
			while( iter1 != list_.end() ) {
				// for each PDU
				for(unsigned k=0; k < pdulist_[i][qindex].size(); k++)
					iter1->queue_.push(pdulist_[i][qindex][k]);
				iter1++;
			}

//		}
	}

}




void
MOStimer::expire(Event *e) {
	// call the handle function for the local node
	a_->trigger();

	// reschedule
	a_->gettimer().resched(0.050); // a_->interval_ (and define via TCL)
}
