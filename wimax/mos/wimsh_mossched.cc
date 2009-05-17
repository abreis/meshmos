#include <wimsh_mossched.h>
#include <wimsh_scheduler_frr.h>
#include <stat.h>
#include <math.h>

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
}

void
WimshMOSScheduler::statSDU(WimaxSdu* sdu)
{

//	if(stats_.size() == 0) { // first run
//		// create an empty MOSFlowInfo and push it
//		MOSFlowInfo flowstat (sdu->flowId(), HDR_CMN(sdu->ip())->uid() - 1);
//		stats_.push_back(flowstat);
//		fprintf(stderr, "%.9f WMOS::statSDU    [%d] Adding flow %d to stats\n",
//				NOW, mac_->nodeId(), sdu->flowId());
//	} else {
	{
		bool fid_present = FALSE;
		for (unsigned i=0; i < stats_.size(); i++) {
			// find this flowID in the vector, if it doesn't exist create it
			if(stats_[i].fid_ == sdu->flowId()) {
				fid_present = TRUE;
				break;
			}
		}
		if(!fid_present) {
			MOSFlowInfo flowstat (sdu->flowId(), -1); // this assumes the MOS scheduler always starts before the flows
			stats_.push_back(flowstat);
			fprintf(stderr, "%.9f WMOS::statSDU    [%d] Adding flow %d to stats\n",
					NOW, mac_->nodeId(), sdu->flowId());
		}
	}

	// now update the flow statistics

	// navigate to the position of MOSFlowInfo(FlowID)
	unsigned int n;
	for (n=0; n < stats_.size(); n++)
		{ if (stats_[n].fid_ == sdu->flowId()) break; }

	stats_[n].count_++;

//	fprintf(stderr, "\t\tstats_[n].lastuid_ %d sdu->seqnumber() %d\n", stats_[n].lastuid_, sdu->seqnumber());
 	if( sdu->nHops() != 0) // ignore the first node, no estimates yet
 	{
 		// if the packet uids received are not sequential, assume missing packets
		if( (stats_[n].lastuid_ + 1) != (int)sdu->seqnumber())
			stats_[n].lostcount_ += sdu->seqnumber() - (stats_[n].lastuid_ + 1);
		// update the last UID received
		stats_[n].lastuid_ = sdu->seqnumber();
 	}

	// update packet loss estimate
	stats_[n].loss_ = (float)stats_[n].lostcount_ / (float)(stats_[n].lostcount_ + stats_[n].count_);

	// debug timestamps
//	fprintf(stderr, "\t DEBUG\n\t\tNOW %f\n\t\tTimestamp %f\n",
//			NOW, sdu->timestamp());

	// delay estimates
 	if( sdu->nHops() != 0) // ignore the first node, no estimates yet
 		stats_[n].delay_ = stats_[n].delay_*0.75 + (NOW - sdu->timestamp())*0.25;

 	// debug delay
// 	fprintf(stderr, "\t DEBUG\n\t\tOld %f\n\t\tNew %f\n",
// 			stats_[n].delay_, NOW - sdu->timestamp() );

}

void
WimshMOSScheduler::dropPDU(WimaxPdu* pdu)
{
	// get the sdu
	WimaxSdu* sdu = pdu->sdu();

	/* we assume the first packet of a new flow is never dropped,
	 * so no checks are made to the stats_ vector
	 */

	// navigate to the position of MOSFlowInfo(FlowID)
	unsigned int n;
	for (n=0; n < stats_.size(); n++)
		{ if (stats_[n].fid_ == sdu->flowId()) break; }

	// increase the lost packet count
	stats_[n].lostcount_++;
	// update packet loss estimate
	stats_[n].loss_ = (float)stats_[n].lostcount_ / (float)(stats_[n].lostcount_ + stats_[n].count_);

}

float
WimshMOSScheduler::audioMOS(double delay, float loss)
{
	/* data from:
	 * Improving Quality of VoIP Streams over WiMax
	 */

	float mos = 0;
	int R = 0;
	bool H = FALSE;
	float Id = 0, Ie = 0;

	// effects of delay
	if( (delay-177.3) >= 0 ) H = TRUE;
	Id = 0.024*delay + 0.11*(delay - 177.3)*H;

	// effects of loss
	// for G.711
	int gamma1 = 0;
	int gamma2 = 30;
	int gamma3 = 15;
	Ie = gamma1 + gamma2*log(1+gamma3*loss);

	// R-factor
	R = 94.2 - Ie - Id;

	// MOS
	mos = 1 + 0.035*R + 0.000007*R*(R-60)*(100-R);

	// debug
	fprintf(stderr, "\taudioMOS delay %f loss %f mos %f", delay, loss, mos);

	return mos;

}

float
WimshMOSScheduler::dataMOS (float loss, float rate)
{
	/* data from:
	 * MOS-Based Multiuser Multiapplication Cross-Layer
	 * Optimization for Mobile Multimedia Communication
	 */

	float mos = 0;

	// a=2.1 & b=0.3 will fit the curve for MOS=1 at rate 10kbps and MOS 4.5 at rate 450kbps
	float data_a = 2.1;
	float data_b = 0.3;

	mos = data_a * log10(data_b*rate*(1-loss));

	// truncate the mos at maximum value
	if(mos>4.5) mos=4.5;

	return mos;
}

float
WimshMOSScheduler::videoMOS (void)
{
	float mos = 0;
	float psnr = 0;
	float mse = 0;

	psnr = 10*log10(255*255/mse);

	// linear mapping from PSNR 20dB (MOS 1) to 40dB (MOS 5)
	mos = psnr*0.20 - 3;

	return mos;
}

void
WimshMOSScheduler::trigger(void)
{
	fprintf(stderr, "%.9f WMOS::trigger    [%d] MOS Scheduler timer fired\n",
			NOW, mac_->nodeId());

	// run the buffer algorithms
	bufferMOS();

	// print some statistics, for debugging
	fprintf (stderr,"\tflow statistics:\n");
	for (unsigned i=0; i < stats_.size(); i++) {
		fprintf (stderr,"\t\tfid %d lastuid %d count %d lost %d lossrate %f delay %f\n",
			stats_[i].fid_, stats_[i].lastuid_, stats_[i].count_, stats_[i].lostcount_, stats_[i].loss_,
			stats_[i].delay_);
	}

}

void
WimshMOSScheduler::bufferMOS(void)
{
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
					if(pdulist_[i][k][l]->sdu()->ip()->datalen()) {
						if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOD_DATA) {
							VideoData* vodinfo_ = (VideoData*)pdulist_[i][k][l]->sdu()->ip()->userdata();
							fprintf (stderr, "\t\tVOD_DATA\tfid %d ndx %d distortion %f\n", pdulist_[i][k][l]->sdu()->flowId(), i, vodinfo_->distortion());
						} else if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOIP_DATA) {
							fprintf (stderr, "\t\tVOIP_DATA\tfid %d ndx %d\n", pdulist_[i][k][l]->sdu()->flowId(), i);
						}
					} else {
							fprintf (stderr, "\t\tFTP_DATA\tfid %d ndx %d\n", pdulist_[i][k][l]->sdu()->flowId(), i);
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
			buffusage, sched_->maxBufSize(), ((float)buffusage/(float)sched_->maxBufSize())*100 );

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
