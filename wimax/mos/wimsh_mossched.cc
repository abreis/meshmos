#include <wimsh_mossched.h>
#include <wimsh_scheduler_frr.h>
#include <stat.h>
#include <string>

#include <float.h>
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

//	fprintf(stderr, "\tDEBUG got sdu fid %d id %d\n", sdu->flowId(), sdu->seqnumber() );

	// the following code fills the stats_ vector with MOSFlowInfo elements, as new flowIDs are found
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

			if(sdu->ip()->datalen()) {
				if(sdu->ip()->userdata()->type() == VOD_DATA) {
					flowstat.traffic_ = M_VOD;
				} else if(sdu->ip()->userdata()->type() == VOIP_DATA) {
					flowstat.traffic_ = M_VOIP;
				}
			} else {
				flowstat.traffic_ = M_FTP;
			}

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

 	// distortion tracking
 	// being done on drops
//	if(sdu->ip()->userdata()->type() == VOD_DATA) {
//		VideoData* vodinfo_ = (VideoData*)sdu->ip()->userdata();
//		// store distortion info
//		stats_[n].mse_.push_back(vodinfo_->distortion());
//		// store packet ID
//		stats_[n].vod_id_.push_back(sdu->seqnumber());
//
////		fprintf(stderr, "\tDEBUG storing fid %d id %d\n", sdu->flowId(), sdu->seqnumber());
//
//		// keep vectors at length 30
//		if(stats_[n].mse_.size() > 30)
//			stats_[n].mse_.erase(stats_[n].mse_.begin()); // pop front
//		if(stats_[n].vod_id_.size() > 30)
//			stats_[n].vod_id_.erase(stats_[n].vod_id_.begin()); // pop front
//
//		fprintf(stderr, "\tDEBUG vod_id first %d last %d\n", stats_[n].vod_id_[0],
//				stats_[n].vod_id_[stats_[n].vod_id_.size()-1] );
//	}


	// re-evaluate the flow's MOS
 	switch(stats_[n].traffic_) {
 	case M_VOIP:
 		stats_[n].mos_ = audioMOS(stats_[n].delay_, stats_[n].loss_);
 		break;

 	case M_VOD:
 		stats_[n].mos_ = videoMOS( &(stats_[n].mse_), stats_[n].loss_ );
 		break;

 	default:
 		break;
 	}

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

 	// distortion tracking
	if(sdu->ip()->userdata()->type() == VOD_DATA) {
		VideoData* vodinfo_ = (VideoData*)sdu->ip()->userdata();
		// store MSE
		stats_[n].mse_.push_back(vodinfo_->distortion());
		// store ID
		stats_[n].vod_id_.push_back(sdu->seqnumber());

		// inform global statistics of the lost packet
        Stat::put ("rd_vod_lost_mse", sdu->flowId(), vodinfo_->distortion());
        Stat::put ("rd_vod_lost_frames", sdu->flowId(), 1);
//		fprintf(stderr, "\tDEBUG STAT fid %d dist %f\n", sdu->flowId(), vodinfo_->distortion());

		// keep vectors at length 30
//		if(stats_[n].mse_.size() > 30)
//			stats_[n].mse_.erase(stats_[n].mse_.begin()); // pop front

//		fprintf(stderr, "\tDEBUG vod_id first %d last %d\n", stats_[n].vod_id_[0],
//				stats_[n].vod_id_[stats_[n].vod_id_.size()-1] );
	}


	// re-evaluate the flow's MOS
 	switch(stats_[n].traffic_) {
 	case M_VOIP:
 		stats_[n].mos_ = audioMOS(stats_[n].delay_, stats_[n].loss_);
 		break;

 	case M_VOD:
 		stats_[n].mos_ = videoMOS( &(stats_[n].mse_), stats_[n].loss_ );
 		break;

 	default:
 		break;
 	}

	fprintf(stderr, "%.9f WMOS::dropPDU    [%d] Received PDU fid %d seq %d \n",
			NOW, mac_->nodeId(), sdu->flowId(), sdu->seqnumber());

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
	delay *= 1000; // to ms
	delay += 25; // 25ms -> encoding time for G.711
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
//	fprintf(stderr, "\t[%d] audioMOS delay %f Id %f loss %f Ie %f mos %f\n",
//			mac_->nodeId(), delay, Id, loss, Ie, mos);

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
WimshMOSScheduler::videoMOS (vector<float>* mse, float loss)
{
	/* data from:
	 * Real-Time Monitoring of Video Quality in IP Networks
	 */

	// if no frames were lost, return MOS 4.5
	if(mse->size() == 0) {
//		fprintf(stderr, "\t[%d] videoMOS no lost frames, returning MOS 4.5\n", mac_->nodeId());
		return 4.5;
	}

	// define some parameters
	float gamma = 0; // attenuation factor
	float T = 10; // P-frames per GOP (10 in a 30-GOP)
	float alpha = ( pow(gamma, T+1) - (T+1)*gamma + T) / ( T*(1-gamma)*(1-gamma) ); // (5)


	// get sigma for (5), which is the average MSE of all dropped frames
	float mse_total = 0;
	for(unsigned i=0; i < mse->size(); i++)
	{
		mse_total += (*mse)[i];
	}

	float sigma = mse_total / mse->size();

	// obtain D1 (5)
	float D1 = alpha*sigma;

	// now obtain distortion D; we assume H.264 and D=s*n\*Pe*L*D1
//	float s = 8; 	// 8 slices per frame, CIF
	float s = 1; 	// since we have the MSE per frame, we consider 1 frame=1 slice
	float L = 1; 	// 1 frame per packet, encoding
	float n = 1;	// bernoulli, each loss affects 1 packet
	float Pe = loss;	// packet loss rate

	float D = s*n*Pe*L*D1;

	// map distortion D to PSNR
	float psnr = 10*log10(255*255/D);

	// map PSNR to quality
	float b1 = 0.5, b2 = 30;
	float ql = 1 / ( 1 + exp(b1*(psnr-b2)) );

	// map quality to MOS = Quality*(-3.5)+4.5
	float mos = ql * 3.5 + 4.5;

//	// linear mapping from PSNR 20dB (MOS 1) to 40dB (MOS 5)
//	float mos = psnr*0.20 - 3;

//	fprintf(stderr, "\t[%d] videoMOS:\n"
//			"\t\tgamma %f T %f alpha %f\n"
//			"\t\tlost %zd mse_total %f sigma %f D1 %f\n"
//			"\t\ts %f L %f n %f Pe %f\n"
//			"\t\tD %f PSNR %f MOS %f\n",
//			mac_->nodeId(),
//			gamma, T, alpha,
//			mse->size(), mse_total, sigma, D1,
//			s, L, n, Pe,
//			D, psnr, mos);

	return mos;
}

float
WimshMOSScheduler::mseVideoMOS (float mse, unsigned int nlost, float loss)
{
	if(mse== 0)
		return 4.5;

	// define some parameters
	float gamma = 0; // attenuation factor
	float T = 10; // P-frames per GOP (10 in a 30-GOP)
	float alpha = ( pow(gamma, T+1) - (T+1)*gamma + T) / ( T*(1-gamma)*(1-gamma) ); // (5)

	// get sigma for (5), which is the average MSE of all dropped frames
	float sigma = mse / nlost;

	// obtain D1 (5)
	float D1 = alpha*sigma;

	// now obtain distortion D; we assume H.264 and D=s*n\*Pe*L*D1
	//	float s = 8; 	// 8 slices per frame, CIF
	float s = 1; 	// since we have the MSE per frame, we consider 1 frame=1 slice
	float L = 1; 	// 1 frame per packet, encoding
	float n = 1;	// bernoulli, each loss affects 1 packet
	float Pe = loss;	// packet loss rate

	float D = s*n*Pe*L*D1;

	// map distortion D to PSNR
	float psnr = 10*log10(255*255/D);

	// map PSNR to quality
	float b1 = 0.5, b2 = 30;
	float ql = 1 / ( 1 + exp(b1*(psnr-b2)) );

	// map quality to MOS = Quality*(-3.5)+4.5
	float mos = ql * 3.5 + 4.5;

//	// linear mapping from PSNR 20dB (MOS 1) to 40dB (MOS 5)
//	float mos = psnr*0.20 - 3;

	return mos;
}


float
WimshMOSScheduler::deltaVideoMOS (vector<float>* mse, vector<float>* dropdist, MOSFlowInfo* flowinfo)
{
	if(mse->size() == 0) {
		// calculate MOS for a single drop of a frame, then for all drops, and delta

		// get the smallest mse
		float lowestmse=FLT_MAX;
		for(unsigned i=0; i < dropdist->size(); i++)
		{
			if( (*dropdist)[i] < lowestmse )
				lowestmse = (*dropdist)[i];
		}

		// estimate video MOS for dropping that mse
		std::vector<float> tempmse;
		tempmse.push_back(lowestmse);
		float oldloss = (float)(1) / (float)(1 + flowinfo->count_);
		float oldMOS = videoMOS(&tempmse, oldloss);


		// now estimate video MOS for all the MSEs
		if(dropdist->size() > 1) // if there's more than one packet being dropped
		{
			float newloss = (float)(dropdist->size()) / (float)(dropdist->size() + flowinfo->count_);
			float newMOS = videoMOS(mse, newloss); // new MOS

			// 4.5 is the estimate for no packet losses, if dropping this frame doesn't create a MOS > 4.5, then use that as a reference
			if(newMOS < 4.5)
				return (newMOS - 4.5);

			return (newMOS - oldMOS);
		}
		else
		{
			// 4.5 is the estimate for no packet losses, if dropping this frame doesn't create a MOS > 4.5, then use that as a reference
			if(oldMOS < 4.5)
				return (oldMOS - 4.5);

			tempmse.push_back(lowestmse); // again
			float newloss = (float)(2) / (float)(2 + flowinfo->count_); // two of these packets
			float newMOS = videoMOS(&tempmse, newloss); // new MOS
			return (newMOS - oldMOS);
		}
	}

	float oldloss = flowinfo->loss_;
	float newloss = (float)(flowinfo->lostcount_ + dropdist->size()) /
					(float)(flowinfo->lostcount_ + dropdist->size() + flowinfo->count_);

	std::vector<float> totalmse;

	// push original
	for(unsigned i=0; i < mse->size(); i++)
		totalmse.push_back( (*mse)[i] );

	// combination drop increase
	for(unsigned i=0; i < dropdist->size(); i++)
		totalmse.push_back( (*dropdist)[i] );

	float oldMOS=videoMOS(mse, oldloss);
	float newMOS=videoMOS(&totalmse, newloss);

	return (newMOS - oldMOS);
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

	// scheduler code enters here
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
							fprintf (stderr, "\t\tVOD_DATA\tfid %d ndx %d id %d size %d\tdistortion %f\n",
									pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
									pdulist_[i][k][l]->size() ,vodinfo_->distortion());
						} else if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOIP_DATA) {
							fprintf (stderr, "\t\tVOIP_DATA\tfid %d ndx %d id %d size %d\n",
									pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
									pdulist_[i][k][l]->size());
						}
					} else {
							fprintf (stderr, "\t\tFTP_DATA\tfid %d ndx %d id %d size %d\n",
									pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
									pdulist_[i][k][l]->size());
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

	// print flow ID information
	for(unsigned int i=0; i < stats_.size(); i++)
	{
		string traffic;
		switch(stats_[i].traffic_) {
		case M_VOD: traffic="VOD "; break;
		case M_VOIP: traffic="VOIP"; break;
		default: traffic="DATA";
		}
		fprintf (stderr, "\t\tfid %d traffic %s mos %f\n", stats_[i].fid_, traffic.c_str(), stats_[i].mos_);
	}

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


	// get the number of packets in the buffers
	unsigned int npackets = 0;
	for(unsigned i=0; i < pdulist_.size(); i++)
//		for(unsigned j=0; j < pdulist_[i].size(); j++)
			for(unsigned k=0; k < pdulist_[i].size(); k++)
				npackets += pdulist_[i][k].size();

	fprintf (stderr, "\t%d packets in the buffers\n", npackets);

	if(npackets > 0 && ((float)buffusage/(float)sched_->maxBufSize()) > 0.70) // temp, remove when this routine is called on buffer overflow
	{

		/* here, we have:
		 * pdulist_[i][k][l] -> vector with all packets in the buffers
		 * npackets
		 * flowids_[i] -> vector of all flowIDs
		 * stats_[i] -> flowinfo stats
		 */

		// need to get all packet combinations that reduce the buffer by X
		// combinations represent binary values as to whether the packet has been chosen or not
		// combinations (2^npackets)
		long int ncombs = pow(2, npackets);
		fprintf (stderr, "\tevaluating %ld combinations for size match\n", ncombs);

		unsigned lbound = 1500;
		unsigned rbound = 2000;

		fprintf (stderr, "\t\tcombination matches for [%d,%d]:\n\t\tID:", lbound, rbound);
//		std::vector<bool> binComb(npackets, 0);
		std::vector<long> validCombs;
		for (long combID = 0; combID < ncombs; combID++)
		{
			std::vector<bool> binComb;
			dec2bin(combID, &binComb);

			// fill the binComb for missing zeroes
			for(unsigned k=binComb.size(); k<npackets; k++)
				binComb.push_back(0);

			unsigned int combsize = 0;
			unsigned int packetid = 0;
			for(unsigned i=0; i < pdulist_.size(); i++)
		//		for(unsigned j=0; j < pdulist_[i].size(); j++)
					for(unsigned k=0; k < pdulist_[i].size(); k++)
						for(unsigned l=0; l < pdulist_[i][k].size(); l++)
						{
							if(combsize > rbound) break; // premature break in order to speed up the process
							if(binComb[packetid] == TRUE)
								combsize += pdulist_[i][k][l]->size();
							packetid++;
						}

			if(combsize > lbound && combsize < rbound) {
				validCombs.push_back(combID);
				fprintf (stderr, " %ld", combID);
			}
		}
		fprintf (stderr, "\n");

		// here we have: validCombs, with all matching combinations
		// use dec2bin(combID, &binComb) to get packet references

		fprintf (stderr, "\t%zd combinations matched, processing...\n", validCombs.size());

		// process all combinations
		// vector to hold CombInfo elements
		std::vector<CombInfo> combstats_;
		for(unsigned p=0; p<validCombs.size(); p++)
		{
			fprintf (stderr, "\t\tcombID %ld:\n", validCombs[p]);

			std::vector<bool> binComb;
			dec2bin(validCombs[p], &binComb);
			for(unsigned k=binComb.size(); k<npackets; k++)
				binComb.push_back(0);
			unsigned int packetid = 0;
			std::vector<WimaxPdu*> combPdu;
			for(unsigned i=0; i < pdulist_.size(); i++)
//				for(unsigned j=0; j < pdulist_[i].size(); j++)
					for(unsigned k=0; k < pdulist_[i].size(); k++)
						for(unsigned l=0; l < pdulist_[i][k].size(); l++)
						{
							if(binComb[packetid] == TRUE)
							{
							// print the packet info
								if(pdulist_[i][k][l]->sdu()->ip()->datalen()) {
									if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOD_DATA) {
										VideoData* vodinfo_ = (VideoData*)pdulist_[i][k][l]->sdu()->ip()->userdata();
										fprintf (stderr, "\t\t\tVOD_DATA\tfid %d ndx %d id %d size %d\tdistortion %f\n",
												pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
												pdulist_[i][k][l]->size() ,vodinfo_->distortion());
									} else if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOIP_DATA) {
										fprintf (stderr, "\t\t\tVOIP_DATA\tfid %d ndx %d id %d size %d\n",
												pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
												pdulist_[i][k][l]->size());
									}
								} else {
										fprintf (stderr, "\t\t\tFTP_DATA\tfid %d ndx %d id %d size %d\n",
												pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
												pdulist_[i][k][l]->size());
								}
							// store the packet for MOS processing
								combPdu.push_back(pdulist_[i][k][l]);
							}
							packetid++;
						}


			// calculate MOS and impact here
			// vector combPdu is where it's at
			{
				// fill a vector with the flowIDs in this packet combination
				std::vector <int> combfIDs;
				for(unsigned l=0; l < combPdu.size(); l++)
				{
					// vector empty, push
					if(combfIDs.size() == 0) {
						combfIDs.push_back(combPdu[l]->sdu()->flowId());
					} else {
						// see if the fid is already in the list
						bool inlist_ = false;
						for(unsigned int m=0; m < combfIDs.size(); m++) {
							if(combfIDs[m] == combPdu[l]->sdu()->flowId()) {
								inlist_ = true;
							}
						}
						// not on the list, push
						if(!inlist_)
							combfIDs.push_back(combPdu[l]->sdu()->flowId());
					}
				}

				// for each flow, aggregate all of its packets and estimate MOS impact
				vector<float> combMOSdrop; // stores accumulated MOS variations for this combID
				for(unsigned l=0; l < combfIDs.size(); l++)
				{
					// get the FlowInfo of this flowID
					unsigned k;
					for(k=0; stats_[k].fid_ != combfIDs[l]; k++);

					if(stats_[k].traffic_ == M_VOD)
					{
						// VOD, sum the distortion of all packets and pass it to deltaVideoMOS
						vector<float> dropdist; unsigned nCombPackets=0;
						for(unsigned n=0; n<combPdu.size(); n++) // for each PDU
						{
							if(combPdu[n]->sdu()->flowId() == combfIDs[l]) // if it belongs to the current flowID
							{
								if(combPdu[n]->sdu()->ip()->userdata()->type() == VOD_DATA) { // and is VOD
									VideoData* vodinfo_ = (VideoData*)combPdu[n]->sdu()->ip()->userdata();
									dropdist.push_back(vodinfo_->distortion()); // accumulate distortion
									nCombPackets++;
								}
							}
						}

						// {fid, dist, packets}: {combfIDs[l], totalDist, nCombPackets}
						assert(nCombPackets>0);

						// nCombPackets, estimate drop percentage increase, get new MOS
						float deltaMOS = deltaVideoMOS(&(stats_[k].mse_), &dropdist, &(stats_[k]));

						// associate the deltaMOS to this combID for later processing
						combMOSdrop.push_back(deltaMOS);

						fprintf(stderr, "\t\t\t(-) VOD  fid %d drop %d delta %f\n",
								combfIDs[l], nCombPackets, deltaMOS);
					}
					else if(stats_[k].traffic_ == M_VOIP)
					{
						// VOIP, get the number of packets for error percentage
						unsigned nCombPackets=0;
						for(unsigned n=0; n<combPdu.size(); n++) // for each PDU
							if(combPdu[n]->sdu()->flowId() == combfIDs[l]) // if it belongs to the current flowID
								nCombPackets++;

						assert(nCombPackets>0);

						// nCombPackets, estimate drop percentage increase, get new MOS
						float newloss = (float)(stats_[k].lostcount_ + nCombPackets) /
										(float)(stats_[k].lostcount_ + nCombPackets + stats_[k].count_);

						float newMOS = audioMOS(stats_[k].delay_, newloss);
						float deltaMOS = newMOS - stats_[k].mos_;

						// associate the deltaMOS to this combID for later processing
						combMOSdrop.push_back(deltaMOS);

						fprintf(stderr, "\t\t\t(-) VOIP fid %d drop %d delta %f oldMOS %f newMOS %f\n",
								combfIDs[l], nCombPackets, deltaMOS, stats_[k].mos_, newMOS);

					}
					else
					{
						// FTP
						// TODO: nothing for now
					}

				} // end of for each fID of this comb


				// process list of distortion impacts, combMOSdrop[]
				{
					// total drop
					float totalMOSdrop = 0;
					for(unsigned i=0; i<combMOSdrop.size(); i++)
						totalMOSdrop+=combMOSdrop[i];

					// average drop
					float avgMOSdrop = totalMOSdrop / combMOSdrop.size();

					// standard deviation
					float stdMOSdrop = stddev(&combMOSdrop);

					// store info on vector combstats_
					CombInfo tempcombinfo(validCombs[p], totalMOSdrop, avgMOSdrop, stdMOSdrop);
					combstats_.push_back(tempcombinfo);

					fprintf(stderr, "\t\t\t(>) combID %ld flows impacted %zd MOSimpact total %f avg %f std %f\n",
							validCombs[p], combMOSdrop.size(), totalMOSdrop, avgMOSdrop, stdMOSdrop);
				}

			} // end of MOS calculations

		} // end of combination packet listing

		// algorithm to choose the best combID from combstats_ TODO
		long dropCombID = 0; // best combination for dropping
		float combImpact = -500; // valor

		float stddevcoeff = 1;

		for(unsigned r=0; r < combstats_.size(); r++)
		{
			float impact = combstats_[r].total_ - stddevcoeff*combstats_[r].std_;
			combstats_[r].impact_ = impact;

//			fprintf(stderr, "\tDEBUG impact %f total %f std %f\n",
//					impact, combstats_[r].total_, stddevcoeff*combstats_[r].std_);

			if(impact > combImpact) // we're working with negative values
			{
				combImpact = impact;
				dropCombID = combstats_[r].combID_;
			}

		}

		// sum up results
		fprintf(stderr, "\tsynopsis of combinations:\n");
		for(unsigned j=0; j<combstats_.size(); j++)
		{
			fprintf(stderr, "\t\tcombID %ld:\t total %f avg %f std %f impact %f\n",
					combstats_[j].combID_, combstats_[j].total_, combstats_[j].avg_, combstats_[j].std_,
					combstats_[j].impact_);
		}

		fprintf(stderr, "\tSelecting combID %ld for drop, impact %f\n", dropCombID, combImpact);

		// packet dropping of chosen combID
		std::vector<bool> binComb;
		dec2bin(dropCombID, &binComb);

		// fill the binComb for missing zeroes
		for(unsigned k=binComb.size(); k<npackets; k++)
			binComb.push_back(0);

		// kill the packets
		unsigned int packetid = 0;
		for(unsigned i=0; i < pdulist_.size(); i++)
	//		for(unsigned j=0; j < pdulist_[i].size(); j++)
				for(unsigned k=0; k < pdulist_[i].size(); k++)
				{
					vector<WimaxPdu*>::iterator iter1 = pdulist_[i][k].begin();
					while( iter1 != pdulist_[i][k].end())
					{
						if(binComb[packetid] == TRUE)
						{
							WimaxPdu* pdu = *iter1;
							fprintf (stderr, "\t\tDropping packet fid %d ndx %d id %d size %d\n",
									pdu->sdu()->flowId(), i, pdu->sdu()->seqnumber(), pdu->size());

							dropPDU(pdu);
//							pdu->sdu()->freePayload();
//							delete pdu->sdu();
//							delete pdu;
							// reduce buffer occupancy
							sched_->setBufSize(sched_->bufSize() - pdu->size());
							// erase packet from the pdu list
							pdulist_[i][k].erase(iter1);
						}
						packetid++;
						++iter1;
					 }
				}

	} // end of buffer size check


	// re-print the packet list, for tests
	for(unsigned i=0; i < pdulist_.size(); i++)
//		for(unsigned j=0; j < pdulist_[i].size(); j++)
			for(unsigned k=0; k < pdulist_[i].size(); k++)
				for(unsigned l=0; l < pdulist_[i][k].size(); l++) {
					if(pdulist_[i][k][l]->sdu()->ip()->datalen()) {
						if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOD_DATA) {
							VideoData* vodinfo_ = (VideoData*)pdulist_[i][k][l]->sdu()->ip()->userdata();
							fprintf (stderr, "\t\tVOD_DATA\tfid %d ndx %d id %d size %d\tdistortion %f\n",
									pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
									pdulist_[i][k][l]->size() ,vodinfo_->distortion());
						} else if(pdulist_[i][k][l]->sdu()->ip()->userdata()->type() == VOIP_DATA) {
							fprintf (stderr, "\t\tVOIP_DATA\tfid %d ndx %d id %d size %d\n",
									pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
									pdulist_[i][k][l]->size());
						}
					} else {
							fprintf (stderr, "\t\tFTP_DATA\tfid %d ndx %d id %d size %d\n",
									pdulist_[i][k][l]->sdu()->flowId(), i, pdulist_[i][k][l]->sdu()->seqnumber(),
									pdulist_[i][k][l]->size());
					}
				}


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

void
WimshMOSScheduler::dec2bin(long decimal, vector<bool>* binary)
{
	long remain;
	vector<bool> temp;
	do	{
		remain = decimal % 2;
		// whittle down the decimal number
		decimal = decimal / 2;
		// converts digit 0 or 1 to character '0' or '1'
		binary->push_back(remain);
	} while (decimal > 0);
}

float
WimshMOSScheduler::stddev(vector<float>* values)
{
	float total = 0, mean = 0, stddev = 0;

	// get the total
	for(unsigned i=0; i < values->size(); i++)
		total += (*values)[i];

	// get the mean
	mean = total / values->size();

	// get the squares
	float sqraccum = 0;
	for(unsigned i=0; i < values->size(); i++)
		sqraccum += pow( ((*values)[i] - mean), 2);

	// get the avg square
	float avgsqr = sqraccum / values->size();

	// get the stddev
	stddev = sqrt(avgsqr);

	return stddev;
}
