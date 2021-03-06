/*
 *  Copyright (C) 2007 Dip. Ing. dell'Informazione, University of Pisa, Italy
 *  http://info.iet.unipi.it/~cng/ns2mesh80216/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA, USA
 */

#include <wimsh_scheduler_frr.h>

#include <wimsh_mac.h>
#include <wimsh_packet.h>
#include <wimsh_bwmanager.h>
#include <wimsh_mossched.h>

#include <stat.h>
#include <ip.h>

#include <videodata.h>
#include <ns-process.h>

WimshSchedulerFairRR::WimshSchedulerFairRR (WimshMac* m) :
	WimshScheduler (m), timer_(this)
{
	roundDuration_     = 0;
	bufferSharingMode_ = SHARED;
	interval_          = 0;
	for ( unsigned int i = 0 ; i < WimaxMeshCid::MAX_PRIO ; i++ )
		prioWeights_[i] = 1.0;
}

int
WimshSchedulerFairRR::command (int argc, const char * const* argv)
{
	if ( argc == 2 && strcmp (argv[0], "round-duration") == 0 ) {
		if ( atoi (argv[1]) <= 0 ) {
			fprintf (stderr, "Invalid round duration '%s'. "
					"Choose a positive number, in bytes\n", argv[1]);
			return TCL_ERROR;
		}
		roundDuration_ = (unsigned int) atoi (argv[1]);
		return TCL_OK;
	} else if ( argc == 2 && strcmp (argv[0], "buffer-sharing") == 0 ) {
		if ( strcmp (argv[1], "shared") == 0 ) {
			bufferSharingMode_ = SHARED;
		} else if ( strcmp (argv[1], "per-link") == 0 ) {
			bufferSharingMode_ = PER_LINK;
		} else if ( strcmp (argv[1], "per-flow") == 0 ) {
			bufferSharingMode_ = PER_FLOW;
		} else {
			fprintf (stderr, "Invalid buffer sharing strategy '%s'. "
					"Choose 'shared' or 'per-link' or 'per-flow'\n", argv[1]);
			return TCL_ERROR;
		}
		return TCL_OK;
   } else if ( argc == 2 && strcmp (argv[0], "weight-timeout") == 0 ) {
      interval_ = atof (argv[1]);
      return TCL_OK;
   } else if ( argc == 3 && strcmp (argv[0], "prio-weight") == 0 ) {
      unsigned int i = atoi (argv[1]);
      double x = atof (argv[2]);
      if ( i >= WimaxMeshCid::MAX_PRIO ) {
         fprintf (stderr, "priority '%d' exceeds the maximum value '%d'\n",
               i, WimaxMeshCid::MAX_PRIO - 1 );
         return TCL_ERROR;
      }
      if ( x <= 0 ) {
         fprintf (stderr,
               "priority weight must be >= 0 ('%f' not allowed)\n", x);
         return TCL_ERROR;
      }
      prioWeights_[i] = x;
      return TCL_OK;
	}

	return WimshScheduler::command (argc, argv);
}

void
WimshSchedulerFairRR::initialize ()
{
	// resize the vector of links and unfinished vectors
	link_.resize (mac_->nneighs());
	unfinishedRound_.resize (mac_->nneighs());

	// start the timer to remove stale flows
	if ( interval_ > 0 ) timer_.start ( interval_ );
}

void
WimshSchedulerFairRR::addPdu (WimaxPdu* pdu)
{
	if ( WimaxDebug::trace("WSCH::addPdu") ) fprintf (stderr,
			"%.9f WSCH::addPdu     [%d] %s seq %d\n",
			NOW, mac_->nodeId(), WimaxDebug::format(pdu), pdu->sdu()->seqnumber());

	// index used to identify the next-hop neighbor node
	const unsigned int ndx = mac_->neigh2ndx (pdu->hdr().meshCid().dst());

	// priority of this PDU
	const unsigned char prio = pdu->hdr().meshCid().priority();

	// call the RD scheduler on buffer overflow
	if(bufSize_ + pdu->size() > maxBufSize_)
		mac_->mosscheduler()->trigger(pdu->size());

	// if the size of this PDU overflows the buffer size, drop the PDU/SDU/IP
	if ( ( bufferSharingMode_ == SHARED && bufSize_ + pdu->size() > maxBufSize_ ) ||
		  ( bufferSharingMode_ == PER_LINK && link_[ndx].size_ + pdu->size() > maxBufSize_) )
	{
		// stat the lost packet for the RD sched
		Stat::put ("rd_packet_lost_overflow", pdu->sdu()->flowId(), 1);
		// notify the MAC layer of the pdu about to be dropped
		mac_->dropPDU(pdu);
		drop (pdu);
		return;
	}

	// check if there are already PDUs buffered from this traffic flow
	// to do so, we create a new FlowDesc with empty queue, which is
	// used only to this purpose. If there are no buffered PDUs that
	// belong to this traffic flow, then this descriptor is added
	// to the round-robin list

	// get the source/destination addresses of this PDU
	const WimaxNodeId src = (WimaxNodeId) HDR_IP(pdu->sdu()->ip())->saddr();
	const WimaxNodeId dst = (WimaxNodeId) HDR_IP(pdu->sdu()->ip())->daddr();

	// temporary flow descriptor used to search into the active-list
	FlowDesc tmp (src, dst, prio);

	// flag used to check whether such a flow already exists
	bool valid;
	FlowDesc& desc = link_[ndx].rr_.find (tmp, valid);

	// check for per-flow buffer overflow
	if ( bufferSharingMode_ == PER_FLOW &&
		  desc.size_ + pdu->size() > maxBufSize_ ) {
		drop (pdu);
		return;
	}

	// add the PDU to the flow descriptor
	desc.queue_.push (pdu);

	// call the RD scheduler
//	RDscheduler(pdu, ndx);

	// update the cumulative/flow/link buffer occupancy
	// this includes MAC overhead (header/crc), but not the fragmentation
	// subheader. In fact, at the moment we do not know whether it will
	// be added by the fragmentation buffer or not

	desc.size_ += pdu->size();        // per flow
	link_[ndx].size_ += pdu->size();  // per link
	bufSize_ += pdu->size();          // shared

	Stat::put ("wimsh_bufsize_mac_a", mac_->index(), bufSize_ );
	Stat::put ("wimsh_bufsize_mac_d", mac_->index(), bufSize_ );

	// add the flow descriptor into the active list, is not already there
	// if this is the case, then the weights should be recomputed, as well
	if ( ! valid ) {
		link_[ndx].rr_.insert (desc);
		recompute (link_[ndx].rr_);
	}

	// indicate the updated backlog to the bandwidth manager
	mac_->bwmanager()->backlog (
			desc.src_,                        // src node
			desc.dst_,                        // dst node
			pdu->hdr().meshCid().priority(),  // priority
			pdu->hdr().meshCid().dst(),       // next-hop
		   pdu->size());                     // bytes
}

void
WimshSchedulerFairRR::schedule (WimshFragmentationBuffer& frag, WimaxNodeId dst)
{
	// get the link index
	unsigned int ndx = mac_->neigh2ndx(dst);

	// print some debug information
	if ( WimaxDebug::trace("WSCH::schedule") ) {
		fprintf (stderr,
			"%.9f WSCH::schedule   [%d] dst %d remaining %d",
			NOW, mac_->nodeId(), dst, frag.size());

		CircularList<FlowDesc>& rr = link_[ndx].rr_;
		fprintf (stderr, " backlog %d", link_[ndx].size_);
		if ( rr.size() == 0 ) {
			fprintf (stderr, " qsize 0");
		} else {
			fprintf (stderr, " qsize %d", rr.size());
			for ( unsigned int i = 0 ; i < rr.size() ; i++ ) {
				FlowDesc& flow = rr.current();
				fprintf (stderr, " (%d,%d,%d; %d; %d,%d)",
						flow.src_, flow.dst_, flow.prio_,
						flow.size_,
						flow.quantum_, flow.deficit_);
				rr.move ();
			}
		}

		fprintf (stderr, "\n");
	}

	// true until there is spare room into the burst
	bool spare = true;

	// if there is an unfinished round, then the latter must be terminated
	// before new flows get scheduled, even though the priority of the
	// current flow is lower than that of other active flows
	// however, then we do not update the deficit counter with a new quantum

	if ( unfinishedRound_[ndx] ) spare = serve (frag, ndx, true);

	// get the round-robin list
	CircularList<FlowDesc>& rr = link_[ndx].rr_;

	// schedule PDUs directed to the specified neighbor on a DRR fashion
	// until there is room into the fragmentation buffer
	// and the output queue towards that node is backlogged
	while ( ! rr.empty() && spare ) spare = serve (frag, ndx, false);
}

void
WimshSchedulerFairRR::drop (WimaxPdu* pdu)
{
	pdu->sdu()->freePayload();
	delete pdu->sdu();
	delete pdu;

	Stat::put ("wimsh_drop_overflow", mac_->index(), 1.0);
}

bool
WimshSchedulerFairRR::serve (WimshFragmentationBuffer& frag,
		unsigned int ndx, bool unfinished)
{
	// get the current round-robin list
	CircularList<FlowDesc>& rr = link_[ndx].rr_;

	// get the current traffic flow
	FlowDesc& flow = rr.current();

	// update the deficit counter, unless unfinished
	if ( ! unfinished ) {
		flow.deficit_ += flow.quantum_;
		flow.deficit_ =
			( flow.deficit_ > flow.size_ ) ? flow.size_ : flow.deficit_;
	}

	// becomes false if there is not anymore capacity available
	bool spare = true;

	// try to grant up to flow.deficit_ bytes from this queue
	while ( spare && flow.deficit_ > 0 ) {

		// get the head-of-line PDU
		WimaxPdu* pdu = flow.queue_.front ();

		// check whether this PDU fits into the remaining deficit
		// if not, then it is not possible to transmit it
		if ( pdu->size() > flow.deficit_ ) break;

		// otherwise, remove the PDU from queue
		// note that it may happen that this PDU will not be entirely
		// transmitted by the fragmentation buffer. This is not an issue
		// from the point of view of the scheduler, since the MAC is
		// responsible for ensuring that pending fragments are transmitted
		// in order before other PDUs have any chance to interfere
		flow.queue_.pop ();

		// update the buffer occupancies
		flow.size_ -= pdu->size();          // flow
		link_[ndx].size_ -= pdu->size();    // link

		if(bufSize_ < pdu->size())	// TODO: this is a hack, something's missing
			bufSize_ = 0;
		else
			bufSize_ -= pdu->size();            // MAC

		Stat::put ("wimsh_bufsize_mac_a", mac_->index(), bufSize_ );
		Stat::put ("wimsh_bufsize_mac_d", mac_->index(), bufSize_ );

		// update the deficit counter
////////////////////////////////////////////////////////////////////////////////
// :XXX: -= sdu->size()
		flow.deficit_ -= pdu->size();
////////////////////////////////////////////////////////////////////////////////

		// add the PDU to the fragmentation buffer
		spare = frag.addPdu (pdu);
	}

	// if this round terminated because there is no more spare room
	// in the burst of PDUs, then we neither move the round-robin pointer
	// nor we remove the current element. In this way, next scheduling
	// phase will restart exactly from where we leave now

	if ( ! spare ) {
		unfinishedRound_[ndx] = true;
		return spare;
	}

	unfinishedRound_[ndx] = false;

	// if the current element does not have any more PDUs, remove it
	if ( flow.size_ == 0 )
		rr.erase ();

	// else, move the round-robin pointer to the next element in the list
	else
		rr.move ();

	return spare;
}

void
WimshSchedulerFairRR::recompute (CircularList<FlowDesc>& rr)
{
	// Number of elements in the round-robin list.
	const unsigned int N = rr.size();
	// Sum of the priority weights of all the elements.
	double sum = 0;

	// Compute the sum of the priority weights.
	for ( unsigned int i = 0 ; i < N ; i++ ) {
		sum += prioWeights_[rr.current().prio_];
		rr.move ();
	}

	// Update the quanta. The minimum quantum value is 1 byte.
	for ( unsigned int i = 0 ; i < N ; i++ ) {
		unsigned int quantum = (unsigned int) ( roundDuration_ *
				( prioWeights_[rr.current().prio_] / (double)sum ) );
		rr.current().quantum_ = ( quantum > 0 ) ? quantum : 1;
		rr.move ();
	}
}

void
WimshSchedulerFairRR::handle ()
{
	if ( WimaxDebug::trace("WSCH::handle") ) fprintf (stderr,
			"%.9f WSCH::handle     [%d]\n", NOW, mac_->nodeId());

	// if the user specified a zero or negative value for the interval
	// then he/she wants the weights never to expire => return immediately
	// without restarting the timer, hence this function will not be
	// called again
	if ( interval_ <= 0 ) return;

	// remove stale flows
	std::vector<LinkDesc>::iterator link;
	for ( link = link_.begin() ; link != link_.end() ; ++link ) {
		CircularList<FlowDesc>& rr = link->rr_;

		bool changed = false;

		for ( unsigned int i = 0 ; i < rr.size() ; i++ ) {
		   // get the current traffic flow
		   FlowDesc& flow = rr.current();

			// remove the traffic flow if it has to be considered as inactive, i.e.
			// there are not any enqueued PDUs and the user-specified flow timeout
			// interval has expired since the last PDU has been received
			if ( NOW - flow.last_ > interval_ && flow.queue_.empty() ) {
				rr.erase ();
				changed = true;
			} else {
				rr.move ();
			}
		}
		// recompute weights, if needed
		if ( changed ) recompute (rr);
	}

	// restart the timer
	timer_.start ( interval_ );
}

void
WimshSchedulerFairRR::RDscheduler (WimaxPdu* pdu, unsigned int ndx) {

	if ( WimaxDebug::trace("WSCH::RDscheduler") ) fprintf (stderr,
			"%.9f WSCH::RDscheduler[%d]\n", NOW, mac_->nodeId());

	// this scheduler won't work with per-flow or per-link buffer sharing
	if (bufferSharingMode_ == PER_FLOW || bufferSharingMode_ == PER_LINK) return;

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

	if(pdu->sdu()->ip()->datalen())
		if(pdu->sdu()->ip()->userdata()->type() == VOD_DATA) {
			VideoData* vodinfo_ = (VideoData*)pdu->sdu()->ip()->userdata();
			fprintf (stderr, "\tGot a VOD_DATA packet, distortion %f\n", vodinfo_->distortion());
		}

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
			std::list<FlowDesc> list_ = link_[i].rr_.list();
			list<FlowDesc>::iterator iter1 = list_.begin();

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
	fprintf (stderr, "\tBuffer usage %d/%d, %.2f%%\n", buffusage, maxBufSize_, ((float)buffusage/(float)maxBufSize_)*100 );

	// reconstruct queues
	for(unsigned i=0; i < mac_->nneighs(); i++) {
		// run all services
//		for(unsigned j=0; j < wimax::N_SERV_CLASS; j++) {
			std::list<FlowDesc> list_ = link_[i].rr_.list();
			list<FlowDesc>::iterator iter1 = list_.begin();

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
