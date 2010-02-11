/*
 * Copyright (C) 2009 Andre Braga Reis
 * andrebragareis at gmail dot com
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

// TODO: doc

#ifndef ns_videodata_h
#define ns_videodata_h

class VideoData : public AppData {
private:
	float dist_;	// distortion importance
//	char ftype_;	// frame type
public:
//	VideoData(float distortion, char type) : AppData(VOD_DATA), dist_(0), ftype_(0)
	VideoData() : AppData(VOD_DATA), dist_(0) {};
	VideoData(float distortion) : AppData(VOD_DATA)
		{ dist_ = distortion; }
	~VideoData() {}
	float distortion() { return dist_; }
	int size() const { return sizeof(VideoData); }
	AppData* copy() {
		AppData *dup = new VideoData(dist_);
		return dup;
	}
	VideoData& operator= (const VideoData &vodinfo) {
		if(&vodinfo != this){
			dist_ = vodinfo.dist_;
//			ftype_ = vodinfo.ftype_;
		}
		return *this;
	}
};

#endif
