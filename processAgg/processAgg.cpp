#include "stdafx.h"
#include "Point.h"
#include "Rectangle.h"

#include "SFCConversion.h"
//#include "OutputSchema.h"
#include "QueryBySFC.h"

#include "SFCPipeline.h"
#include "SFCConversion2.h"
#include "OutputSchema2.h"

#include "RandomLOD.h"

#include "hiredis/hircluster.h"

#include "tbb/task_scheduler_init.h"

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


void print_ranges(char * str, vector<sfc_bigint>& ranges)
{
	sfc_bigint ntotal_len = 0;
	if (str == NULL) return;

	cout << str << endl;
	for (int i = 0; i < ranges.size(); i = i + 2)
	{
		//printf("\n");

		//printf("%lld---%lld\n", ranges[i], ranges[i + 1]);
		//cout << ranges[i] << "----" << ranges[i + 1] <<endl;

		ntotal_len += (ranges[i + 1] - ranges[i] + 1);
	}

	cout << "total ranges len:  " << ntotal_len << endl;
}

template<int nDims, int mBits>
struct ptdecode
{
	vector<sfc_bigint>& _vec_pts;
	Point<double, nDims>* _vec_outpts;
	CoordTransform<double, long, nDims>& _cotrans;

	int _sfctype;

	ptdecode(vector<sfc_bigint>& _pts, Point<double, nDims>* _outpts, CoordTransform<double, long, nDims>& cotrans_, int type_) : \
		_vec_pts(_pts), _vec_outpts(_outpts), _cotrans(cotrans_), _sfctype(type_) {}

	void operator( )(const blocked_range<size_t> range)const
	{
		SFCConversion<nDims, mBits> sfcgen;

		for (size_t i = range.begin(); i != range.end(); ++i)
		{
			if (_sfctype == 0)
				_vec_outpts[i] = _cotrans.InverseTransform(sfcgen.MortonDecode(_vec_pts[i]));
			else
				_vec_outpts[i] = _cotrans.InverseTransform(sfcgen.HilbertDecode(_vec_pts[i]));
		}//end for
	}//end functioner
};


int integerHash(double delta[], double x, double y)
{
	return (x - delta[2]) * 1000000000 + (y - delta[1]) * 10000;

}

string antiHash(double delta[], int key)
{
	double val[2] = { 0 };
	val[0] = (key / 100000) / 10000. + delta[2];
	val[1] = (key % 100000) / 10000. + delta[1];


	ostringstream oss;
	oss << "{\"x\":" << val[0] << ",\"y\":" << val[1] << ",\"v\":";
	string res = oss.str();

	return res;
}

int main(int argc, char* argv[])
{
	const int ndims = 3;
	const int mbits = 13;

	int nsfc_type = 0;
	int nencode_type = 0;
	///////////////////////////////////////////////////////////////////////
	//here the SFCQuery tool
	//-i 347068810/347068850/-73.96/-73.91/40.5/41/-73.99/-73.90/40.5/41 -s 1 -e 0 -t ./cttaxi.txt -n 2000 -k 4 -o range.sql
	//85999.42,446266.47,-1.65,9,651295384353375995169439
	//-i 85999.0/85999.5/446266/446266.5/-2.0/-1.5/8/9 -s 1 -e 0 -t ct.txt -n 1000 -o qq3.sql
	//85999.1,446250.23,-1.69,9,651295397912973650169147
	//-i 85999.0/85999.5/446250/446250.4/-2.0/-1.5/8/9 -s 1 -e 0 -t ct.txt -n 0 -o qq5.sql
	//-i 85545.3000/85695.3000/446465.6500/446615.6500/-99999999.0000/-99999999.0000/-99999999.0000/-99999999.0000 -s 1 -e 0 -t ct.txt -n 0 -o qq5.sql
	//-i 85545.3000/85695.3000/446465.6500/446615.6500/-2.0000/-1.5000/8.0000/9.0000 -s 1 -e 0 -t ..\SFCLib\ct.txt -n 5000 -o qq5.sql -v -p 1
	bool bstat = false; //control statistics

	int nparallel = 1;
	int nranges = 0; //if nranges =0; means search to the bottom level
	int ktimes = ndims; // the ktimes* nranges is used to control tree traversal depth

	char szinput[1024] = { 0 };//1.xyz
	char szoutput[256] = { 0 };
	char sztransfile[256] = { 0 };
	char lvl[10] = { 0 };

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-i") == 0)//input filter coordinates
		{
			i++;
			strcpy(szinput, argv[i]);
			continue;
		}

		if (strcmp(argv[i], "-o") == 0)//output file path
		{
			i++;
			strcpy(szoutput, argv[i]);
			continue;
		}

		if (strcmp(argv[i], "-s") == 0)//sfc conversion type: 0 morthon, 1 hilbert
		{
			i++;
			nsfc_type = atoi(argv[i]);
			continue;
		}

		if (strcmp(argv[i], "-e") == 0)//output encoding type: 0 number 1 base32 2 base64
		{
			i++;
			nencode_type = atoi(argv[i]);
			continue;
		}

		if (strcmp(argv[i], "-t") == 0)//coordinates transformation file, two lines: translation and scale, comma separated
		{
			i++;
			strcpy(sztransfile, argv[i]);
			continue;
		}
		if (strcmp(argv[i], "-n") == 0)//number of return ranges
		{
			i++;
			nranges = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "-k") == 0)//k times of returned ranes for gap merge
		{
			i++;
			ktimes = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "-v") == 0)//k times of returned ranes for gap merge
		{
			i++;
			bstat = true;
			continue;
		}
		if (strcmp(argv[i], "-p") == 0)//if parallel: 0 sequential, 1 max parallel
		{
			i++;
			nparallel = atoi(argv[i]);
			continue;
		}
		if (strcmp(argv[i], "-l") == 0)//give the map level
		{
			i++;
			strcpy(lvl, argv[i]);
			continue;
		}
	}

	///////////////////////////////////////////////////
	////get the coordinates transfomration file--one more for lod value
	double delta[ndims + 1] = { 0 }; // 526000, 4333000, 300
	double  scale[ndims + 1] = { 1 }; //100, 100, 1000

	for (int i = 1; i < ndims + 1; i++)
	{
		delta[i] = 0;
		scale[i] = 1;
	}

	if (strlen(sztransfile) != 0)
	{
		FILE* input_file = NULL;
		input_file = fopen(sztransfile, "r");
		if (input_file)
		{
			int j;
			char buf[1024];
			char * pch, *lastpos;
			char ele[128];

			//////translation
			memset(buf, 0, 1024);
			fgets(buf, 1024, input_file);

			j = 0;
			lastpos = buf;
			pch = strchr(buf, ',');
			while (pch != NULL)
			{
				memset(ele, 0, 128);
				strncpy(ele, lastpos, pch - lastpos);
				//printf("found at %d\n", pch - str + 1);
				delta[j] = atof(ele);
				j++;

				lastpos = pch + 1;
				pch = strchr(lastpos, ',');
			}
			delta[j] = atof(lastpos); //final part

									  //////scale
			memset(buf, 0, 1024);
			fgets(buf, 1024, input_file);

			j = 0;
			lastpos = buf;
			pch = strchr(buf, ',');
			while (pch != NULL)
			{
				memset(ele, 0, 128);
				strncpy(ele, lastpos, pch - lastpos);
				//printf("found at %d\n", pch - str + 1);
				scale[j] = atof(ele);
				j++;

				lastpos = pch + 1;
				pch = strchr(lastpos, ',');
			}
			scale[j] = atof(lastpos); //final part

			fclose(input_file);
		}//end if input_file
	}//end if strlen

	CoordTransform<double, long, ndims> cotrans;
	cotrans.SetTransform(delta, scale);
	////////////////////////////////////////////////
	//get the input filter
	double pt1[ndims] = { 0.0f };
	double pt2[ndims] = { 0.0f };

	unsigned int dim_valid[ndims] = { 0 };

	memset(pt1, 0, sizeof(double)*ndims);
	memset(pt2, 0, sizeof(double)*ndims);
	memset(dim_valid, 0, sizeof(unsigned int)*ndims);

	char * pch, *lastpos;
	char ele[128];

	lastpos = szinput;
	for (int i = 0; i < ndims; i++)
	{
		///////min
		memset(ele, 0, 128);

		pch = strchr(lastpos, '//');
		strncpy(ele, lastpos, pch - lastpos);

		if (strcmp(ele, "-99999999.0000") != 0)//if "-99999999.0000", not set
		{
			pt1[i] = atof(ele);
			dim_valid[i] = 1;
		}
		else
		{
			pt1[i] = 0; ///this min value is not set,just assign 0
		}

		lastpos = pch + 1;
		///////max
		if (i != ndims - 1)
		{
			memset(ele, 0, 128);

			pch = strchr(lastpos, '//');
			strncpy(ele, lastpos, pch - lastpos);

			if (strcmp(ele, "-99999999.0000") != 0) //if "-99999999.0000", not set
			{
				pt2[i] = atof(ele);
				dim_valid[i] = 1;
			}
			else
			{
				pt2[i] = 0; ///this max value is not set,just assign 2^mbits -11 <  < mbits - 1
			}

			lastpos = pch + 1;
		}
		else
		{
			if (strcmp(lastpos, "-99999999.0000") != 0) //if "-99999999.0000", not set
			{
				pt2[i] = atof(lastpos);
				dim_valid[i] = 1;
			}
			else
			{
				pt2[i] = 0; ///this max value is not set,just assign 2^mbits -1 1 << mbits - 1
			}
		}

	}
	///////////////////////////////////////////////
	//point transfomration
	Point<double, ndims> MinPt1(pt1);
	Point<double, ndims> MaxPt1(pt2);

	Point<long, ndims> MinPt2 = cotrans.Transform(MinPt1);
	Point<long, ndims> MaxPt2 = cotrans.Transform(MaxPt1);

	///to check if any dim is not set
	for (int i = 0; i < ndims; i++)
	{
		if (dim_valid[i] == 0)// this dim is not set
		{
			MinPt2[i] = 0;
			MaxPt2[i] = 1 << mbits - 1;
		}
	}

	/////////////////////////////////////////////////////
	////query
	Rect<long, ndims> rec(MinPt2, MaxPt2);
	QueryBySFC<long, ndims, mbits> querytest;

	tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());

	tbb::tick_count t0 = tbb::tick_count::now();

	vector<sfc_bigint> vec_res2;

	//////output
	if (nparallel == 0)
		vec_res2 = querytest.RangeQueryByRecursive_LNG(rec, (SFCType)nsfc_type, nranges, ktimes);
	else
		vec_res2 = querytest.RangeQueryByRecursive_LNG_P(rec, (SFCType)nsfc_type, nranges, ktimes);

	tbb::tick_count t1 = tbb::tick_count::now();

	if (bstat)
		cout << "ranges time = " << (t1 - t0).seconds() << endl;

	///////////////////////////////////////////////////////////////////////
	//////here the redis query tool
	//char lvl[10] = { 0 };
	tbb::tick_count t3 = tbb::tick_count::now();

	/*redisClusterContext *cc = redisClusterContextInit();
	redisClusterSetOptionAddNodes(cc, "192.168.0.20:7006,192.168.0.27:7000,192.168.0.28:7001,192.168.0.29:7002,192.168.0.31:7004,192.168.0.32:7005");

	redisClusterConnect2(cc);*/

	redisClusterContext *cc = redisClusterConnect("192.168.0.20:7006,192.168.0.27:7000,192.168.0.28:7001,192.168.0.29:7002,192.168.0.31:7004,192.168.0.32:7005", HIRCLUSTER_FLAG_NULL);
	if (cc != NULL && cc->err)
	{
		printf("Error: %s\n", cc->errstr);
		// handle error
	}

	int count = 0;

	for (int i = 0; i < vec_res2.size(); i = i + 2)
	{
		redisClusterAppendCommand(cc, "zrangebyscore %s %s %s", lvl, vec_res2[i].str().c_str(), vec_res2[i + 1].str().c_str()); ///r1,r2

		count++;
	}

	char sfc_c[128];
	char sfc_v[128];
	redisReply * reply;

	vector<sfc_bigint> vec_sfc;
	vector<int> vec_cnt;

	vec_sfc.reserve(count);
	vec_cnt.reserve(count);

	for (int i = 0; i < count; i++)
	{
		redisClusterGetReply(cc, (void **)&reply);
		for (int j = 0; j < reply->elements; j++)
		{
			char* pch = strchr(reply->element[j]->str, ',');

			memset(sfc_c, 0, 128);
			memset(sfc_v, 0, 128);

			strncpy(sfc_c, reply->element[j]->str, pch - reply->element[j]->str);
			strcpy(sfc_v, pch + 1);
			sfc_v[strlen(sfc_v) - 1] = '\0';

			vec_sfc.push_back(sfc_bigint(sfc_v));
			vec_cnt.push_back(atoi(sfc_c));
		}
		freeReplyObject(reply);
	}
	redisClusterReset(cc);

	redisClusterFree(cc);

	tbb::tick_count t4 = tbb::tick_count::now();

	if (bstat)
		cout << "redis query time = " << (t4 - t3).seconds() << endl;

	///////////////////////////////////////////////////////////////////////
	////here is the parallel decoding
	tbb::tick_count t5 = tbb::tick_count::now();

	long pts_size = vec_sfc.size();
	Point<double, ndims>* p_outpts = new Point<double, ndims>[pts_size];

	parallel_for(blocked_range<size_t>(0, pts_size, 500), ptdecode<ndims, mbits>(vec_sfc, p_outpts, cotrans, nsfc_type));


	tbb::tick_count t6 = tbb::tick_count::now();

	if (bstat)
		cout << "decode time = " << (t6 - t5).seconds() << endl;
	////////////////////////////////////////////////////////////////////////
	////here is the merge step
	tbb::tick_count t7 = tbb::tick_count::now();

	std::map<int, int> map;
	std::map<int, int>::iterator it;

	for (long i = 0; i < pts_size; i++)
	{
		/*ostringstream oss;
		oss << p_outpts[i][ndims - 1] << ",\"y\":" << p_outpts[i][ndims - 2];
		string key = oss.str();*/

		int key = integerHash(delta, p_outpts[i][ndims - 1], p_outpts[i][ndims - 2]);

		it = map.find(key);
		if (it == map.end()) {
			map[key] = vec_cnt.at(i);
		}
		else {
			map[key] += vec_cnt.at(i);
		}
	}

	for (it = map.begin(); it != map.end(); ++it)
		cout << antiHash(delta, it->first) << it->second << "},";


	tbb::tick_count t8 = tbb::tick_count::now();

	if (bstat)
		cout << "merge time = " << (t8 - t7).seconds() << endl;

	////////////////////////////////
	delete[] p_outpts;

	return 0;

}