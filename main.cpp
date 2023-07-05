
/*
 * Peng Jiang: pengj@alumni.princeton.edu, May 19, 2023
 * The current version only accept 10x data with regular layout ST detection spots
 * */

extern "C" {
#include "MoranI.h"
}

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

// tool function to read a position value and check its format correctness
template <class T>
T load_positive_value(const string s, const string fieldname, const T minvalue, const T maxvalue)
{
	istringstream iss(s);
	T n;

	if(s[0] == '-') {
		cerr << "Please input a positive value for field \"" << fieldname << "\"." << endl;
		cerr << s << " is your input." << endl;
		exit(1);
	}

	iss >> n;

	if(iss.fail())
	{
		cerr << "Please input a valid value for field \"" << fieldname << "\"." << endl;
		cerr << s << " is your input." << endl;
		exit(1);
	}

	if(n < minvalue) {
		cerr << "The minimum value for field \"" << fieldname << "\" is " << minvalue << '.' << '\n';
		cerr << s << " is your input." << endl;
		exit(1);
	}

	if(n > maxvalue) {
		cerr << "The maximum value for field \"" << fieldname << "\" is " << maxvalue << '.' << '\n';
		cerr << s << " is your input." << endl;
		exit(1);
	}

	return n;
}


// read data from R like matrix file
double *read_matrix(const string &file, vector<string> &colnames, vector<string> &rownames)
{
	// temporary variables for parsing
	uint i, counter;
	string line;

	double *data;

	// read in matrix to temporary space
	ifstream fin(file.c_str());

	if(fin.fail())
	{
		cerr << "Error: Cannot open \"" << file << "\"." << endl;
		exit(1);
	}

	getline(fin, line, '\n');

	istringstream iss(line);

	// load column names
	for(;getline(iss,line,'\t');){
		colnames.push_back(line);
	}

	iss.clear();

	if(colnames.size()==0)
	{
		cerr << "Error: No columns in matrix \"" << file << "\"." << endl;
		exit(1);
	}

	// read row names, and get matrix size
	for(i=0; getline(fin,line,'\n'); i++)
	{
		// jump the last '\r'
		line.erase(line.find_last_not_of("\r") + 1);

		iss.str(line);
		getline(iss,line,'\t');

		if(iss.fail())
		{
			cerr << "Error: Empty line " << i << " in \"" << file << "\"." << endl;
			exit(1);
		}

		rownames.push_back(line);
		iss.clear();
	}

	fin.close();

	cout << colnames.size() << " columns\t" << rownames.size() << " rows" << endl;

	// total data size
	counter = colnames.size() * rownames.size();

	// allocate row major double array space
	data = new double[counter];
	memset(data, 0, counter * sizeof(double));

	// second round reading of float numbers
	fin.open(file.c_str());

	getline(fin, line, '\n');

	for(counter=0; getline(fin,line,'\n');)
	{
		// jump the last '\r'
		line.erase(line.find_last_not_of("\r") + 1);

		iss.str(line);
		getline(iss,line,'\t');

		for(i=0; i<colnames.size(); i++, counter++){
			iss >> data[counter];
		}

		iss.clear();
	}

	fin.close();

	return data;
}


int main(int argc, char *argv[])
{
	string input, output, value, type;
	vector<string> gene_names, spots;
	double *data, *distance, val;

	istringstream iss;

	uint i, j, n, ncol, max_radius = 5, platform = VISIUM, parseCnt = (argc - 1)/2;	// filter over rare Y out comes
	uint *spot_x, *spot_y, max_spot_x, max_spot_y, *spot_index_map;

	// for progress bar
	uint step, total, counter, nstep;
	
	// Whether calculate the 1st gene only
	bool allGenes_1stGene = true;

	// whether consider the same spot for computation
	bool flag_samespot = true;

	////////////////////////////////////////////////////////////////////////////////////////
	// Part 0: parameter input and check
	if ( argc < 5 )
	{
		if(argc == 2 && (value = argv[1]) == "-h")
		{
			cout << "\nCompute pairwise moran between all genes in 10x Spatial Transcriptomics\n" << endl;
			cout << "Usage: pairwise_moran_I -i input -o output [OPTIONS]\n"<<endl;

			cout << "Options:" <<endl;

			cout << "\t-r\t\tMaximum grid radius for computing distances. Default: "<< max_radius << endl;
			cout << "\t-p\t\t10x platform (" << VISIUM << ": VISIUM, " << OLD << ": OLD ST). Default: " << platform <<endl;
			cout << "\t-g\t\tWhether calculate all genes (1) or only 1st gene (0). Default: " << allGenes_1stGene << endl;
			cout << "\t-s\t\tWhether consider the same spot. Default: " << flag_samespot << endl;

			cout << "\nReport bugs to pengj@alumni.princeton.edu\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"pairwise_moran_I -h\" for help."<<endl;
			exit(1);
		}
	}

	// read in all parameters
	for(i=0;i<parseCnt;i++)
	{
		type = argv[2*i+1];
		value = argv[2*i+2];

		if(type == "-i"){
			input = value;

		}else if (type == "-o"){
			output = value;

		}else if (type == "-p"){
			platform = load_positive_value(value, type, (uint)0, (uint)10);

		}else if (type == "-r"){
			max_radius = load_positive_value(value, type, (uint)0, (uint)100);
		
		}else if (type == "-g"){
			allGenes_1stGene = (value[0]!='0');
		
		}else if (type == "-s"){
			flag_samespot = (value[0]!='0');

		}else if (type == "-h"){
			cerr << "Please don't use \"-h\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	if(input.empty())
	{
		cerr << "Cannot input file." << endl;
		exit(1);
	}

	if(output.empty())
	{
		cerr << "Cannot output file." << endl;
		exit(1);
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// Part 1: normalize data matrix, and create column spot ID index
	data = read_matrix(input, spots, gene_names);

	z_normalize(data, gene_names.size(), spots.size());

	// parse spot index
	spot_x = new uint[spots.size()];
	spot_y = new uint[spots.size()];

	for(i=0, max_spot_x=max_spot_y=0; i<spots.size(); i++)
	{
		iss.str(spots[i]);

		getline(iss, value, 'x');
		spot_x[i] = atoi(value.c_str());

		getline(iss, value, 'x');
		spot_y[i] = atoi(value.c_str());

		max_spot_x = MAX(max_spot_x, spot_x[i]);
		max_spot_y = MAX(max_spot_y, spot_y[i]);

		iss.clear();
	}

	// make sure all spots are sorted, later algorithm will need this
	for(i=1; i<spots.size(); i++)
	{
		assert(spot_x[i-1] <= spot_x[i]);

		if(spot_x[i-1] == spot_x[i]){
			assert(spot_y[i-1] <= spot_y[i]);
		}
	}

	// assume zero index, turn the index to max size
	max_spot_x++;
	max_spot_y++;

	// create map: spot x,y to i in data column
	n = max_spot_x * max_spot_y;
	spot_index_map = new uint[n];

	// means out of bound, position not exist
	for(i=0; i<n; i++){spot_index_map[i] = UINT_MAX;}

	for(i=0; i<spots.size(); i++){
		spot_index_map[spot_index(spot_x[i], spot_y[i], max_spot_y)] = i;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// Part 2: pre-compute distance matrices
	n = 2*max_radius*max_radius;

	distance = new double[n];
	memset(distance, 0, n * sizeof(double));

	create_distance(distance, max_radius, platform);

	// same spot has weight zero, if not considering
	if(!flag_samespot) distance[0] = 0;

	////////////////////////////////////////////////////////////////////////////////////////
	// Part 3: compute pairwise moran-I

	n = gene_names.size();
	ncol = spots.size();

	ofstream fout(output.c_str());

	total = n * (n+1) / 2;
	step = 100000;

	cout << total << " total computations" << endl;
	cout << step << " step size in time estimations" << endl;

	nstep = total/step;
	cout << nstep << " total steps" << endl;

	clock_t start_time = clock();
	
	if(allGenes_1stGene)
	{
		for(counter=i=0;i<n;i++)
		{
			for(j=0;j<=i;j++, counter++)
			{
				if(counter % step == 0)
				{
					cout << "Progress: " << 100.0 * counter / total << "%";
	
					if(counter > 0){
						// time elapsed until now
						val = (clock() - start_time) / (double)CLOCKS_PER_SEC;
	
						// time left from now
						val = val * (total - counter)/counter;
						cout << "\tEstimated time left: " << val/3600 << " hours";
					}
	
					cout << endl;
				};
	
				fout << moran_I(
					data + i*ncol, data + j*ncol, ncol,
					distance, max_radius,
					spot_x, spot_y, max_spot_x, max_spot_y,
					spot_index_map
					) << (j==i? '\n': '\t');
			}
		}
		
	}else{
		for(counter=i=0;i<n;i++)
		{
			for(j=0;j<=0;j++, counter++)
			{
				if(counter % step == 0)
				{
					cout << "Progress: " << 100.0 * counter / total << "%";
	
					if(counter > 0){
						// time elapsed until now
						val = (clock() - start_time) / (double)CLOCKS_PER_SEC;
	
						// time left from now
						val = val * (total - counter)/counter;
						cout << "\tEstimated time left: " << val/3600 << " hours";
					}
	
					cout << endl;
				};
	
				fout << moran_I(
					data + i*ncol, data + j*ncol, ncol,
					distance, max_radius,
					spot_x, spot_y, max_spot_x, max_spot_y,
					spot_index_map
					) << (j==i? '\n': '\t');
			}
		}
		
	}

	fout.close();

	delete [] data;
	delete [] distance;
	delete [] spot_x;
	delete [] spot_y;
	delete [] spot_index_map;

	return 0;
}
