#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <time.h>
using namespace std;

// Definition for a 2D pixel point.
// Add additional fields if necessary.
// 
class Coordinate
{  
public:
  Coordinate() :row(0), col(0) {};
  double row, col;
};

// Ground truth nformation for training or testing image
class GTImage
{
public:
  string image_filename;
  Coordinate coords[6];
};

// Code for computing unary potential functions
//  Feel free to modify this if you want, but 
//  it should work as-is.
#include <UnaryPotential.h>

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlay rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// Load ground truth coordinates and image names from file
//
vector<GTImage> load_groundtruth(const string &gt_filename)
{
  // read in training images, one at a time, and put the ground truth coordinates
  //  into a vector.
  ifstream ifs(gt_filename.c_str());
  
  vector<GTImage> training_images;
  while(ifs.good())
    {
      string image_number;
      ifs >> image_number;
      if(!ifs.good())	break;
      
      GTImage gt;
      gt.image_filename = string("imgs/image_") + image_number + ".png";
      for(int ii=0; ii<6; ii++)
	ifs >> gt.coords[ii].col >> gt.coords[ii].row;

      training_images.push_back(gt);
    }

  return training_images;
}


// Draw a given part configuration on an image, using simple
//  colored rectangles.
// Color code is: 
// left eye: bright green
// right eye: dark green
// nose: bright red
// left mouth: bright blue
// right mouth: dark blue
// chin: dark red
void draw_configuration(const char *output_fname, const char *sample_filename, const vector<Coordinate> &configuration)
{
  int g_code[] = {255, 128,   0,   0,   0, 0};
  int r_code[] = {0, 0,     255,   0,   0, 128};
  int b_code[] = {0, 0,       0, 255, 128, 0};

  SDoublePlane img = SImageIO::read_png_file(sample_filename);
  SDoublePlane r = img, g = img, b = img;
  for(int ii=0; ii<configuration.size(); ii++)
    {
      overlay_rectangle(g, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, g_code[ii], 3);
      overlay_rectangle(r, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, r_code[ii], 3);
      overlay_rectangle(b, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, b_code[ii], 3);
    }

  SImageIO::write_png_file(output_fname, r, g, b);
}


int main(int argc, char *argv[])
{
  try {
    if(!(argc == 3)) {
      cerr << "usage: " << argv[0] << " training_gt_file testing_gt_file" << endl;
      return 1;
    }

    string train_gt_filename(argv[1]), test_gt_filename(argv[2]);

    // This is a parameter used to compute the unary potentials. It specifies
    //  how big the part appearance models should be.
    const int template_size=15;

    // Load in the training data
    vector<GTImage> training_data = load_groundtruth(train_gt_filename);
	// Load the test data
	vector<GTImage> test_data = load_groundtruth(test_gt_filename);

    // Now set up the unary potential functions
    UnaryPotential unary_potentials(training_data, template_size);
    cerr << "Done learning unary potentials" << endl << endl;

	//open an output file
	ofstream ofs;
	ofs.open("output/output.txt");
	
	float total_accuracy =0;
	float total_accuracy_part[6] = {};
	int total_correct = 0;
	int total_correct_part[6] = {};
	
	for (int sample = 0;sample<test_data.size();sample = sample+1)
	{
	string sample_filename = test_data[sample].image_filename;
	
	cout << "the sample_filename is " <<sample_filename<<endl;
	
    // Now run markov network inference on one of the sample images
    //  
	
    //string sample_filename = "imgs/image_0001.png";
    cerr << "Running naive inference on " << sample_filename << "..." << endl;

    // Compute unary potentials for this image
    vector<SDoublePlane> my_potentials = unary_potentials.compute_unary_logpotential(sample_filename);
	

	//learn the parameters mu and sigma

	Coordinate mean[6][6]; //lower triangle matrix
	Coordinate variance[6][6]; //lower triangle matrix
	
	
	const int mt=training_data.size();
	const int mtest = test_data.size();
	
	// the dimension of the picture for each L_i
	const int m = my_potentials[0].rows();
	const int n = my_potentials[0].cols();
		
    // Naive inference just involves choosing best coordinate in each individual
    //  unary potential function.

	const int part_count = 6;
    vector<Coordinate> best_configuration(part_count);
    for(int ii=0; ii<part_count; ii++)
      {
	     float max_val=-1e100;
	for(int i=0; i<my_potentials[ii].rows(); i++)
	  for(int j=0; j<my_potentials[ii].cols(); j++)
	    if(my_potentials[ii][i][j] > max_val)
	      {
		best_configuration[ii].row = i;
		best_configuration[ii].col = j;
		max_val = my_potentials[ii][i][j];
	      }
      }

    // Visualize result by drawing detected parts on the image, and 
    //  outputing to a new image file.
    string output_filename = string("Naive") + sample_filename;
    draw_configuration( output_filename.c_str(), sample_filename.c_str(), best_configuration);
	
	//calculate the accuracy for graph (b)	

	int correct = 0; // case b, case d , naive inference
	int correct_part[6] = {};	
	float accuracy = 0;
	float accuracy_part[6] ={} ; // save the accuracy for each part for each method

		//count for naive 
		for(int k=0 ;k < part_count; k++)
		{
			float d = pow( best_configuration[k].row - test_data[sample].coords[k].row, 2);
				  d = d + pow( best_configuration[k].col - test_data[sample].coords[k].col, 2);
				  d = sqrt(d);
				  
				  if(d <= 20)
				  {
					  correct= correct + 1;//MN (b)		
					  correct_part[k] = correct_part[k]+1		  ;
				  }		  
		}
		
		//calculate the accuracy for b
		 accuracy = (float)correct/ (float)part_count;
		 cout<<"correct answers is "<<correct<<endl;
		 cout<<"accuracy is "<<accuracy<<endl;
					
		//output the accuracy for each sample
		ofs <<"performance for each sample:"<<endl;
		ofs <<"sample_number   Naive_inference "<<endl;
		ofs << sample<<"  "<<accuracy <<endl;
		
		total_correct = total_correct+correct;
		for(int k=0;k<part_count;k++)
		{
		total_correct_part[k] = total_correct_part[k]+correct_part[k];
		}
		
	}//end the loof for sample
	
	
	int mtest = test_data.size();
	//compute the overall accuracy	
	total_accuracy = (float)total_correct / (float)(mtest*6.0);
	
	for(int k=0;k<6;k++)
	{
		total_accuracy_part[k] = total_correct_part[k] /(float)mtest;
	}
		
		//output the overall accuracy
	ofs <<"-------------------------------"<<endl;
	ofs <<"performance for all the sample:"<<endl;
	ofs <<"    Naive_inference"<<endl;
	ofs<<total_accuracy <<"  "<<endl;
	
		ofs << "accuracy for each part using Naive inference"<<endl;
		for(int k =0 ; k< 6; k++)
		{
		ofs << total_accuracy_part[k]<<"  ";	
		}
		ofs<<endl;
		
		
  }catch(std::string &str) {
    cerr << str << endl;
  }

  return 0;
}
