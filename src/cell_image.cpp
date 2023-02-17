/*
* Copyright (C) 2023 Rishvanth Prabakar
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageRegionIteratorWithIndex.h"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

constexpr size_t Dimension = 2;

using IntegerPixelType = uint8_t;
using IntegerImageType = itk::Image<IntegerPixelType, Dimension>;
using RGBPixelType = itk::RGBPixel<IntegerPixelType>;
using RGBImageType = itk::Image<RGBPixelType, Dimension>;


int
main (int argc, char* argv[]) {
  try {
    if (argc != 10) {
      cerr << argv[0]
           << " <frame_dir> <out_dir> <cell_list_file>" << endl
           << "\t<ch_use_vector> <red_start> <green_start>"
           << " <blue_start> <white_start> <out_frame_size>" << endl;
      return EXIT_FAILURE;
    }

    const string frame_dir = argv[1];
    const string out_dir = argv[2];
    const string cell_list_file = argv[3];

    const string ch_to_use = argv[4];
    const size_t red_start = std::stoi(argv[5]);
    const size_t green_start = std::stoi(argv[6]);
    const size_t blue_start = std::stoi(argv[7]);
    const size_t white_start = std::stoi(argv[8]);
    const size_t frame_len = std::stoi(argv[9]); 

    const string in_format = "Tile%06d.tif";


  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
