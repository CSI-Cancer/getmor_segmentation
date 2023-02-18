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


/******************************************************************************/
/* Type definitions                                                           */
/******************************************************************************/
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

using NameGeneratorType = itk::NumericSeriesFileNames;

using ReaderType = itk::ImageFileReader<IntegerImageType>;
using WriterType = itk::ImageFileWriter<IntegerImageType>;
using RGBWriterType = itk::ImageFileWriter<RGBImageType>;

/******************************************************************************/
/* Reader and writer                                                          */
/******************************************************************************/
void
read_frame(ReaderType::Pointer &reader,
           const string &file,
           IntegerImageType::Pointer &in) {

  reader->SetFileName(file);
  try {
    reader->Update();
    in = reader->GetOutput();
  }
  catch (const itk::ExceptionObject &e) {
    cerr << "ITK EXCEPTION:" << endl
         << e << endl;
  }
}


void
write_frame(WriterType::Pointer &writer,
            const string &file,
            const IntegerImageType::Pointer &out) {

  writer->SetFileName(file + ".tif");
  writer->SetInput(out);
  try {
    writer->Update();
  }
  catch (const itk::ExceptionObject &e) {
    cerr << "ITK EXCEPTION:" << endl
         << e << endl;
  }
}

/******************************************************************************/
/* Misc.                                                                      */
/******************************************************************************/
// conver a csv string to vector of ints
static void
csv_to_uint(const string &in,
            vector<size_t> &tokens) {

  tokens.clear();
  size_t pos = 0;
  const char delim = ',';
  size_t next_pos = in.find(delim);
  while (next_pos != string::npos) {
    tokens.push_back(std::stoi(in.substr(pos, next_pos - pos)));
    pos = ++next_pos;
    next_pos = in.find(delim, pos);
  }
  tokens.push_back(std::stoi(in.substr(pos)));
}

/******************************************************************************/
/* Main                                                                       */
/******************************************************************************/
int
main (int argc, char* argv[]) {
  try {
    /**************************************************************************/
    /* Parse arguments                                                        */
    /**************************************************************************/
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

    const string ch_to_use_str = argv[4];
    vector<size_t> ch_to_use;
    csv_to_uint(ch_to_use_str, ch_to_use);
    if (ch_to_use.size() != 4) {
      cerr << "ERROR: ch_use_vector must of length 4" << endl;
      return EXIT_FAILURE;
    }
    

    const size_t red_start = std::stoi(argv[5]);
    const size_t green_start = std::stoi(argv[6]);
    const size_t blue_start = std::stoi(argv[7]);
    const size_t white_start = std::stoi(argv[8]);
    const size_t frame_len = std::stoi(argv[9]); 

    const string in_format = "Tile%06d.tif";

    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    RGBWriterType::Pointer rgb_writer = RGBWriterType::New();
    
    NameGeneratorType::Pointer name_gen = NameGeneratorType::New();

    // calculate number of channels
    size_t n_channels = 0;
    for (auto it = ch_to_use.begin(); it != ch_to_use.end(); ++it) {
      if (*it != 0) ++n_channels;
    } 
    cout << "number of ch: " << n_channels << endl;

    // Allocate output tile image
    IntegerImageType::Pointer tile_img = IntegerImageType::New();
    
    IntegerImageType::RegionType::SizeType tile_size;
    tile_size[0] = frame_len * n_channels;
    tile_size[1] = frame_len; 
     
    IntegerImageType::RegionType::IndexType tile_start;
    tile_start[0] = 0;
    tile_start[1] = 0;

    IntegerImageType::RegionType tile_region;
    tile_region.SetSize(tile_size);
    tile_region.SetIndex(tile_start);

    IntegerImageType::PointType tile_origin;
    tile_origin[0] = 0;
    tile_origin[1] = 0;

    IntegerImageType::SpacingType tile_spacing;
    tile_spacing[0] = 1;
    tile_spacing[1] = 1;

    tile_img->SetOrigin(tile_origin);
    tile_img->SetSpacing(tile_spacing);
    tile_img->SetRegions(tile_region);
    tile_img->Allocate();

    // Allocate output RGB image
    RGBImageType::Pointer rgb_img = RGBImageType::New();
    
    RGBImageType::RegionType::SizeType rgb_size;
    rgb_size[0] = frame_len;
    rgb_size[1] = frame_len;
    
    RGBImageType::RegionType::IndexType rgb_start;
    rgb_start[0] = 0;
    rgb_start[1] = 0;

    RGBImageType::RegionType rgb_region;
    rgb_region.SetSize(rgb_size);
    rgb_region.SetIndex(rgb_start);

    RGBImageType::PointType rgb_origin;
    rgb_origin[0] = 0;
    rgb_origin[1] = 0;

    RGBImageType::SpacingType rgb_spacing;
    rgb_spacing[0] = 1;
    rgb_spacing[1] = 1;

    rgb_img->SetOrigin(rgb_origin);
    rgb_img->SetSpacing(rgb_spacing);
    rgb_img->SetRegions(rgb_region);
    rgb_img->Allocate();
      

 
 

    /**************************************************************************/
    /* Process frames                                                         */
    /**************************************************************************/



  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
