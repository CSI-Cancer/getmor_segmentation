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

using IntegerImageConstIteratorType =
  itk::ImageRegionConstIteratorWithIndex<IntegerImageType>;
using IntegerImageIteratorType =
  itk::ImageRegionIteratorWithIndex<IntegerImageType>;
using RGBImageIteratorType =
  itk::ImageRegionIteratorWithIndex<RGBImageType>;

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

void
write_frame_rgb(RGBWriterType::Pointer &writer,
                const string &file,
                const RGBImageType::Pointer &out) {

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

static void
split_string(const string &in,
             vector<string> &tokens,
             const char delim = '\t') {

  tokens.clear();
  size_t pos = 0;
  size_t next_pos = in.find(delim);
  while (next_pos != string::npos) {
    tokens.push_back(in.substr(pos, next_pos - pos));
    pos = ++next_pos;
    next_pos = in.find(delim, pos);
  }
  tokens.push_back(in.substr(pos));
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
    if (argc != 7) {
      cerr << argv[0]
           << " <frame_dir> <out_dir> <cell_list_file>" << endl
           << "\t<ch_use_vector> <ch_start_offset> <out_frame_size>" << endl;
      return EXIT_FAILURE;
    }

    const string frame_dir = argv[1];
    const string out_dir = argv[2];
    const string cell_list_file = argv[3];

    const string ch_to_use_str = argv[4];
    vector<size_t> ch_to_use;
    csv_to_uint(ch_to_use_str, ch_to_use);

    const string ch_start_str = argv[5];
    vector<size_t> ch_start;
    csv_to_uint(ch_start_str, ch_start);
    if ((ch_to_use.size()) != 4 || (ch_start.size() != 4)) {
      cerr << "ERROR: ch_use_vector and ch_start_offset must of length 4"
           << endl;
      return EXIT_FAILURE;
    }

    const size_t frame_len = std::stoi(argv[6]);

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
    std::ifstream cell_list(cell_list_file);
    if (!cell_list) {
      cerr << "Cannot open " + cell_list_file << endl;
      return EXIT_FAILURE;
    }

    string line;
    while (getline(cell_list, line)) {
      // parse cell info
      vector<string> tokens;
      split_string(line, tokens);
      size_t frame_id = std::stoi(tokens[0]);
      size_t cell_id = std::stoi(tokens[1]);
      float cell_x = std::round(std::stof(tokens[2]));
      float cell_y = std::round(std::stof(tokens[3]));

      // zero out the output images
      tile_img->FillBuffer(128);
      rgb_img->FillBuffer(RGBPixelType{});

      size_t tile_count = 0;

      for (size_t i = 0; i < 4; ++i) {
        if (ch_to_use[i]) {
          // read input image
          name_gen->SetSeriesFormat(in_format);
          name_gen->SetStartIndex(frame_id + ch_start[i] - 1);
          name_gen->SetEndIndex(frame_id + ch_start[i] - 1);

          IntegerImageType::Pointer in_frame;
          read_frame(reader, frame_dir + name_gen->GetFileNames()[0],
                     in_frame);


          // define region in inpput image
          IntegerImageType::RegionType crop_region;
          IntegerImageType::RegionType::IndexType crop_index;
          IntegerImageType::RegionType::SizeType crop_size;

          // set the crop start index
          if ((cell_x - (frame_len / 2)) > 0) {
            crop_index[0] = cell_x - (frame_len / 2);
          }
          else {
            crop_index[0] = 0;
          }
          if ((cell_y - (frame_len / 2) > 0)) {
            crop_index[1] = cell_y - (frame_len / 2);
          }
          else {
            crop_index[1] = 0;
          }

          // set the crop size
          size_t in_frame_xlen =
            in_frame->GetRequestedRegion().GetSize()[0];
          if (crop_index[0] + frame_len <= in_frame_xlen) {
            crop_size[0] = frame_len;
          }
          else {
            crop_size[0] = in_frame_xlen - crop_index[0];
          }
          size_t in_frame_ylen =
            in_frame->GetRequestedRegion().GetSize()[1];
          if (crop_index[1] + frame_len <= in_frame_ylen) {
            crop_size[1] = frame_len;
          }
          else {
            crop_size[1] = in_frame_ylen - crop_index[1];
          }

          crop_region.SetIndex(crop_index);
          crop_region.SetSize(crop_size);

          // create parts of the output images
          // Iterators for input and output images
          IntegerImageConstIteratorType in_it(in_frame, crop_region);
          IntegerImageIteratorType tile_it(tile_img,
                                           tile_img->GetRequestedRegion());
          RGBImageIteratorType rgb_it(rgb_img,
                                      rgb_img->GetRequestedRegion());

          // zero offset to properly index output
          in_it.GoToBegin();
          IntegerImageType::OffsetType zero_offset;
          zero_offset[0] = in_it.GetIndex()[0];
          zero_offset[1] = in_it.GetIndex()[1];

          // tile image offset
          IntegerImageType::OffsetType tile_offset;
          tile_offset[0] = frame_len * tile_count;
          tile_offset[1] = 0;


          while (!in_it.IsAtEnd()) {
            // tile image
            tile_it.SetIndex(in_it.GetIndex() - zero_offset + tile_offset);
            tile_it.Set(in_it.Get());

            // rgb image
            rgb_it.SetIndex(in_it.GetIndex() - zero_offset);
            RGBPixelType rgb_pixel;
            if (i == 0) {
              rgb_pixel.SetRed(in_it.Get());
            }
            else if (i == 1) {
              rgb_pixel.SetGreen(in_it.Get());
            }
            else if (i == 2) {
              rgb_pixel.SetBlue(in_it.Get());
            }
            else if (i == 3) {
              rgb_pixel.SetRed(in_it.Get());
              rgb_pixel.SetGreen(in_it.Get());
              rgb_pixel.SetBlue(in_it.Get());
            }

            rgb_it.Set(rgb_pixel + rgb_it.Get());

            ++in_it;
          }
          ++tile_count;

        }

      }
      // write images to file
      write_frame(writer, out_dir + std::to_string(frame_id) + "_" +
                  std::to_string(cell_id) + "_tile", tile_img);

      write_frame_rgb(rgb_writer, out_dir + std::to_string(frame_id) +
                      "_" + std::to_string(cell_id) + "_rgb", rgb_img);

    }

    cell_list.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
