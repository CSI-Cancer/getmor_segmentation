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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNumericSeriesFileNames.h"


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

using ReaderType = itk::ImageFileReader<IntegerImageType>;
using WriterType = itk::ImageFileWriter<IntegerImageType>;

using NameGeneratorType = itk::NumericSeriesFileNames;

using ConstImageIteratorType =
  itk::ImageRegionConstIterator<IntegerImageType>;
using ImageIteratorType =
  itk::ImageRegionIterator<IntegerImageType>;

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
/* Main                                                                       */
/******************************************************************************/
int
main (int argc, char* argv[]) {
  try {
    /**************************************************************************/
    /* Parse arguments                                                        */
    /**************************************************************************/
    if (argc != 6) {
      cerr << argv[0]
           << " <in_dir> <out_dir> <start_offset> <num_frames>"
           << " <gain>" << endl;
      return EXIT_FAILURE;
    }

    const string in_dir = argv[1];
    const string out_dir = argv[2];
    if (in_dir == out_dir) {
      cerr << "Input and outpur directories cannot be the same." << endl;
      return EXIT_FAILURE;
    }

    const size_t start_offset = std::stoi(argv[3]);
    const size_t n_frames = std::stoi(argv[4]);

    const float gain = std::stof(argv[5]);

    const string frame_format = "Tile%06d.tif";

    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    const size_t PixelMax =
      std::numeric_limits<IntegerPixelType>::max();

    // Readers and writers
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();

    // generate file names
    NameGeneratorType::Pointer name_gen = NameGeneratorType::New();
    name_gen->SetStartIndex(start_offset);
    name_gen->SetEndIndex(start_offset + n_frames - 1);
    name_gen->SetIncrementIndex(1);

    vector<string> frame_files;
    name_gen->SetSeriesFormat(frame_format);
    frame_files = name_gen->GetFileNames();

    // dummy read first image to set output image dimensions
    IntegerImageType::Pointer in_img;
    read_frame(reader, in_dir + frame_files[0], in_img);

    // Allocate output image
    IntegerImageType::Pointer out_img = IntegerImageType::New();

    IntegerImageType::RegionType::SizeType out_size;
    out_size[0] = in_img->GetRequestedRegion().GetSize()[0];
    out_size[1] = in_img->GetRequestedRegion().GetSize()[1];

    IntegerImageType::RegionType::IndexType out_start;
    out_start[0] = 0;
    out_start[1] = 0;

    IntegerImageType::RegionType out_region;
    out_region.SetSize(out_size);
    out_region.SetIndex(out_start);

    out_img->SetOrigin(in_img->GetOrigin());
    out_img->SetSpacing(in_img->GetSpacing());
    out_img->SetRegions(out_region);
    out_img->Allocate();


    /**************************************************************************/
    /* Process frames                                                         */
    /**************************************************************************/
    for (size_t i = 0; i < n_frames; ++i) {

      // read image
      read_frame(reader, in_dir + frame_files[i], in_img);

      // iterate over image and set gain
      ConstImageIteratorType in_it(in_img, in_img->GetRequestedRegion());
      ImageIteratorType out_it(out_img, out_img->GetRequestedRegion());

      in_it.GoToBegin();
      while (!in_it.IsAtEnd()) {

        float out_val = round(in_it.Get() * gain);
        if (out_val > PixelMax)
          out_val = PixelMax;
        out_it.Set(out_val);

        ++out_it;
        ++in_it;
      }
      write_frame(writer, out_dir + frame_files[i], out_img);
    }


  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
