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
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"

#include "itkCastImageFilter.h"

#include "itkLabelOverlayImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"

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

using LargeIntegerPixelType = uint16_t;
using LargeIntegerImageType = itk::Image<LargeIntegerPixelType, Dimension>;

using RealPixelType = float;
using RealImageType = itk::Image<RealPixelType, Dimension>;

using RGBPixelType = itk::RGBPixel<uint8_t>;
using RGBImageType = itk::Image<RGBPixelType, Dimension>;

using NameGeneratorType = itk::NumericSeriesFileNames;
using ReaderType = itk::ImageFileReader<IntegerImageType>;
using WriterType = itk::ImageFileWriter<IntegerImageType>;
using LargeIntegerWriterType = itk::ImageFileWriter<LargeIntegerImageType>;
using RGBWriterType = itk::ImageFileWriter<RGBImageType>;

using CastToRealType =
  itk::CastImageFilter<IntegerImageType, RealImageType>;
using CastToIntegerType =
  itk::CastImageFilter<RealImageType, IntegerImageType>;

using DistMapType =
  itk::DanielssonDistanceMapImageFilter<RealImageType, RealImageType>;
using InvertIntensityType =
  itk::InvertIntensityImageFilter<RealImageType, RealImageType>;
using WatershedType =
  itk::MorphologicalWatershedImageFilter<RealImageType,
                                         LargeIntegerImageType>;
using WatershedMaskType = itk::MaskImageFilter<LargeIntegerImageType,
                                               RealImageType,
                                               LargeIntegerImageType>;

using RescalerType =
  itk::RescaleIntensityImageFilter<RealImageType, RealImageType>;
using LabelOverlayType = itk::LabelOverlayImageFilter<RealImageType,
                                                      LargeIntegerImageType,
                                                      RGBImageType>;

/******************************************************************************/
/* Reader and writer                                                          */
/******************************************************************************/
// image reader
static void
read_frame(ReaderType::Pointer &reader,
           const string &file,
           IntegerImageType::Pointer &in) {

  reader->SetFileName(file);
  try {
    reader->Update();
    in = reader->GetOutput();
  }
  catch (const itk::ExceptionObject &e) {
    cerr << "ITK EXCEPTION:"  << endl
         << e << endl;
  }
}

static void
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


static void
write_frame_large_int(LargeIntegerWriterType::Pointer &writer,
                      const string &file,
                      const LargeIntegerImageType::Pointer &out) {

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

static void
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
/* Watershed segmentation                                                     */
/******************************************************************************/
static void
label_cells(const RealImageType::Pointer &in,
            DistMapType::Pointer &dist_map,
            InvertIntensityType::Pointer &in_inverter,
            InvertIntensityType::Pointer &dist_inverter,
            WatershedType::Pointer &watershed,
            WatershedMaskType::Pointer &watershed_mask,
            LargeIntegerImageType::Pointer &out,
            const size_t VERBOSE = 0,
            const string out_file_prefix = "") {

  // compute the distance map
  in_inverter->SetInput(in);
  dist_map->SetInput(in_inverter->GetOutput());
  dist_inverter->SetInput(dist_map->GetOutput());
  dist_inverter->Update();

  // watershed the distance map
  watershed->SetInput(dist_inverter->GetOutput());

  // mask out the background regions of watershed image
  watershed_mask->SetInput(watershed->GetOutput());
  watershed_mask->SetMaskImage(in);
  watershed_mask->Update();

  out = watershed_mask->GetOutput();

  // write and distance and label map on high verbosity
  if (VERBOSE >= 2) {
    constexpr size_t PixelMax =
      std::numeric_limits<IntegerPixelType>::max();
    constexpr size_t PixelMin =
      std::numeric_limits<IntegerPixelType>::min();

    RescalerType::Pointer debug_rescaler = RescalerType::New();
    CastToIntegerType::Pointer debug_to_int =
      CastToIntegerType::New();
    WriterType::Pointer debug_writer = WriterType::New();

    // write distance map
    debug_rescaler->SetOutputMinimum(PixelMin);
    debug_rescaler->SetOutputMaximum(PixelMax);

    debug_rescaler->SetInput(dist_inverter->GetOutput());
    debug_to_int->SetInput(debug_rescaler->GetOutput());
    write_frame(debug_writer, out_file_prefix + "_dist",
                debug_to_int->GetOutput());

    // write rgb label map
    LabelOverlayType::Pointer debug_overlay = LabelOverlayType::New();
    debug_overlay->SetInput(in);
    debug_overlay->SetLabelImage(out);
    debug_overlay->SetOpacity(0.3);

    RGBWriterType::Pointer debug_rgb_writer = RGBWriterType::New();
    write_frame_rgb(debug_rgb_writer, out_file_prefix + "_rgb_label",
                    debug_overlay->GetOutput());
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
    if (argc < 7) {
      cerr << argv[0]
           << " <in_dir> <out_dir>"
           << " <start_offset> <num_frames> <watershed_level>"
           << " <verbose>" << endl;
      return EXIT_FAILURE;
    }

    const string in_dir = argv[1];
    const string out_dir = argv[2];

    const size_t start_offset = std::stoi(argv[3]);
    const size_t n_frames = std::stoi(argv[4]);

    const float watershed_level = std::stof(argv[5]);

    const size_t VERBOSE = std::stoi(argv[6]);

    const string in_format = "mask%06d.tif";
    const string outfile_prefix = "label%06d";

    /**************************************************************************/
    /* Initialize filters                                                     */
    /**************************************************************************/
    if (VERBOSE >= 1)
      cerr << "Initializing ITK filters" << endl;

    // constants
    constexpr size_t PixelMax =
      std::numeric_limits<IntegerPixelType>::max();

    // Name genreators, readers, and writers
    NameGeneratorType::Pointer name_gen = NameGeneratorType::New();
    ReaderType::Pointer reader = ReaderType::New();
    LargeIntegerWriterType::Pointer writer = LargeIntegerWriterType::New();

    // casters
    CastToRealType::Pointer to_real = CastToRealType::New();

    // watershed related
    DistMapType::Pointer ws_dist_map = DistMapType::New();
    InvertIntensityType::Pointer ws_in_inverter = InvertIntensityType::New();
    InvertIntensityType::Pointer ws_dist_inverter = InvertIntensityType::New();
    ws_in_inverter->SetMaximum(PixelMax);
    ws_dist_inverter->SetMaximum(PixelMax);

    WatershedType::Pointer ws_watershed = WatershedType::New();
    ws_watershed->SetLevel(watershed_level);
    ws_watershed->SetFullyConnected(true);
    ws_watershed->SetMarkWatershedLine(false);

    WatershedMaskType::Pointer ws_mask = WatershedMaskType::New();

    /**************************************************************************/
    /* Generate file names                                                    */
    /**************************************************************************/
    // input file names
    name_gen->SetSeriesFormat(in_format);
    name_gen->SetStartIndex(start_offset);
    name_gen->SetEndIndex(start_offset + n_frames);
    name_gen->SetIncrementIndex(1);

    vector<string> in_frame_files;
    in_frame_files = name_gen->GetFileNames();

    // output file names
    name_gen->SetSeriesFormat(outfile_prefix);
    name_gen->SetStartIndex(start_offset);
    name_gen->SetEndIndex(start_offset + n_frames);
    name_gen->SetIncrementIndex(1);

    vector<string> out_frame_files;
    out_frame_files = name_gen->GetFileNames();

    /**************************************************************************/
    /* Process frames                                                         */
    /**************************************************************************/
    if (VERBOSE >= 1)
      cerr << "Processing frames" << endl;


    for (size_t i = 0; i < n_frames; ++i) {
      if (VERBOSE >= 1) {
        if (!(i % 100)) {
          cerr << "\tProcessed " << i << " frames" << endl;
        }
      }

      // read input frame
      IntegerImageType::Pointer in_frame;
      read_frame(reader, in_dir + in_frame_files[i], in_frame);

      // convert to float
      to_real->SetInput(in_frame);

      // TODO: merge masks

      // perfrom watershed
      LargeIntegerImageType::Pointer label_mask;
      label_cells(to_real->GetOutput(), ws_dist_map, ws_in_inverter,
                  ws_dist_inverter, ws_watershed, ws_mask, label_mask,
                  VERBOSE, out_dir + out_frame_files[i]);

      // write label mask
      write_frame_large_int(writer, out_dir + out_frame_files[i], label_mask);
    }
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
