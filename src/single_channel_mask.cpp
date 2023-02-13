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

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

#include "itkMaskImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMaskedImageToHistogramFilter.h"
#include "itkBinaryReconstructionByDilationImageFilter.h"

/******************************************************************************/
/* Type definitions                                                           */
/******************************************************************************/

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;

constexpr size_t Dimension = 2;

using IntegerPixelType = uint8_t;
using IntegerImageType = itk::Image<IntegerPixelType, Dimension>;

using RealPixelType = float;
using RealImageType = itk::Image<RealPixelType, Dimension>;

using NameGeneratorType = itk::NumericSeriesFileNames;
using ReaderType = itk::ImageFileReader<IntegerImageType>;
using WriterType = itk::ImageFileWriter<IntegerImageType>;

using CastToRealType =
  itk::CastImageFilter<IntegerImageType, RealImageType>;
using CastToIntegerType =
  itk::CastImageFilter<RealImageType, IntegerImageType>;

using BlurType =
  itk::CurvatureAnisotropicDiffusionImageFilter<RealImageType,
                                                RealImageType>;
using HistogramType =
  itk::Statistics::ImageToHistogramFilter<RealImageType>;
using MaskedHistogramType =
  itk::Statistics::MaskedImageToHistogramFilter<RealImageType,
                                                RealImageType>;
using BinaryThType =
  itk::BinaryThresholdImageFilter<RealImageType, RealImageType>;
using ReconstructType =
  itk::BinaryReconstructionByDilationImageFilter<RealImageType>;

/******************************************************************************/
/* Reader and writer                                                          */
/******************************************************************************/
// image reader
// blocks until a frame is read
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
    cerr << "ITK EXCEPTION:"  << endl
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
/* Blur frame                                                                 */
/******************************************************************************/
void
blur_frame(BlurType::Pointer &smoother,
           const RealImageType::Pointer &in,
           RealImageType::Pointer &out,
           const size_t VERBOSE = 0,
           const string out_file_prefix = "") {

  smoother->SetInput(in);
  smoother->Update();

  out = smoother->GetOutput();

  // write debug output
  if (VERBOSE >= 2) {
    CastToIntegerType::Pointer debug_to_int = CastToIntegerType::New();
    WriterType::Pointer debug_writer = WriterType::New();

    debug_to_int->SetInput(out);
    write_frame(debug_writer, out_file_prefix + "_blur",
                debug_to_int->GetOutput());
  }
}


/******************************************************************************/
/* Double threshold                                                           */
/******************************************************************************/
void
double_threshold(HistogramType::Pointer &low_th_hist,
                 MaskedHistogramType::Pointer &high_th_hist,
                 const float cell_frac,
                 const float low_th_offset,
                 const float high_th_ratio,
                 const float high_th_base_quantile,
                 const RealImageType::Pointer &in,
                 BinaryThType::Pointer &binary_low_th,
                 BinaryThType::Pointer &binary_high_th,
                 ReconstructType::Pointer &reconstructor,
                 RealImageType::Pointer &out,
                 ofstream &stats_file,
                 const size_t VERBOSE = 0,
                 const string out_file_prefix = "") {

  constexpr size_t PixelMax =
    std::numeric_limits<IntegerPixelType>::max();

  // compute low threshold
  low_th_hist->SetInput(in);
  low_th_hist->Update();

  float low_quantile = 1 - cell_frac + low_th_offset;
  if (low_quantile >= 1)
    low_quantile = 0.995;
  else if (low_quantile < 0)
    low_quantile = 0;

  size_t low_th =
    std::round(low_th_hist->GetOutput()->Quantile(0, low_quantile));

  stats_file << "\t" << low_th;


  // binary threshold based on low threshold quantile
  binary_low_th->SetLowerThreshold(low_th);
  binary_low_th->SetInput(in);
  binary_low_th->Update();

  // histogram of masked image
  high_th_hist->SetInput(in);
  high_th_hist->SetMaskImage(binary_low_th->GetOutput());
  high_th_hist->Update();

  // compute upper threshold
  size_t high_th;
  float base_quantile =
    high_th_hist->GetOutput()->Quantile(0, high_th_base_quantile);
  float quantile =
    high_th_hist->GetOutput()->Quantile(0, high_th_base_quantile);
  float quantile_ratio = 0;
  const float increment = 0.001;
  float x = high_th_base_quantile - increment;

  while ((quantile_ratio < high_th_ratio) & (x < 1.0)) {
    x += increment;
    // to avoid rounding error at 1. Ignoring rounding error at
    // other values.
    if ((x > (1.0 - increment/2)) & (x < (1.0 + increment/2)))
      x = 1.0;
    quantile = high_th_hist->GetOutput()->Quantile(0, x);
    quantile_ratio = (quantile / base_quantile) - 1;
  }

  // Nothing to segment
  if (x == 1.0)
    high_th = PixelMax * 2;
  else
    high_th = std::round(quantile);

  stats_file << "\t" << x << "\t" << high_th;


  // binary threshold based on upper threshold
  binary_high_th->SetLowerThreshold(high_th);
  binary_high_th->SetInput(in);
  binary_high_th->Update();

  // reconstruct by dilation the upper thresholded images from the
  // low thresholded image
  reconstructor->SetMaskImage(binary_low_th->GetOutput());
  reconstructor->SetMarkerImage(binary_high_th->GetOutput());
  reconstructor->Update();

  // write output
  out = reconstructor->GetOutput();


  // write debug output
  if (VERBOSE >= 2) {
    CastToIntegerType::Pointer debug_to_int = CastToIntegerType::New();
    WriterType::Pointer debug_writer = WriterType::New();

    debug_to_int->SetInput(binary_low_th->GetOutput());
    write_frame(debug_writer, out_file_prefix + "_low_th",
                debug_to_int->GetOutput());

    debug_to_int->SetInput(binary_high_th->GetOutput());
    write_frame(debug_writer, out_file_prefix + "_high_th",
                debug_to_int->GetOutput());
  }
}

int
main (int argc, char *argv[]) {
  try {
    /**************************************************************************/
    /* Parse arguments                                                        */
    /**************************************************************************/
    if (argc < 12) {
      cerr << argv[0]
           << " <in_dir> <out_dir> <sample_name>"
           << " <start_offset> <num_frames> " << endl
           << "\t<blur_iterations> <cell_frac> <low_th_offset>"
           << " <high_th_base_quantile> <high_th_ratio>" << endl
           << "\t<verbose>" << endl;
      return EXIT_FAILURE;
    }

    const string in_dir = argv[1];
    const string out_dir = argv[2];
    const string sample_name = argv[3];
    const size_t start_offset = std::stoi(argv[4]);
    const size_t n_frames = std::stoi(argv[5]);

    const size_t blur_iterations = std::stoi(argv[6]);

    const float cell_frac_method = std::stof(argv[7]);
    const float low_th_offset = std::stof(argv[8]);
    const float high_th_base_quantile = std::stof(argv[9]);
    const float high_th_ratio = std::stof(argv[10]);

    const size_t VERBOSE = std::stoi(argv[11]);

    const string in_format = "Tile%06d.tif";
    const string outfile_prefix = "mask%06d";

    /**************************************************************************/
    /* Initialize filters                                                     */
    /**************************************************************************/

    if (VERBOSE >= 1)
      cerr << "Initializing ITK filters" << endl;

    // constants
    constexpr size_t PixelMax =
      std::numeric_limits<IntegerPixelType>::max();
    constexpr size_t PixelMin =
      std::numeric_limits<IntegerPixelType>::min();

    // Name genreators, readers, and writers
    NameGeneratorType::Pointer name_gen = NameGeneratorType::New();
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();

    // casters
    CastToRealType::Pointer to_real = CastToRealType::New();
    CastToIntegerType::Pointer to_int = CastToIntegerType::New();

    // blur filter
    const float blur_time_step = 0.125;
    const float blur_conductance = 3;
    BlurType::Pointer blur = BlurType::New();
    blur->SetNumberOfIterations(blur_iterations);
    blur->SetTimeStep(blur_time_step);
    blur->SetConductanceParameter(blur_conductance);

    // cell fraction

    // double threshold
    HistogramType::HistogramSizeType hist_size(1);
    HistogramType::HistogramMeasurementVectorType hist_min(1);
    HistogramType::HistogramMeasurementVectorType hist_max(1);
    hist_size.Fill(PixelMax);
    hist_min.Fill(-0.5);
    hist_max.Fill(PixelMax + 0.5);

    HistogramType::Pointer th_low_hist = HistogramType::New();
    th_low_hist->SetHistogramSize(hist_size);
    th_low_hist->SetHistogramBinMinimum(hist_min);
    th_low_hist->SetHistogramBinMaximum(hist_max);
    th_low_hist->SetMarginalScale(10);

    MaskedHistogramType::Pointer th_high_hist = MaskedHistogramType::New();
    th_high_hist->SetHistogramSize(hist_size);
    th_high_hist->SetHistogramBinMinimum(hist_min);
    th_high_hist->SetHistogramBinMaximum(hist_max);
    th_high_hist->SetMarginalScale(10);
    th_high_hist->SetMaskValue(PixelMax);

    BinaryThType::Pointer th_low_binary = BinaryThType::New();
    th_low_binary->SetOutsideValue(PixelMin);
    th_low_binary->SetInsideValue(PixelMax);

    BinaryThType::Pointer th_high_binary = BinaryThType::New();
    th_high_binary->SetOutsideValue(PixelMin);
    th_high_binary->SetInsideValue(PixelMax);

    ReconstructType::Pointer reconstructor = ReconstructType::New();
    reconstructor->SetBackgroundValue(PixelMin);
    reconstructor->SetForegroundValue(PixelMax);

    // file_handlers
    ofstream stats_file(out_dir + sample_name + "_stats.txt");

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

      stats_file << start_offset + i;

      // read input frame
      IntegerImageType::Pointer in_frame;
      read_frame(reader, in_dir + in_frame_files[i], in_frame);

      // convert to float
      to_real->SetInput(in_frame);

      // blur frame
      RealImageType::Pointer blurred_frame;
      blur_frame(blur, to_real->GetOutput(), blurred_frame,
                 VERBOSE, out_dir + out_frame_files[i]);

      // double threshold
      float cell_frac = 1 - cell_frac_method;
      RealImageType::Pointer double_th_frame;
      double_threshold(th_low_hist, th_high_hist,
                       cell_frac, low_th_offset,
                       high_th_ratio, high_th_base_quantile,
                       blurred_frame, th_low_binary, th_high_binary,
                       reconstructor, double_th_frame,
                       stats_file, VERBOSE, out_dir + out_frame_files[i]);

      // write the frame mask
      to_int->SetInput(double_th_frame);
      write_frame(writer, out_dir + out_frame_files[i],
                  to_int->GetOutput());

      stats_file << endl;
    }


    stats_file.close();
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
