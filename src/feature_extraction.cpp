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
#include "itkNumericSeriesFileNames.h"

#include "itkLabelStatisticsImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

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
using LargeIntegerPixelType = uint16_t;
using LargeIntegerImageType = itk::Image<LargeIntegerPixelType, Dimension>;

using NameGeneratorType = itk::NumericSeriesFileNames;
using ReaderType = itk::ImageFileReader<IntegerImageType>;
using LargeIntegerReaderType = itk::ImageFileReader<LargeIntegerImageType>;

using LabelStatsType = itk::LabelStatisticsImageFilter<IntegerImageType,
                                                       LargeIntegerImageType>;
using ShapeLabelObjectType = itk::ShapeLabelObject<LargeIntegerPixelType,
                                                   Dimension>;
using ShapeLabelMapType = itk::LabelMap<ShapeLabelObjectType>;
using LabelShapeType =
  itk::LabelImageToShapeLabelMapFilter<LargeIntegerImageType,
                                       ShapeLabelMapType>;


/******************************************************************************/
/* Readers                                                                    */
/******************************************************************************/
// image reader
static void
read_img(ReaderType::Pointer &reader,
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

// large int frame frame reader
static void
read_img(LargeIntegerReaderType::Pointer &reader,
         const string &file,
         LargeIntegerImageType::Pointer &in) {

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


/******************************************************************************/
/* Feature extraction                                                         */
/******************************************************************************/
// Cell shape features: depends only on the mask
// For each cell, add the the cell id, cell_x, cell_y, number of pixels,
// elongation, and roundness.
static void
shape_features(const LargeIntegerImageType::Pointer &label_mask,
               LabelShapeType::Pointer &shape_stats,
               vector<vector<float>> &feature_vector) {

  shape_stats->SetInput(label_mask);
  shape_stats->Update();

  ShapeLabelMapType::Pointer shape_label_map = shape_stats->GetOutput();

  size_t n_objects = shape_label_map->GetNumberOfLabelObjects();
  feature_vector.resize(n_objects + 1);

  for (size_t i = 0; i < n_objects; ++i) {
    ShapeLabelObjectType::Pointer label_object =
      shape_label_map->GetNthLabelObject(i);

    feature_vector[i+1].push_back(i+1);
    feature_vector[i+1].push_back(label_object->GetCentroid()[0]);
    feature_vector[i+1].push_back(label_object->GetCentroid()[1]);
    feature_vector[i+1].push_back(label_object->GetNumberOfPixels());
    // feature_vector[i+1].push_back(label_object->GetNumberOfPixelsOnBorder());
    feature_vector[i+1].push_back(label_object->GetElongation());
    feature_vector[i+1].push_back(label_object->GetRoundness());
    // feature_vector[i+1].push_back(label_object->GetFlatness());
  }
}


// Cell intensity features: computed for each channel based on the mask
// For each cell and for each channel, add the the mean intensity,
// standard deviation of intensity, minimum intensity, and maximum
// intensity.
static void
intensity_features(const IntegerImageType::Pointer &in,
                   const LargeIntegerImageType::Pointer &label_mask,
                   LabelStatsType::Pointer &label_stats,
                   vector<vector<float>> &feature_vector) {

  label_stats->SetLabelInput(label_mask);
  label_stats->SetInput(in);
  label_stats->Update();

  // cout << feature_vector.size() << endl;
  // cout << label_stats->GetNumberOfLabels() << endl;
  if (label_stats->GetNumberOfLabels() > 1)
    assert(feature_vector.size() == label_stats->GetNumberOfLabels());

  for (auto it = label_stats->GetValidLabelValues().begin();
       it != label_stats->GetValidLabelValues().end(); ++it) {

    feature_vector[*it].push_back(label_stats->GetMean(*it));
    feature_vector[*it].push_back(label_stats->GetSigma(*it));
    feature_vector[*it].push_back(label_stats->GetMinimum(*it));
    feature_vector[*it].push_back(label_stats->GetMaximum(*it));
  }

}


// format for writing to file
static string
format_feature_vector(const size_t frame_id,
                      const vector<vector<float>> feature_vec) {
  std::ostringstream oss;
  for (auto it = feature_vec.begin() + 1; it != feature_vec.end(); ++it) {
    oss << frame_id << "\t";
    for (auto jt = (*it).begin(); jt != (*it).end(); ++jt) {
      oss << *jt << "\t";
    }
    oss << endl;
  }
  return oss.str();
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

int
main (int argc, char* argv[]) {
  try {
    /**************************************************************************/
    /* Parse arguments                                                        */
    /**************************************************************************/
    if (argc < 7) {
      cerr << argv[0]
           << " <frame_dir> <label_dir> <sample_name>"
           << " <start_offset> <num_frames> <channel_start>" << endl;
      return EXIT_FAILURE;
    }

    const string frame_dir = argv[1];
    const string label_dir = argv[2];
    const string sample_name = argv[3];
    const size_t start_offset = std::stoi(argv[4]);
    const size_t n_frames = std::stoi(argv[5]);

    const string channel_start_str = argv[6];
    vector<size_t> channel_start;
    csv_to_uint(channel_start_str, channel_start);
    const size_t n_channels = channel_start.size();

    const string label_format = "label%06d.tif";
    const string frame_format = "Tile%06d.tif";

    /**************************************************************************/
    /* Initialize filters                                                     */
    /**************************************************************************/

    // Name genreators, readers, and writers
    NameGeneratorType::Pointer name_gen = NameGeneratorType::New();
    ReaderType::Pointer reader = ReaderType::New();
    LargeIntegerReaderType::Pointer reader_large_int =
      LargeIntegerReaderType::New();

    // feature extraction
    LabelStatsType::Pointer feature_intensity_stats = LabelStatsType::New();
    LabelShapeType::Pointer feature_shape_stats = LabelShapeType::New();

    // file handlers
    ofstream features_file(label_dir + sample_name + "_feature_vec.txt");

    /**************************************************************************/
    /* Generate file names                                                    */
    /**************************************************************************/
    // label file names
    name_gen->SetSeriesFormat(label_format);
    name_gen->SetStartIndex(start_offset);
    name_gen->SetEndIndex(start_offset + n_frames - 1);
    name_gen->SetIncrementIndex(1);

    vector<string> label_files;
    label_files = name_gen->GetFileNames();

    // scanned frame file names
    vector<vector<string>> in_frame_files;
    for (size_t i = 0; i < n_channels; ++i) {
      name_gen->SetSeriesFormat(frame_format);
      name_gen->SetStartIndex(channel_start[i]);
      name_gen->SetEndIndex(channel_start[i] + n_frames - 1);
      name_gen->SetIncrementIndex(1);

      in_frame_files.push_back(name_gen->GetFileNames());
    }

    /**************************************************************************/
    /* Process frames                                                         */
    /**************************************************************************/
    for (size_t i = 0; i < n_frames; ++i) {

      // read label
      LargeIntegerImageType::Pointer in_label;
      read_img(reader_large_int, label_dir + label_files[i], in_label);

      // generate shape features
      // feature extraction
      vector<vector<float>> feature_vector;
      shape_features(in_label, feature_shape_stats, feature_vector);

      // generate intensity features for each channel
      for (size_t j = 0; j < n_channels; ++j) {
        IntegerImageType::Pointer in_frame;
        read_img(reader, frame_dir + in_frame_files[j][i], in_frame);
        intensity_features(in_frame, in_label,
                           feature_intensity_stats, feature_vector);

      }

      // write feature vector to file
      features_file << format_feature_vector(start_offset + i, feature_vector);
    }

    /**************************************************************************/
    /* Clean-up and exit                                                      */
    /**************************************************************************/
    features_file.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
