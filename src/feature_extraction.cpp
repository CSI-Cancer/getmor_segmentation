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

using NameGeneratorType = itk::NumericSeriesFileNames;
using ReaderType = itk::ImageFileReader<IntegerImageType>;
using LargeIntegerReaderType = itk::ImageFileReader<LargeIntegerImageType>;

/******************************************************************************/
/* Readers                                                                    */
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
           << " <in_dir> <out_dir> <sample_name>"
           << " <start_offset> <num_frames> <channel_start>" << endl;
      return EXIT_FAILURE;
    }

    const string in_dir = argv[1];
    const string out_dir = argv[2];
    const string sample_name = argv[3];
    const size_t start_offset = std::stoi(argv[4]);
    const size_t n_frames = std::stoi(argv[5]);

    const string channel_start_str = argv[6];
    vector<size_t> channel_start;
    csv_to_uint(channel_start_str, channel_start);
    const size_t n_channels = channel_start.size();

    const string label_format = "watershed%04d";
    const string frame_format = "Tile%06d";

    /**************************************************************************/
    /* Initialize filters                                                     */
    /**************************************************************************/
    
    // Name genreators, readers, and writers
    NameGeneratorType::Pointer name_gen = NameGeneratorType::New();
    ReaderType::Pointer reader = ReaderType::New();
    LargeIntegerReaderType::Pointer reader_large_int = 
      LargeIntegerReaderType::New(); 
 
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

    }
  
    /**************************************************************************/
    /* Clean-up and exit                                                      */
    /**************************************************************************/

  } 
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } 
  return EXIT_SUCCESS;
}
