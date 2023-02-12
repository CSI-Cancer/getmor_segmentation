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

/******************************************************************************/
/* Type definitions                                                           */
/******************************************************************************/
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;


/******************************************************************************/
/* Main                                                                       */
/******************************************************************************/
int
main (int argc, char* argv[]) {
  try {
    /**************************************************************************/
    /* Parse arguments                                                        */
    /**************************************************************************/
    if (argc < 8) {
      cerr << argv[0] 
           << " <in_dir> <out_dir> <sample_name>"
           << " <start_offset> <num_frames>" << endl
           << "\t<num_channels> <channel_start>" << endl;
      return EXIT_FAILURE;
    }

    const string in_dir = argv[1];
    const string out_dir = argv[2];
    const string sample_name = argv[3];
    const size_t start_offset = std::stoi(argv[4]);
    const size_t n_frames = std::stoi(argv[5]);

    const size_t n_channels = std::stoi(argv[6]);
    const string channel_start = argv[7];

    const string in_format = "watershed%04d";
  } 
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } 
  return EXIT_SUCCESS;
}
