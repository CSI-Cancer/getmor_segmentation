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

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;


int
main (int argc, char* argv[]) {
  try {
    if (argc != 9) {
      cerr << argv[0]
           << " <frame_dir> <out_dir> <feature_vector>" << endl
           << "\t<ch_use_vector> <red_start> <green_start>"
           << " <blue_start> <white_start> <out_frame_size>" << endl;
      return EXIT_FAILURE;
    }
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
