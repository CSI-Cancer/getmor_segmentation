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
/* Misc.                                                                      */
/******************************************************************************/
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
    if (argc < 4) {
      cerr << argv[0]
           << " <in_features_file> <out_features_file>"
           << " <col(1-based),min,max>*" << endl;
      return EXIT_FAILURE;
    }

    const string in_file = argv[1];
    const string out_file = argv[2];

    vector<size_t> feature_col;
    vector<float> feature_min;
    vector<float> feature_max;


    for (int i = 3; i < argc; ++i) {
      vector<string> tokens;
      split_string(argv[i], tokens, ',');
      if (tokens.size() != 3) {
        cerr << "ERROR: not in <col,min,max> format" << endl;
        return EXIT_FAILURE;
      }
      feature_col.push_back(std::stoi(tokens[0]));
      feature_min.push_back(std::stof(tokens[1]));
      feature_max.push_back(std::stof(tokens[2]));
    }

    /**************************************************************************/
    /* Process features file                                                  */
    /**************************************************************************/
    // open input file
    std::ifstream in_features(in_file);
    if (!in_features) {
      cerr << "ERROR: Cannot open " + in_file << endl;
      return EXIT_FAILURE;
    }

    // open output file
    std::ofstream out_features(out_file);
    if (!out_features) {
      cerr << "ERROR: Cannot open " + out_file << endl;
      return EXIT_FAILURE;
    }

    string line;
    while (getline(in_features, line)) {
      vector<string> tokens;
      split_string(line, tokens);

      bool discard = false;
      size_t i = 0;
      while (!discard && i < feature_col.size()) {
        float val = std::stof(tokens[feature_col[i] - 1]);
        if ((val < feature_min[i]) || val > feature_max[i]) {
          discard = true;
        }
        ++i;
      }

      if (!discard) {
        out_features << line << endl;
      }

    }

    // close files
    in_features.close();
    out_features.close();

  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
