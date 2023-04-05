// -----------------------------------------------------------------------------
//
// Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
//
// -----------------------------------------------------------------------------

#ifndef DATA_PARSER_H_
#define DATA_PARSER_H_

#include <string>
#include "TXMLEngine.h"
#include "biodynamo.h"

namespace bdm {

struct VesselSegment {
  Double3 start_position;
  Double3 end_position;
  double radius;
};

class DataParser {
 public:
  virtual ~DataParser() = default;

  virtual void ParseData(const std::string& filename) = 0;

  VesselSegment& operator[](size_t index) { return data[index]; }

  void PlotHistograms() const;

  // Member variables
  std::vector<VesselSegment> data;
};

class DataParserVTP : public DataParser {
 public:
  DataParserVTP() = default;
  ~DataParserVTP() = default;

  void ParseData(const std::string& filename) override;

  // Getters (for testing purposes); return a const reference to the
  // corresponding member variable
  const size_t& GetNumLines() const { return num_lines_; }
  const size_t& GetNumPoints() const { return num_points_; }
  const std::vector<double>& GetPressure() const { return pressure_; }
  const std::vector<double>& GetRadii() const { return radii_; }
  const std::vector<double>& GetG() const { return g_; }
  const std::vector<double>& GetMu() const { return mu_; }
  const std::vector<Double3>& GetPoints() const { return points_; }
  const std::vector<int>& GetConnectivity() const { return connectivity_; }
  const std::vector<int>& GetOffsets() const { return offsets_; }

  // Set the starting_lines vector
  void SetStartingLines(const std::vector<int>& starting_lines) {
    starting_lines_ = starting_lines;
  }

  void PostProcessData();

 private:
  // ROOT based function to parse VTP file
  void RecursivelyParseVTPFile(TXMLEngine& xml, XMLNodePointer_t node,
                               Int_t level);

  std::vector<double> pressure_;
  std::vector<double> radii_;
  std::vector<double> g_;
  std::vector<double> mu_;
  std::vector<Double3> points_;
  std::vector<int> connectivity_;
  std::vector<int> offsets_;
  // std::vector<size_t> starting_points_;
  std::vector<int> starting_lines_;
  size_t num_lines_{0};
  size_t num_points_{0};
  double x_min_{0};
  double x_max_{0};
  double y_min_{0};
  double y_max_{0};
  double z_min_{0};
  double z_max_{0};
};

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

// Helper function template that takes a string "x y z .." and returns as
// std::vector of type T
template <typename T>
std::vector<T> ParseString(const std::string& str) {
  std::vector<T> result;
  // Split the string into a vector of strings using a space as a delimiter
  std::vector<std::string> tokens;
  std::istringstream iss(str);
  std::string token;
  while (std::getline(iss, token, ' ')) {
    // Remove any leading or trailing whitespace
    token.erase(0, token.find_first_not_of(' '));
    token.erase(token.find_last_not_of(' ') + 1);
    // Remove any newline characters from the string
    token.erase(std::remove(token.begin(), token.end(), '\n'), token.end());
    // Add to the vector of tokens if the token is not empty
    if (!token.empty()) {
      tokens.push_back(token);
    }
  }
  // Convert the tokens to type T and add them to the result vector
  for (auto& token : tokens) {
    result.push_back(std::stod(token));
  }
  return result;
};

// Template helper function that takes a string "some text 123 some text" and
// returns the number 123 as type T
template <typename T>
T ParseStringForNumber(const std::string& str) {
  std::string result;
  // Split the string into a vector of strings using a space as a delimiter
  std::vector<std::string> tokens;
  std::istringstream iss(str);
  std::string token;
  while (std::getline(iss, token, ' ')) {
    tokens.push_back(token);
  }
  // Convert the tokens to type T and add them to the result vector
  for (auto& token : tokens) {
    if (std::isdigit(token[0])) {
      result = token;
    }
  }
  return std::stod(result);
};

void ConstructUniquePoints(const std::vector<Double3>& points,
                           std::vector<Double3>& unique_points,
                           const std::vector<int>& connectivity,
                           std::vector<int>& unique_connectivity);

std::vector<std::pair<int, int>> ConstructLines(
    const std::vector<int>& connectivity);

std::vector<bool> VerifyStartingLines(
    const std::vector<int>& starting_lines,
    const std::vector<std::pair<int, int>>& lines);

bool AdjustStartingLines(const std::vector<int>& starting_lines,
                         std::vector<std::pair<int, int>>& connectivity);

void RestructureToTree(const std::vector<int>& starting_lines,
                       const std::vector<std::pair<int, int>>& connectivity,
                       std::vector<std::pair<int, int>>& tree,
                       std::vector<int>& permutation);

}  // namespace bdm

#endif  // DATA_PARSER_H_
