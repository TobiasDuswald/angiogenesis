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

/// @brief  A struct that contains the start and end position of a vessel
/// segment, as well as the radius of the vessel segment.
struct VesselSegment {
  Double3 start_position;
  Double3 end_position;
  double radius;
};

/// @brief Abstract base class for parsing data from a file.
class DataParser {
 public:
  virtual ~DataParser() = default;

  /// Pure virtual function that parses the data from the file.
  virtual void ParseData(const std::string& filename) = 0;

  /// Accessor function for the data.
  VesselSegment& operator[](size_t index) { return data[index]; }

  /// Plot a histogram of the radii and the lengths of the vessel segments.
  void PlotHistograms() const;

  // Member variables
  std::vector<VesselSegment> data;
};

/// @brief A class that parses data from a VTP file.
/// The parser expects data in a VTP file format. Specifically, the data is
/// expected to be in the following format:
/// <?xml version="1.0"?>
/// <VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
///   <PolyData>
///     <Piece NumberOfLines="X" NumberOfPoints="2X">
///       <PointData Scalars="pressure [mmHg]">
///         <DataArray type="Float32" Name="pressure [mmHg]"
///         NumberOfComponents="1" format="ascii">
///           p1 p2 p3 ... p2X
///         </DataArray>
///       </PointData>
///       <CellData Scalars="R">
///         <DataArray type="Float32" Name="R" NumberOfComponents="1"
///         format="ascii">
///           r1 r2 r3 ... rX
///         </DataArray>
///         <DataArray type="Float32" Name="G" NumberOfComponents="1"
///         format="ascii">
///           g1 g2 g3 ... gX
///         </DataArray>
///         <DataArray type="Float32" Name="mu" NumberOfComponents="1"
///         format="ascii">
///           mu1 mu2 mu3 ... muX
///         </DataArray>
///       </CellData>
///       <Points>
///         <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3"
///         format="ascii">
///           s1_1 s1_2 s1_3 e1_1 e1_2 e1_3 s2_1 s2_2 s2_3 e2_1 e2_2 e2_3 ...
///         </DataArray>
///       </Points>
///       <Lines>
///         <DataArray type="Int32" Name="connectivity" NumberOfComponents="1"
///         format="ascii">
///           0 1 2 3 4 5 ... 2X-2 2X-1
///         </DataArray>
///         <DataArray type="Int32" Name="offsets" NumberOfComponents="1"
///         format="ascii">
///           2 4 6 ... 2X-2 2X
///         </DataArray>
///       </Lines>
///     </Piece>
///   </PolyData>
/// </VTKFile>
class DataParserVTP : public DataParser {
 public:
  DataParserVTP() = default;
  ~DataParserVTP() = default;

  /// Parses the data from the VTP file with ROOT's TXMLEngine. The data is
  /// stored in the member variables of the class. After this step, the vectors
  /// pressure_, radii_, g_, mu_, points_, connectivity_, and offsets_ contain
  /// the raw data from the VTP file. No postprocessing is done at this point.
  void ParseData(const std::string& filename) override;

  /// This function defines the origins of each (unconnected) vessels Segment.
  /// This needs to be called before the data is postprocessed because the
  /// problem of mapping from the VTP file to the simulation is not well
  /// posed and requires the user to define the origins of each vessel segment.
  /// Otherwise, the data will be postprocessed incorrectly.
  void SetStartingLines(const std::vector<int>& starting_lines) {
    starting_lines_ = starting_lines;
  }

  /// This function postprocesses the data (points, connectivity, radii) that
  /// were previously parsed from the VTP file. The postprocessing step changes
  /// the affected vectors. Note that other data (pressure, offsets, g, mu) are
  /// not affected by this step, but will not match the other 3 after this step.
  /// This function restructures the data such that:
  /// 1. The points contain only unique points (no duplicates)
  /// 2. The connectivity refers to the unique points
  /// 3. The connectivity is sorted such that it begins with the starting lines
  ///    and then contains batches of the next "layer" of segments. This way we
  ///    ensure that when iterating over the segments, we always have the
  ///    previous segment defined. (Important for the pointer relationship)
  /// 4. The radii are sorted such that they correspond to the sorted
  ///    connectivity.
  void PostProcessData();

  // -------------------------------------------------------------------------
  // Getters (for testing purposes); return a const reference to the
  // corresponding member variable
  // -------------------------------------------------------------------------
  const size_t& GetNumLines() const { return num_lines_; }
  const size_t& GetNumPoints() const { return num_points_; }
  const std::vector<double>& GetPressure() const { return pressure_; }
  const std::vector<double>& GetRadii() const { return radii_; }
  const std::vector<double>& GetG() const { return g_; }
  const std::vector<double>& GetMu() const { return mu_; }
  const std::vector<Double3>& GetPoints() const { return points_; }
  const std::vector<int>& GetConnectivity() const { return connectivity_; }
  const std::vector<int>& GetOffsets() const { return offsets_; }

 private:
  /// Core of the information extraction routine. Called by ParseData. Based on
  /// ROOT's TXMLEngine.
  void RecursivelyParseVTPFile(TXMLEngine& xml, XMLNodePointer_t node,
                               Int_t level);

  std::vector<double> pressure_;     // pressure at each point
  std::vector<double> radii_;        // radius at each segment
  std::vector<double> g_;            // shear modulus at each segment
  std::vector<double> mu_;           // viscosity at each segment
  std::vector<Double3> points_;      // points of each segment
  std::vector<int> connectivity_;    // connectivity of points
  std::vector<int> offsets_;         // offsets of connectivity (useless)
  std::vector<int> starting_lines_;  // starting lines of each vessel
  size_t num_lines_{0};              // number of lines
  size_t num_points_{0};             // number of points
  double x_min_{0};                  // bounding box of the data
  double x_max_{0};                  // bounding box of the data
  double y_min_{0};                  // bounding box of the data
  double y_max_{0};                  // bounding box of the data
  double z_min_{0};                  // bounding box of the data
  double z_max_{0};                  // bounding box of the data
};

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

// Helper function template that takes a string "x y z .." and returns as
// std::vector of type T: std::vector<T> = {x, y, z, ...}
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
// returns the number 123 as type T.
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

/// @brief This function constructs a vector of unique points from a vector of
/// points and a vector of connectivity. The connectivity vector is also
/// modified to refer to the unique points and written to the
/// unique_connectivity.
/// @param points The vector of points
/// @param unique_points The vector of unique points (output)
/// @param connectivity The vector of connectivity
/// @param unique_connectivity The vector of connectivity referring to the
/// unique points (output)
void ConstructUniquePoints(const std::vector<Double3>& points,
                           std::vector<Double3>& unique_points,
                           const std::vector<int>& connectivity,
                           std::vector<int>& unique_connectivity);

/// @brief This function constructs a vector of lines from a vector of
/// connectivity. Each line is represented by a pair of integers, where the
/// first integer is the index of the first point of the line and the second
/// integer is the index of the second point of the line. This function simply
/// iterates over the connectivity vector and constructs the lines.
std::vector<std::pair<int, int>> ConstructLines(
    const std::vector<int>& connectivity);

/// @brief This function verifies if the lines labeled by the starting_lines
/// vector are indeed starting lines. A starting line is a line whose first
/// point is not connected to any other line, but whose second point is
/// connected to at least one other line. This function returns a vector of
/// booleans, where the i-th element is true if the i-th line is a starting
/// line, and false otherwise.
std::vector<bool> VerifyStartingLines(
    const std::vector<int>& starting_lines,
    const std::vector<std::pair<int, int>>& lines);

/// @brief  This function adjusts the starting lines vector to make sure that
/// the lines are indeed starting lines. It does so by swapping the first and
/// second point of each line that is not a starting line. This function calls
/// back to VerifyStartingLines to verify that the lines are indeed starting
/// lines. The function prints the which lines are starting points, which are
/// modified, and which are starting lines after modification. If all lines are
/// starting lines, the function returns true. If not, the function returns
/// false.
/// @param starting_lines
/// @param connectivity
/// @return
bool AdjustStartingLines(const std::vector<int>& starting_lines,
                         std::vector<std::pair<int, int>>& connectivity);

/// @brief This function constructs a tree data structure from a vector of
/// starting lines and a vector of connectivity. By tree, we mean a vector of
/// connections, but the connections are ordered in a way that the first
/// connection is the root of the tree. All following connections have a
/// starting point that is the end point of a previous connection. Additionally,
/// we track the permutation of the connectivity vector that is required to
/// obtain the tree. The permutation can be used to transform other line based
/// data structures (e.g. radii, shear modulus, viscosity) to the tree.
void RestructureToTree(const std::vector<int>& starting_lines,
                       const std::vector<std::pair<int, int>>& connectivity,
                       std::vector<std::pair<int, int>>& tree,
                       std::vector<int>& permutation);

}  // namespace bdm

#endif  // DATA_PARSER_H_
