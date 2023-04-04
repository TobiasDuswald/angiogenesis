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

#include "data_parser.h"
#include <experimental/filesystem>
#include "analysis.h"
#include "util/vector_operations.h"

namespace bdm {

// Inline helper function that checks if the second string is contained in the
// first string
inline bool CheckIfContained(const std::string& str,
                             const std::string& sub_str) {
  return str.find(sub_str) != std::string::npos;
}

void DataParser::PlotHistograms() const {
  // std::string path = Simulation::GetActive()->GetOutputDir();
  std::string path = std::experimental::filesystem::current_path();
  std::string filename_radius = "radius_histogram";
  std::string filename_length = "length_histogram";

  // Create a std::vector of the radii
  std::vector<double> radii;
  std::vector<double> lengths;
  radii.reserve(data.size());
  lengths.reserve(data.size());
  for (auto& segment : data) {
    radii.push_back(segment.radius);
    lengths.push_back((segment.start_position - segment.end_position).Norm());
  }

  // Create histograms
  PlotAndSaveHistogram(radii, filename_radius, path);
  PlotAndSaveHistogram(lengths, filename_length, path);
}

// This function parses VTK file of type polydata. The structure of the file is
// as follows:
//
// <?xml version="1.0"?>
// <VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
//   <PolyData>
//     <Piece NumberOfLines="X" NumberOfPoints="2X">
//       <PointData Scalars="pressure [mmHg]">
//         <DataArray type="Float32" Name="pressure [mmHg]"
//         NumberOfComponents="1" format="ascii">
//           ...
//         </DataArray>
//       </PointData>
//       <CellData Scalars="R">
//         <DataArray type="Float32" Name="R" NumberOfComponents="1"
//         format="ascii">
//           ...
//         </DataArray>
//         <DataArray type="Float32" Name="G" NumberOfComponents="1"
//         format="ascii">
//           ...
//         </DataArray>
//         <DataArray type="Float32" Name="mu" NumberOfComponents="1"
//         format="ascii">
//           ...
//         </DataArray>
//       </CellData>
//       <Points>
//         <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3"
//         format="ascii">
//           ...
//         </DataArray>
//       </Points>
//       <Lines>
//         <DataArray type="Int32" Name="connectivity" NumberOfComponents="1"
//         format="ascii">
//           ...
//         </DataArray>
//         <DataArray type="Int32" Name="offsets" NumberOfComponents="1"
//         format="ascii">
//           ...
//         </DataArray>
//       </Lines>
//     </Piece>
//   </PolyData>
// </VTKFile>
// We are interested in the following data:
// - radius of the vessel segment
// - start and end position of the vessel segment
void DataParserVTP::ParseData(const std::string& filename) {
  // First create engine
  TXMLEngine xml;

  // Parse file (limited syntax for this xml engine)
  XMLDocPointer_t xmldoc = xml.ParseFile(filename.c_str());
  if (!xmldoc) {
    Log::Fatal("DataParserVTP", "Could not parse file ", filename);
    return;
  }

  // Get the main node of the file
  XMLNodePointer_t mainnode = xml.DocGetRootElement(xmldoc);
  RecursivelyParseVTPFile(xml, mainnode, 0);

  // Free memory
  xml.FreeDoc(xmldoc);

  // Print bounding box
  std::cout << "Bounding box: " << std::endl;
  std::cout << "  x: " << x_min_ << " - " << x_max_ << std::endl;
  std::cout << "  y: " << y_min_ << " - " << y_max_ << std::endl;
  std::cout << "  z: " << z_min_ << " - " << z_max_ << std::endl;
}

void DataParserVTP::PostProcessData() {
  // Rescale the data
  constexpr double kSimSpace = 2000;

  // Compute the center of the bounding box
  const double x_center = (x_max_ + x_min_) / 2;
  const double y_center = (y_max_ + y_min_) / 2;
  const double z_center = (z_max_ + z_min_) / 2;

  // Calculate the scaling factor
  const double x_range = x_max_ - x_min_;
  const double y_range = y_max_ - y_min_;
  const double z_range = z_max_ - z_min_;
  const double max_range = std::max(x_range, std::max(y_range, z_range));
  double scaling_factor = kSimSpace / max_range;

  // Rescale and center the data
  for (auto& point : points_) {
    point[0] = (point[0] - x_center) * scaling_factor;
    point[1] = (point[1] - y_center) * scaling_factor;
    point[2] = (point[2] - z_center) * scaling_factor;
  }

  // Determine the segments data.reserve(points_.size());
  double min_length = std::numeric_limits<double>::max();
  for (size_t i = 0; i < points_.size(); i += 2) {
    // Get the start and end point of the segment
    VesselSegment vs;
    vs.start_position = points_[i];
    vs.end_position = points_[i + 1];

    // Get the length of the segment
    double length = (vs.start_position - vs.end_position).Norm();
    min_length = std::min(min_length, length);

    // Get the radius of the segment
    vs.radius = radii_[i / 2] * scaling_factor;

    // Create the segment
    data.push_back(vs);
  }
  std::cout << "Minimum segment length: " << min_length / scaling_factor
            << std::endl;

  // Restructure the data
  std::vector<Double3> unique_points;
  for (size_t i = 0; i < points_.size(); i++) {
    // Check if the point is already in the list
    bool is_unique = true;
    size_t unique_index{0};
    for (size_t j = 0; j < unique_points.size(); j++) {
      // if (unique_points[j] == points_[i]) {
      if ((unique_points[j] - points_[i]).Norm() / scaling_factor < 2e-6) {
        is_unique = false;
        unique_index = j;
        break;
      }
    }

    // If the point is unique, add it to the list
    if (is_unique) {
      unique_points.push_back(points_[i]);
      unique_index = unique_points.size() - 1;
    }

    // Replace all occurences of i in connectivity_ with unique_index
    for (auto& connection : connectivity_) {
      if (connection == i) {
        connection = unique_index;
      }
    }
  }
  // Rescale all radii
  for (auto& radius : radii_) {
    radius *= scaling_factor;
  }

  // Restucture the the connectivity to a vector of pairs
  std::vector<std::pair<int, int>> connectivity;
  for (size_t i = 0; i < connectivity_.size(); i += 2) {
    auto start_index = connectivity_[i];
    auto end_index = connectivity_[i + 1];
    if (start_index == end_index) {
      Log::Fatal("DataParserVTP", "Start and end index are the same");
    }
    if (start_index > end_index) {
      connectivity.push_back(
          std::make_pair(connectivity_[i + 1], connectivity_[i]));
    } else {
      connectivity.push_back(
          std::make_pair(connectivity_[i], connectivity_[i + 1]));
    }
  }

  // Sort the connectivity vector by the start index
  auto p = GetSortPermutation(connectivity, [](std::pair<int, int> const& a,
                                               std::pair<int, int> const& b) {
    return a.first < b.first;
  });
  ApplyPermutationInPlace(connectivity, p);
  ApplyPermutationInPlace(radii_, p);

  // std::sort(connectivity.begin(), connectivity.end(),
  //           [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
  //             return a.first < b.first;
  //           });

  // Write the connectivity back to the connectivity_ vector
  connectivity_.clear();
  for (const auto& connection : connectivity) {
    connectivity_.push_back(connection.first);
    connectivity_.push_back(connection.second);
  }

  // Print the first 10 entries of the connectivity vector
  std::cout << "Connectivity: " << std::endl;
  for (size_t i = 0; i < 20; i++) {
    std::cout << connectivity_[i] << " ";
  }
  std::cout << std::endl;

  // Verify connectivity_.
  std::vector<int> start_indices;
  std::vector<int> end_indices;
  int max_index = 0;
  for (size_t i = 0; i < connectivity_.size(); i += 2) {
    start_indices.push_back(connectivity_[i]);
    end_indices.push_back(connectivity_[i + 1]);
    // Start and end indices must be different
    if (connectivity_[i] == connectivity_[i + 1]) {
      Log::Fatal("DataParserVTP", "Start and end index are the same");
    }
    max_index = std::max(max_index, connectivity_[i]);
    max_index = std::max(max_index, connectivity_[i + 1]);
  }
  // Verify that the start index is always lower than the end index
  for (size_t i = 0; i < start_indices.size(); i++) {
    if (start_indices[i] >= end_indices[i]) {
      Log::Fatal("DataParserVTP", "Start index is higher than end index.");
    }
  }
  std::vector<int> start_indices_count = std::vector<int>(max_index, 0);
  std::vector<int> end_indices_count = std::vector<int>(max_index, 0);
  for (size_t i = 0; i < start_indices.size(); i++) {
    start_indices_count[start_indices[i]]++;
    end_indices_count[end_indices[i]]++;
  }
  // Print all start and end indices with a count of 2 or more
  for (size_t i = 0; i < start_indices_count.size(); i++) {
    if (start_indices_count[i] > 1) {
      std::cout << "Start index " << i << " has a count of "
                << start_indices_count[i] << std::endl;
    }
    if (end_indices_count[i] > 1) {
      std::cout << "End index " << i << " has a count of "
                << end_indices_count[i] << std::endl;
    }
  }
  // Print the first 10 entries of the start_indices and end_indices vectors
  std::cout << "Start indices: " << std::endl;
  for (size_t i = 0; i < 10; i++) {
    std::cout << start_indices[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "End indices: " << std::endl;
  for (size_t i = 0; i < 10; i++) {
    std::cout << end_indices[i] << " ";
  }
  std::cout << std::endl;

  // For each index < max_index, we determine the first start and end index.
  // The start index must be lower than the end index.
  for (int i = 0; i < max_index; i++) {
    // find the first appearance of i in start_indices with std::find
    auto found_in_start_indices =
        std::find(start_indices.begin(), start_indices.end(), i);
    // find the first appearance of i in end_indices with std::find
    auto found_in_end_indices =
        std::find(end_indices.begin(), end_indices.end(), i);
    // If i is found in both vectors, the start index must be lower than the end
    // index
    if (found_in_start_indices != start_indices.end() &&
        found_in_end_indices != end_indices.end()) {
      // Get the indices where i is found in start_indices and end_indices
      auto start_index =
          std::distance(start_indices.begin(), found_in_start_indices);
      auto end_index = std::distance(end_indices.begin(), found_in_end_indices);
      if (start_index < end_index) {
        Log::Fatal(
            "DataParserVTP",
            "Start index for given index is higher than end index. Index ", i,
            " Start index: ", start_index, " End index: ", end_index);
      }
    }
  }

  std::cout << "Number of unique points: " << unique_points.size() << std::endl;
  std::cout << "Number of original points: " << points_.size() << std::endl;
  points_ = unique_points;
}

void DataParserVTP::RecursivelyParseVTPFile(TXMLEngine& xml,
                                            XMLNodePointer_t node,
                                            Int_t level) {
  // this function display all accessible information about xml node and its
  // children

  // printf("%*c node: %s\n", level, ' ', xml.GetNodeName(node));

  // display namespace
  // XMLNsPointer_t ns = xml.GetNS(node);
  // if (ns != 0)
  //   printf("%*c namespace: %s refer: %s\n", level + 2, ' ',
  //   xml.GetNSName(ns),
  //          xml.GetNSReference(ns));

  // Parse attributes
  XMLAttrPointer_t attr = xml.GetFirstAttr(node);
  while (attr != 0) {
    // printf("%*c attr: %s value: %s\n", level + 2, ' ', xml.GetAttrName(attr),
    //        xml.GetAttrValue(attr));
    std::string attr_name = xml.GetAttrName(attr);
    std::string attr_value = xml.GetAttrValue(attr);
    if (attr_name == "NumberOfLines") {
      std::cout << " Found number of lines!" << std::endl;
      num_lines_ = ParseStringForNumber<size_t>(attr_value);
    } else if (attr_name == "NumberOfPoints") {
      std::cout << "Found number of points!" << std::endl;
      num_points_ = ParseStringForNumber<size_t>(attr_value);
    } else if (attr_name == "Name") {
      if (attr_value == "pressure [mmHg]") {
        std::cout << "Found pressure!" << std::endl;
        pressure_ = ParseString<double>(xml.GetNodeContent(node));
      } else if (attr_value == "R") {
        std::cout << "Found radius!" << std::endl;
        radii_ = ParseString<double>(xml.GetNodeContent(node));
      } else if (attr_value == "G") {
        std::cout << "Found G!" << std::endl;
        g_ = ParseString<double>(xml.GetNodeContent(node));
      } else if (attr_value == "mu") {
        std::cout << "Found mu!" << std::endl;
        mu_ = ParseString<double>(xml.GetNodeContent(node));
      } else if (attr_value == "Coordinates") {
        std::cout << "Found coordinates!" << std::endl;
        auto tmp = ParseString<double>(xml.GetNodeContent(node));
        for (size_t i = 0; i < tmp.size(); i += 3) {
          x_min_ = std::min(x_min_, tmp[i]);
          x_max_ = std::max(x_max_, tmp[i]);
          y_min_ = std::min(y_min_, tmp[i + 1]);
          y_max_ = std::max(y_max_, tmp[i + 1]);
          z_min_ = std::min(z_min_, tmp[i + 2]);
          z_max_ = std::max(z_max_, tmp[i + 2]);
          points_.push_back({tmp[i], tmp[i + 1], tmp[i + 2]});
        }
      } else if (attr_value == "connectivity") {
        std::cout << "Found connectivity!" << std::endl;
        connectivity_ = ParseString<int>(xml.GetNodeContent(node));
      } else if (attr_value == "offsets") {
        std::cout << "Found offsets!" << std::endl;
        offsets_ = ParseString<int>(xml.GetNodeContent(node));
      } else {
        // Do nothing
        ;
      }
    } else {
      // Do nothing
      ;
    }
    attr = xml.GetNextAttr(attr);
  }

  // display all child nodes
  XMLNodePointer_t child = xml.GetChild(node);
  while (child != 0) {
    RecursivelyParseVTPFile(xml, child, level + 2);
    child = xml.GetNext(child);
  }
}

}  // namespace bdm
