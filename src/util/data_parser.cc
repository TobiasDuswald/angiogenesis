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
  // ---------------------------------------------------------------------
  // 1 . Rescale & shift the data
  // ---------------------------------------------------------------------
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

  // Rescale all radii
  for (auto& radius : radii_) {
    radius *= scaling_factor;
  }

  // ---------------------------------------------------------------------
  // 2 . Create the vessel segments (old method)
  // ---------------------------------------------------------------------

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
    vs.radius = radii_[i / 2];

    // Create the segment
    data.push_back(vs);
  }
  std::cout << "Minimum segment length: " << min_length / scaling_factor
            << std::endl;

  // ---------------------------------------------------------------------
  // 3 . Determine the set of unique points; adapt connectivity accordingly
  // ---------------------------------------------------------------------

  // Restructure the data
  std::vector<Double3> unique_points;
  std::vector<int> unique_connectivity;
  ConstructUniquePoints(points_, unique_points, connectivity_,
                        unique_connectivity);
  connectivity_ = unique_connectivity;

  // ---------------------------------------------------------------------
  // 4 . Restructure the connectivity to tree like structure
  // ---------------------------------------------------------------------

  // 1 .Restructure the data format of connectivity to a vector of pairs
  std::cout << "<DataParserVTP> Change data format for connectivity data..."
            << std::endl;
  auto connectivity = ConstructLines(connectivity_);

  // 2. Verify the starting points, e.g. make sure that the orientation is
  // correct
  std::cout << "<DataParserVTP> Verifying starting points..." << std::endl;
  AdjustStartingLines(starting_lines_, connectivity);

  // 3 . Adapt connectivity such that it is a tree like structure, e.g. the
  // starting lines are the root lines and all other lines trace back to them.
  std::cout << "<DataParserVTP> Restructuring connectivity data to tree..."
            << std::endl;

  std::vector<std::pair<int, int>> tree;
  std::vector<int> permutation;
  RestructureToTree(starting_lines_, connectivity, tree, permutation);

  // std::vector<int> root_line_endpoints;
  // root_line_endpoints.clear();
  // for (const auto& line : starting_lines_) {
  //   root_line_endpoints.push_back(connectivity[line].second);
  // }
  // std::vector<int> next_root_line_endpoints;
  // std::vector<std::pair<int, int>> connectivity_restructured;
  // std::vector<bool> connectivity_indices_visited(connectivity.size(), false);
  // std::vector<int> permutation(connectivity.size(), -1);

  // // 3.1 Add the starting lines to the connectivity_restructured and label
  // them
  // // as visited. Keep track of the permutation of the connectivity indices.
  // int ctr = 0;
  // for (const auto& line : starting_lines_) {
  //   connectivity_restructured.push_back(connectivity[line]);
  //   connectivity_indices_visited[line] = true;
  //   permutation[ctr] = line;
  //   ctr++;
  // }

  // // 3.2 Iterate over the connectivity and restructure it
  // while (!root_line_endpoints.empty()) {
  //   for (size_t i = 0; i < connectivity.size(); i++) {
  //     // Ignore already visited lines
  //     if (connectivity_indices_visited[i]) {
  //       continue;
  //     }
  //     // Check if the connection either starts or ends at a root line
  //     endpoint,
  //     // if so, add it to the restructured connectivity. In case it is the
  //     end,
  //     // we swap the start and end point. The updated end point is added to
  //     the
  //     // list of the next root line endpoints for the next iteration.
  //     for (const auto& root_endpoint : root_line_endpoints) {
  //       if (connectivity[i].first == root_endpoint) {
  //         connectivity_restructured.push_back(connectivity[i]);
  //         connectivity_indices_visited[i] = true;
  //         next_root_line_endpoints.push_back(connectivity[i].second);
  //         permutation[ctr] = i;
  //         ctr++;
  //       } else if (connectivity[i].second == root_endpoint) {
  //         connectivity_restructured.push_back(
  //             std::make_pair(connectivity[i].second, connectivity[i].first));
  //         connectivity_indices_visited[i] = true;
  //         next_root_line_endpoints.push_back(connectivity[i].first);
  //         permutation[ctr] = i;
  //         ctr++;
  //       }
  //     }
  //   }
  //   root_line_endpoints = next_root_line_endpoints;
  //   next_root_line_endpoints.clear();
  // }

  // // 3.3 Add all unvisited lines to the restructured connectivity
  // for (size_t i = 0; i < connectivity.size(); i++) {
  //   if (!connectivity_indices_visited[i]) {
  //     connectivity_restructured.push_back(connectivity[i]);
  //     permutation[ctr] = i;
  //     ctr++;
  //   }
  // }

  // ----------------------------------------------------------------------
  // 4. Feedback the restructured connectivity to the connectivity_ vector
  //    and adapt the radii_ vector accordingly
  // ----------------------------------------------------------------------

  // 4.1 Copy the restructured connectivity to the connectivity_ vector
  connectivity_.clear();
  // for (const auto& connection : connectivity_restructured) {
  for (const auto& connection : tree) {
    connectivity_.push_back(connection.first);
    connectivity_.push_back(connection.second);
  }

  // 4.2 Restructure the radii vector according to the permutation
  std::vector<double> radii_restructured;
  for (const auto& index : permutation) {
    radii_restructured.push_back(radii_[index]);
  }
  radii_ = radii_restructured;

  // Print the first 10 entries of the connectivity vector
  std::cout << "Connectivity: " << std::endl;
  for (size_t i = 0; i < 20; i++) {
    std::cout << connectivity_[i] << " ";
  }
  std::cout << std::endl;

  // ----------------------------------------------------------------------
  // 5 . Verify connectivity_
  // ----------------------------------------------------------------------

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

void ConstructUniquePoints(const std::vector<Double3>& points,
                           std::vector<Double3>& unique_points,
                           const std::vector<int>& connectivity,
                           std::vector<int>& unique_connectivity) {
  // Clear the unique points vector
  unique_points.clear();
  // Clear and resize the unique connectivity vector
  unique_connectivity.clear();
  unique_connectivity.resize(connectivity.size());
  // Iterate over all points and add them to the unique points vector if they
  // are not already in it
  for (size_t i = 0; i < points.size(); i++) {
    int unique_index = 0;
    auto already_in_unique_points = false;
    for (size_t j = 0; j < unique_points.size(); j++) {
      if ((points[i] - unique_points[j]).Norm() < 1e-6) {
        already_in_unique_points = true;
        unique_index = j;
        break;
      }
    }
    if (!already_in_unique_points) {
      unique_points.push_back(points[i]);
      unique_index = unique_points.size() - 1;
    }
    for (size_t j = 0; j < connectivity.size(); j++) {
      if (connectivity[j] == i) {
        unique_connectivity[j] = unique_index;
      }
    }
  }
}

std::vector<std::pair<int, int>> ConstructLines(
    const std::vector<int>& connectivity) {
  // Verify that the connectivity vector has an even number of elements
  if (connectivity.size() % 2 != 0) {
    Log::Fatal("ConstructLines",
               "Connectivity vector has an odd number of elements!");
  }
  std::vector<std::pair<int, int>> lines;
  for (size_t i = 0; i < connectivity.size(); i += 2) {
    lines.push_back({connectivity[i], connectivity[i + 1]});
  };
  return lines;
}

std::vector<bool> VerifyStartingLines(
    const std::vector<int>& starting_lines,
    const std::vector<std::pair<int, int>>& lines) {
  // Clear and resize the is_valid vector
  std::vector<bool> is_valid;
  is_valid.resize(starting_lines.size(), true);

  // Get the start and end indices of all starting lines
  std::vector<int> start_indices;
  std::vector<int> end_indices;
  for (const auto& line : starting_lines) {
    start_indices.push_back(lines[line].first);
    end_indices.push_back(lines[line].second);
  }

  // Count the number of occurences of each start and end index
  std::vector<int> start_index_count(start_indices.size(), 0);
  std::vector<int> end_index_count(end_indices.size(), 0);
  for (const auto& line : lines) {
    for (size_t i = 0; i < start_indices.size(); i++) {
      if (line.first == start_indices[i] || line.second == start_indices[i]) {
        start_index_count[i]++;
      }
    }
    for (size_t i = 0; i < end_indices.size(); i++) {
      if (line.first == end_indices[i] || line.second == end_indices[i]) {
        end_index_count[i]++;
      }
    }
  }

  // A line is valid if the start index occurs exactly once and the end index
  // occurs at least twice
  for (size_t i = 0; i < is_valid.size(); i++) {
    if (start_index_count[i] != 1 || end_index_count[i] < 2) {
      is_valid[i] = false;
    }
  }

  return is_valid;
}

bool AdjustStartingLines(const std::vector<int>& starting_lines,
                         std::vector<std::pair<int, int>>& connectivity) {
  // Call back to VerifyStartingLines to verify the starting lines
  auto is_valid = VerifyStartingLines(starting_lines, connectivity);

  // Print validity of starting lines
  std::cout << "<AdjustStartingLines> Validity of starting lines:" << std::endl;
  for (size_t i = 0; i < is_valid.size(); i++) {
    std::cout << "  Line " << starting_lines[i] << ": " << is_valid[i]
              << std::endl;
  }

  // Swap the connectivity of invalid starting lines
  for (size_t i = 0; i < is_valid.size(); i++) {
    if (!is_valid[i]) {
      std::cout << "<AdjustStartingLines> Swapping connectivity of line "
                << starting_lines[i] << std::endl;
      std::swap(connectivity[starting_lines[i]].first,
                connectivity[starting_lines[i]].second);
    }
  }

  // Call back to VerifyStartingLines to verify the starting lines
  is_valid = VerifyStartingLines(starting_lines, connectivity);

  // Print validity of starting lines
  std::cout << "<AdjustStartingLines> Validity of starting lines:" << std::endl;
  for (size_t i = 0; i < is_valid.size(); i++) {
    std::cout << "  Line " << starting_lines[i] << ": " << is_valid[i]
              << std::endl;
  }

  // Return true if all starting lines are valid
  for (const auto& valid : is_valid) {
    if (!valid) {
      return false;
    }
  }
  return true;
}

void RestructureToTree(const std::vector<int>& starting_lines,
                       const std::vector<std::pair<int, int>>& connectivity,
                       std::vector<std::pair<int, int>>& tree,
                       std::vector<int>& permutation) {
  std::vector<int> root_line_endpoints;
  root_line_endpoints.clear();
  for (const auto& line : starting_lines) {
    root_line_endpoints.push_back(connectivity[line].second);
  }
  std::vector<int> next_root_line_endpoints;
  tree.clear();
  std::vector<bool> connectivity_indices_visited(connectivity.size(), false);
  permutation.resize(connectivity.size(), -1);

  // 3.1 Add the starting lines to the tree and label
  // them as visited. Keep track of the permutation of the connectivity indices.
  int ctr = 0;
  for (const auto& line : starting_lines) {
    tree.push_back(connectivity[line]);
    connectivity_indices_visited[line] = true;
    permutation[ctr] = line;
    ctr++;
  }

  std::cout << "<RestructureToTree> Iterating over connectivity...\n";
  int loop_cntr = 0;
  // 3.2 Iterate over the connectivity and restructure it
  while (!root_line_endpoints.empty()) {
    for (size_t i = 0; i < connectivity.size(); i++) {
      // Ignore already visited lines
      if (connectivity_indices_visited[i]) {
        continue;
      }
      // Check if the connection either starts or ends at a root line
      // endpoint, if so, add it to the restructured connectivity. In case it is
      // the end, we swap the start and end point. The updated end point is
      // added to the list of the next root line endpoints for the next
      // iteration.
      for (const auto& root_endpoint : root_line_endpoints) {
        if (connectivity[i].first == root_endpoint) {
          tree.push_back(connectivity[i]);
          connectivity_indices_visited[i] = true;
          next_root_line_endpoints.push_back(connectivity[i].second);
          permutation[ctr] = i;
          ctr++;
        } else if (connectivity[i].second == root_endpoint) {
          tree.push_back(
              std::make_pair(connectivity[i].second, connectivity[i].first));
          connectivity_indices_visited[i] = true;
          next_root_line_endpoints.push_back(connectivity[i].first);
          permutation[ctr] = i;
          ctr++;
        }
      }
    }
    if (loop_cntr % 10 == 0) {
      std::cout << "  Iteration: " << loop_cntr << std::endl;
    }
    loop_cntr++;
    root_line_endpoints = next_root_line_endpoints;
    next_root_line_endpoints.clear();
  }
  std::cout << "  Iteration: " << loop_cntr << "\n  Done." << std::endl;

  // 3.3 Add all unvisited lines to the restructured connectivity
  for (size_t i = 0; i < connectivity.size(); i++) {
    if (!connectivity_indices_visited[i]) {
      tree.push_back(connectivity[i]);
      permutation[ctr] = i;
      ctr++;
    }
  }
}

}  // namespace bdm
