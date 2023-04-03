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
  for (size_t i = 0; i < points_.size(); i += 2) {
    // Get the start and end point of the segment
    VesselSegment vs;
    vs.start_position = points_[i];
    vs.end_position = points_[i + 1];

    // Get the radius of the segment
    vs.radius = radii_[i / 2] * scaling_factor;

    // Create the segment
    data.push_back(vs);
  }
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
