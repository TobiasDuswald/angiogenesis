// -----------------------------------------------------------------------------
//
// Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <array>
#include "core/stdfilesystem.h"
#include "util/data_parser.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

TEST(VTP_Parser, ParseString) {
  std::string line = "1.1 2.2 3.3 4.4 5.5 6.6 7.7 8.8 9.9 10.10 11.11 12.12 \n";
  std::vector<double> expected_double = {1.1, 2.2, 3.3, 4.4,   5.5,   6.6,
                                         7.7, 8.8, 9.9, 10.10, 11.11, 12.12};
  std::vector<int> expected_int = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

  std::vector<double> result_double = ParseString<double>(line);
  std::vector<int> result_int = ParseString<int>(line);

  // Compare sizes of the vectors
  ASSERT_EQ(expected_double.size(), result_double.size());
  ASSERT_EQ(expected_int.size(), result_int.size());

  // Compare all components of the vectors
  for (size_t i = 0; i < expected_double.size(); i++) {
    EXPECT_EQ(expected_double[i], result_double[i]);
  }
  for (size_t i = 0; i < expected_int.size(); i++) {
    EXPECT_EQ(expected_int[i], result_int[i]);
  }
}

TEST(VTP_Parser, ConstructUniquePoints) {
  // Define input points and connectivity
  const std::vector<Double3> points = {{1, 1, 1}, {2, 2, 2}, {2, 2, 2},
                                       {3, 3, 3}, {2, 2, 2}, {4, 4, 4},
                                       {4, 4, 4}, {5, 5, 5}};
  const std::vector<int> connectivity = {0, 1, 2, 3, 4, 5, 6, 7};

  // Define expected output
  const std::vector<Double3> expected_points = {
      {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}, {5, 5, 5}};
  const std::vector<int> expected_connectivity = {0, 1, 1, 2, 1, 3, 3, 4};

  // Construct the unique points and connectivity
  std::vector<Double3> result_points;
  std::vector<int> result_connectivity;
  ConstructUniquePoints(points, result_points, connectivity,
                        result_connectivity);

  // Compare sizes of the vectors
  ASSERT_EQ(expected_points.size(), result_points.size());
  ASSERT_EQ(expected_connectivity.size(), result_connectivity.size());

  // Compare all components of the vectors
  for (size_t i = 0; i < expected_points.size(); i++) {
    for (size_t j = 0; j < 3; j++) {
      EXPECT_DOUBLE_EQ(expected_points[i][j], result_points[i][j]);
    }
  }
}

TEST(VTP_Parser, ConstructLines) {
  // Input
  const std::vector<int> connectivity = {0, 1, 2, 3, 4, 5, 6, 7};

  // Expected output
  const std::vector<std::pair<int, int>> expected_lines = {
      {0, 1}, {2, 3}, {4, 5}, {6, 7}};

  // Construct the lines
  auto result_lines = ConstructLines(connectivity);

  // Compare sizes of the vectors
  ASSERT_EQ(expected_lines.size(), result_lines.size());

  // Compare all components of the vectors
  for (size_t i = 0; i < expected_lines.size(); i++) {
    EXPECT_EQ(expected_lines[i].first, result_lines[i].first);
    EXPECT_EQ(expected_lines[i].second, result_lines[i].second);
  }
}

TEST(VTP_Parser, VerifyStartingLines) {
  // Input
  const std::vector<std::pair<int, int>> lines = {{0, 1}, {1, 2}, {3, 4},
                                                  {5, 6}, {7, 8}, {7, 9}};
  const std::vector<int> start_lines = {0, 2, 4};

  // Expected output
  const std::vector<bool> verification = {true, false, false};

  // Verify the starting lines
  auto result_verification = VerifyStartingLines(start_lines, lines);

  // Compare sizes of the vectors
  ASSERT_EQ(verification.size(), result_verification.size());
}

TEST(VTP_Parser, AdjustStartingLines) {
  // Input
  std::vector<std::pair<int, int>> lines = {{1, 0}, {1, 2}, {3, 4}, {4, 5},
                                            {5, 6}, {7, 8}, {7, 9}};
  const std::vector<int> start_lines = {0, 2, 5};

  // Expected output
  const bool expected_result = true;
  const std::vector<std::pair<int, int>> new_lines = {
      {0, 1}, {1, 2}, {3, 4}, {4, 5}, {5, 6}, {8, 7}, {7, 9}};

  // Verify the starting lines
  auto result_verification = AdjustStartingLines(start_lines, lines);

  // Compare sizes of the vectors
  ASSERT_TRUE(expected_result);
  ASSERT_EQ(new_lines.size(), lines.size());

  // Compare all components of the vectors
  for (size_t i = 0; i < new_lines.size(); i++) {
    EXPECT_EQ(new_lines[i].first, lines[i].first);
    EXPECT_EQ(new_lines[i].second, lines[i].second);
  }
}

TEST(VTP_Parser, RestructureToTree) {
  // Input
  const std::vector<std::pair<int, int>> connectivity = {
      {0, 1},   {1, 2}, {11, 12}, {3, 2}, {98, 99}, {13, 12}, {14, 12},
      {96, 97}, {4, 3}, {14, 15}, {5, 3}, {16, 15}, {5, 6}};
  const std::vector<int> start_lines = {0, 2};

  // Expected output
  const std::vector<std::pair<int, int>> expected_tree = {
      {0, 1}, {11, 12}, {1, 2},   {12, 13}, {12, 14}, {2, 3},  {14, 15},
      {3, 4}, {3, 5},   {15, 16}, {5, 6},   {98, 99}, {96, 97}};
  const std::vector<int> expected_permutation = {0, 2,  1,  5,  6, 3, 9,
                                                 8, 10, 11, 12, 4, 7};

  // Compute the tree
  std::vector<std::pair<int, int>> tree;
  std::vector<int> permutation;
  RestructureToTree(start_lines, connectivity, tree, permutation);

  // Compare sizes of the vectors
  ASSERT_EQ(expected_tree.size(), tree.size());

  // Compare all components of the vectors
  for (size_t i = 0; i < expected_tree.size(); i++) {
    EXPECT_EQ(expected_tree[i].first, tree[i].first);
    EXPECT_EQ(expected_tree[i].second, tree[i].second);
  }

  // Compare sizes of the vectors
  ASSERT_EQ(expected_permutation.size(), permutation.size());

  // Compare all components of the vectors
  for (size_t i = 0; i < expected_permutation.size(); i++) {
    EXPECT_EQ(expected_permutation[i], permutation[i]);
  }
}

TEST(VTP_Parser, ExtractNumericValueFromString) {
  std::string line = "This text should be ignored 123.45 and this too";
  double expected_double = 123.45;
  int expected_int = 123;

  double result_double = ParseStringForNumber<double>(line);
  int result_int = ParseStringForNumber<int>(line);

  EXPECT_EQ(expected_double, result_double);
  EXPECT_EQ(expected_int, result_int);
}

TEST(VTP_Parser, ReadData) {
  // Define expected values (first few numbers of the vectors)
  constexpr size_t kNumLines = 362;
  constexpr size_t kNumPoints = 724;
  constexpr std::array<double, 12> kPr = {29.4699, 0,       0, 0, 60.5525, 0,
                                          0,       59.9974, 0, 0, 27.5873, 0};
  constexpr std::array<double, 12> kR = {7.67469e-06, 7.73683e-06, 1.36467e-05,
                                         1.11149e-05, 1.26339e-05, 1.00491e-05,
                                         5.4898e-06,  1.0232e-05,  1.29626e-05,
                                         6.24953e-06, 1.29063e-05, 1.27413e-05};
  constexpr std::array<double, 12> kG = {3.28699e-11, 2.3735e-11,  1.11621e-11,
                                         1.58124e-10, 1.13404e-10, 2.27996e-11,
                                         2.06982e-12, 1.54049e-11, 2.36403e-10,
                                         2.49727e-12, 6.03061e-11, 2.39504e-11};
  constexpr std::array<double, 12> kMu = {
      0.00399706, 0.00396984, 0.00264072, 0.00299459, 0.00275786, 0.00321928,
      0.00538128, 0.00317653, 0.00271677, 0.00478268, 0.00272359, 0.00274409};
  constexpr std::array<double, 36> kPointsCoordinates = {
      0.0003185,   4.725e-05,   8.75e-07,    0.00031675,  4.8125e-05,
      1.1375e-05,  0.00031675,  4.8125e-05,  1.1375e-05,  0.0003185,
      4.4625e-05,  2.625e-05,   0.000580125, 8.6625e-05,  1.75e-06,
      0.000623,    0.000184625, 0.0004655,   0.000623,    0.000184625,
      0.0004655,   0.000622125, 0.000196875, 0.000469875, 0.000623,
      0.000184625, 0.0004655,   0.000616875, 0.000184625, 0.000497875,
      0.0002905,   0.000382375, 1.75e-06,    0.0002765,   0.00038675,
      5.6e-05};
  constexpr std::array<int, 12> kConnectivity = {0, 1, 2, 3, 4,  5,
                                                 6, 7, 8, 9, 10, 11};
  constexpr std::array<int, 12> kOffsets = {2,  4,  6,  8,  10, 12,
                                            14, 16, 18, 20, 22, 24};
  // Transform kPointsCoordinates to a vector of 3D points
  std::vector<Double3> points_coordinates;
  for (size_t i = 0; i < kPointsCoordinates.size(); i += 3) {
    points_coordinates.push_back({kPointsCoordinates[i],
                                  kPointsCoordinates[i + 1],
                                  kPointsCoordinates[i + 2]});
  }

  // Define parser & parse data
  DataParserVTP parser;
  parser.ParseData("data/network.vtp");

  // Get all the parsed data
  const auto& num_lines = parser.GetNumLines();
  const auto& num_points = parser.GetNumPoints();
  const auto& pressure = parser.GetPressure();
  const auto& radii = parser.GetRadii();
  const auto& g = parser.GetG();
  const auto& mu = parser.GetMu();
  const auto& points = parser.GetPoints();
  const auto& connectivity = parser.GetConnectivity();
  const auto& offsets = parser.GetOffsets();

  // Check if all sizes are correct
  // Line data
  EXPECT_EQ(num_lines, kNumLines);
  EXPECT_EQ(radii.size(), kNumLines);
  EXPECT_EQ(g.size(), kNumLines);
  EXPECT_EQ(mu.size(), kNumLines);
  EXPECT_EQ(offsets.size(), kNumLines);
  // Point Data
  EXPECT_EQ(num_points, kNumPoints);
  EXPECT_EQ(pressure.size(), kNumPoints);
  EXPECT_EQ(points.size(), kNumPoints);
  EXPECT_EQ(connectivity.size(), kNumPoints);

  // Compare the first 12 values of the vectors
  for (size_t i = 0; i < 12; i++) {
    EXPECT_EQ(pressure[i], kPr[i]);
    EXPECT_EQ(radii[i], kR[i]);
    EXPECT_EQ(g[i], kG[i]);
    EXPECT_EQ(mu[i], kMu[i]);
    EXPECT_EQ(points[i][0], points_coordinates[i][0]);
    EXPECT_EQ(points[i][1], points_coordinates[i][1]);
    EXPECT_EQ(points[i][2], points_coordinates[i][2]);
    EXPECT_EQ(connectivity[i], kConnectivity[i]);
    EXPECT_EQ(offsets[i], kOffsets[i]);
  }

  // Not part of test, just some debugging
  parser.PostProcessData();
  parser.PlotHistograms(fs::current_path());
}

}  // namespace bdm
