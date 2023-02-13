#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"

// Execute the sctipt with `root -l exp_fit.cpp` or open it in ROOT and type
// `.x exp_fit.cpp` in the command line

// Enum class to select data. The data is given from treatment on, except for
// kPreTreatment, which is the data before treatment.
enum class DataGroup {
  kPreTreatment,
  kGroup1,
  kGroup2,
  kGroup3,
  kGroup4,
  kGroup5,
  kGroup6
};

// Functions to get the data (data hard coded below)
void get_data(const DataGroup group, const std::string &time_unit, int &n,
              double *&x, double *&y, double *&e);

double get_time_conversion_factor(const std::string &time_unit_in,
                                  const std::string &time_unit_out);

void exp_fit() {
  // ---------------------------------------------------------------------------
  // 1. Define the data points (see growth_data_combined.txt)
  // ---------------------------------------------------------------------------

  // const DataGroup group = DataGroup::kGroup6;
  const DataGroup group = DataGroup::kPreTreatment;
  const std::string time_unit = "min";  // "day", "hour", "min", "sec"

  int n{0};
  double *x{nullptr};
  double *y{nullptr};
  double *e{nullptr};

  get_data(group, time_unit, n, x, y, e);

  // ---------------------------------------------------------------------------
  // 2. Fit the data points
  // ---------------------------------------------------------------------------

  TCanvas *myc = new TCanvas("myc", "Exponential fit of growth data");
  myc->SetGrid();

  // Construct the TGraphErrors object from the data points and draw it
  TGraphErrors *graph_errors = new TGraphErrors(n, x, y, 0, e);
  graph_errors->Draw("a*");
  graph_errors->SetMarkerColor(kBlue);
  graph_errors->SetLineColor(kBlue);
  graph_errors->SetTitle("");

  // Label x and y axis
  const std::string xlabel = "time [" + time_unit + "]";
  graph_errors->GetXaxis()->SetTitle(xlabel.c_str());
  graph_errors->GetYaxis()->SetTitle("volume [mm^3]");

  // Define the fitting function and fit the data points. "expo" is a built-in
  // function in ROOT. It corresponds to a exponential function with two
  // parameters: constant and slope. The fitting function is defined as
  // f(x) = exp(constant + slope * x).
  graph_errors->Fit("expo");

  // Access the fit results:
  TF1 *exp = graph_errors->GetFunction("expo");
  exp->SetName("fit");
  exp->SetLineColor(kBlue);
  exp->SetLineWidth(1);
  const double constant = exp->GetParameter(0);
  const double exp_constant = TMath::Exp(constant);
  const double slope = exp->GetParameter(1);
  const double constant_error = exp->GetParError(0);
  const double slope_error = exp->GetParError(1);

  // Plot the expo function with the fit parameters plus and minus the error
  // on the parameters
  TF1 *exp_plus = new TF1("exp_plus", "exp([0] + [1] * x)", x[0], x[n - 1]);
  exp_plus->SetLineColor(kRed);
  exp_plus->SetLineWidth(1);
  exp_plus->SetParameter(0, constant + constant_error);
  exp_plus->SetParameter(1, slope + slope_error);
  exp_plus->Draw("same");

  TF1 *exp_minus = new TF1("exp_minus", "exp([0] + [1] * x)", x[0], x[n - 1]);
  exp_minus->SetLineColor(kRed);
  exp_minus->SetLineWidth(1);
  exp_minus->SetParameter(0, constant - constant_error);
  exp_minus->SetParameter(1, slope - slope_error);
  exp_minus->Draw("same");

  // Add a legend
  TLegend *leg = new TLegend(0.2, 0.7, 0.55, 0.9);
  leg->AddEntry(graph_errors, "data", "p");
  leg->AddEntry(exp, "fit", "l");
  leg->AddEntry(exp_plus, "fit #pm error", "l");
  leg->Draw();

  // Save the plot
  const std::string filename =
      "exp_fit_" + std::to_string(int(group)) + "_" + time_unit + ".pdf";
  myc->SaveAs(filename.c_str());

  // Print interpretation
  std::cout << "\nThe exponential fit is given by: " << std::endl;
  std::cout << "f(x) = exp(" << constant << " + " << slope << " * x)\n"
            << "     = " << exp_constant << " * exp(" << slope << " * x)"
            << std::endl;

  // Give the fit parameters and their errors
  std::cout << "\nThe fit parameters are: " << std::endl;
  std::cout << "constant                : " << constant << " +- "
            << constant_error << " [1]" << std::endl;
  std::cout << "slope                   : " << slope << " +- " << slope_error
            << " [1/" << time_unit << "]" << std::endl;

  // Convert the constant to the initial volume and give the error via error
  // propagation (see wikipedia for error propagation formula)
  const double initial_volume = TMath::Exp(constant);
  const double initial_volume_error = initial_volume * (constant_error);
  std::cout << "\nThe initial volume is   : " << initial_volume << " +- "
            << initial_volume_error << " [mm^3]" << std::endl;

  // Print predicted doubling time and error via error propagation
  const double doubling_time = TMath::Log(2) / slope;
  const double doubling_time_error = doubling_time * (slope_error / slope);
  // Print the result converted to days
  double time_conversion_factor = get_time_conversion_factor(time_unit, "day");
  std::cout << "\nThe doubling time is    : "
            << doubling_time * time_conversion_factor << " +- "
            << doubling_time_error * time_conversion_factor << " [days]"
            << std::endl;
}

void get_data(const DataGroup group, const std::string &time_unit, int &n,
              double *&x, double *&y, double *&e) {
  // Clear data
  delete x;
  delete y;
  delete e;
  // Set to nullptr to avoid dangling pointers
  x = nullptr;
  y = nullptr;
  e = nullptr;

  // Set number of data points
  switch (group) {
    case DataGroup::kPreTreatment:
      n = 5;
      break;

    default:
      n = 14;
      break;
  }

  // Allocate memory for the data
  x = new double[n];
  y = new double[n];
  e = new double[n];
  double time_offset = 0;

  // Time conversion factor
  double time_conversion_factor = get_time_conversion_factor("day", time_unit);

  // Fill the x data
  switch (group) {
    case DataGroup::kPreTreatment:
      x[0] = 7.0;
      x[1] = 14.0;
      x[2] = 23.0;
      x[3] = 29.0;
      x[4] = 34.0;
      time_offset = 7.0;
      // time_offset = 5.0;
      break;

    default:
      x[0] = 35;
      x[1] = 36;
      x[2] = 37;
      x[3] = 40;
      x[4] = 42;
      x[5] = 44;
      x[6] = 47;
      x[7] = 49;
      x[8] = 51;
      x[9] = 54;
      x[10] = 56;
      x[11] = 61;
      x[12] = 63;
      x[13] = 68;
      time_offset = 35.0;
      break;
  }

  // Remove the time offset and convert to the desired time unit
  for (int i = 0; i < n; i++) {
    x[i] -= time_offset;
    x[i] *= time_conversion_factor;
  }

  // Fill the y and e data
  switch (group) {
    case DataGroup::kPreTreatment:
      y[0] = 42.721428571428575;
      y[1] = 57.53809523809524;
      y[2] = 99.06214285714285;
      y[3] = 203.95428571428567;
      y[4] = 337.21738095238095;
      e[0] = 18.69725544897009;
      e[1] = 26.190921461741254;
      e[2] = 53.76778199500163;
      e[3] = 101.30082863029637;
      e[4] = 150.22145339074225;
      break;
    case DataGroup::kGroup1:
      y[0] = 327.55;
      y[1] = 372.13;
      y[2] = 363.14;
      y[3] = 400.36;
      y[4] = 450.91;
      y[5] = 465.51;
      y[6] = 483.39;
      y[7] = 589.67;
      y[8] = 689.38;
      y[9] = 749.47;
      y[10] = 875.36;
      y[11] = 1194.46;
      y[12] = 1218.30;
      y[13] = 1640.11;
      e[0] = 98.02;
      e[1] = 124.06;
      e[2] = 170.50;
      e[3] = 119.16;
      e[4] = 221.79;
      e[5] = 179.67;
      e[6] = 182.00;
      e[7] = 288.69;
      e[8] = 221.18;
      e[9] = 346.45;
      e[10] = 463.20;
      e[11] = 664.49;
      e[12] = 577.52;
      e[13] = 788.76;
      break;
    case DataGroup::kGroup2:
      y[0] = 253.25;
      y[1] = 291.30;
      y[2] = 287.25;
      y[3] = 308.06;
      y[4] = 369.38;
      y[5] = 395.23;
      y[6] = 497.79;
      y[7] = 574.88;
      y[8] = 599.99;
      y[9] = 655.05;
      y[10] = 846.18;
      y[11] = 1165.65;
      y[12] = 1043.77;
      y[13] = 1238.73;
      e[0] = 91.93;
      e[1] = 137.01;
      e[2] = 127.89;
      e[3] = 160.04;
      e[4] = 175.77;
      e[5] = 182.30;
      e[6] = 236.40;
      e[7] = 257.73;
      e[8] = 311.14;
      e[9] = 314.64;
      e[10] = 387.37;
      e[11] = 462.57;
      e[12] = 782.41;
      e[13] = 860.63;
      break;
    case DataGroup::kGroup3:
      y[0] = 424.15;
      y[1] = 387.86;
      y[2] = 392.03;
      y[3] = 326.11;
      y[4] = 337.42;
      y[5] = 309.13;
      y[6] = 308.30;
      y[7] = 331.70;
      y[8] = 318.42;
      y[9] = 330.26;
      y[10] = 324.89;
      y[11] = 362.43;
      y[12] = 385.02;
      y[13] = 382.49;
      e[0] = 99.11;
      e[1] = 65.67;
      e[2] = 80.06;
      e[3] = 96.22;
      e[4] = 134.63;
      e[5] = 142.16;
      e[6] = 156.32;
      e[7] = 181.88;
      e[8] = 188.58;
      e[9] = 189.56;
      e[10] = 182.62;
      e[11] = 264.85;
      e[12] = 273.41;
      e[13] = 373.35;
      break;
    case DataGroup::kGroup4:
      y[0] = 493.29;
      y[1] = 495.42;
      y[2] = 527.12;
      y[3] = 538.87;
      y[4] = 602.62;
      y[5] = 551.35;
      y[6] = 651.08;
      y[7] = 791.66;
      y[8] = 821.88;
      y[9] = 893.61;
      y[10] = 1315.13;
      y[11] = 1681.84;
      y[12] = 1917.22;
      y[13] = 2571.60;
      e[0] = 208.72;
      e[1] = 224.98;
      e[2] = 232.59;
      e[3] = 204.49;
      e[4] = 207.90;
      e[5] = 189.70;
      e[6] = 220.07;
      e[7] = 303.31;
      e[8] = 288.89;
      e[9] = 340.50;
      e[10] = 185.20;
      e[11] = 251.88;
      e[12] = 299.99;
      e[13] = 414.17;
      break;
    case DataGroup::kGroup5:
      y[0] = 238.00;
      y[1] = 227.29;
      y[2] = 207.24;
      y[3] = 163.55;
      y[4] = 148.91;
      y[5] = 141.74;
      y[6] = 137.84;
      y[7] = 122.13;
      y[8] = 91.34;
      y[9] = 59.33;
      y[10] = 69.93;
      y[11] = 47.73;
      y[12] = 70.08;
      y[13] = 71.78;
      e[0] = 109.67;
      e[1] = 135.08;
      e[2] = 129.63;
      e[3] = 84.69;
      e[4] = 78.78;
      e[5] = 79.90;
      e[6] = 71.89;
      e[7] = 88.98;
      e[8] = 78.96;
      e[9] = 78.71;
      e[10] = 90.59;
      e[11] = 58.84;
      e[12] = 65.56;
      e[13] = 61.66;
      break;
    case DataGroup::kGroup6:
      y[0] = 297.32;
      y[1] = 268.52;
      y[2] = 241.84;
      y[3] = 140.07;
      y[4] = 131.73;
      y[5] = 98.76;
      y[6] = 81.14;
      y[7] = 71.51;
      y[8] = 70.43;
      y[9] = 18.53;
      y[10] = 14.69;
      y[11] = 17.41;
      y[12] = 16.75;
      y[13] = 15.52;
      e[0] = 134.83;
      e[1] = 118.93;
      e[2] = 85.53;
      e[3] = 43.16;
      e[4] = 40.31;
      e[5] = 41.35;
      e[6] = 40.78;
      e[7] = 41.88;
      e[8] = 63.39;
      e[9] = 37.05;
      e[10] = 29.37;
      e[11] = 34.83;
      e[12] = 33.49;
      e[13] = 31.05;
      break;
    default:
      std::cout << "Error: Data group not found!" << std::endl;
      break;
  }
}

// Accepted time units: day, hour, min, sec
double get_time_conversion_factor(const std::string &time_unit_in,
                                  const std::string &time_unit_out) {
  // Time conversion factor
  double time_conversion_factor{1.0};
  if (time_unit_in == "day") {
    if (time_unit_out == "day") {
      time_conversion_factor = 1;
    } else if (time_unit_out == "hour") {
      time_conversion_factor = 24.0;
    } else if (time_unit_out == "min") {
      time_conversion_factor = 24.0 * 60.0;
    } else if (time_unit_out == "sec") {
      time_conversion_factor = 24.0 * 60.0 * 60.0;
    } else {
      std::cout << "Error: Unknown time unit. Please use day, hour, min or sec."
                << std::endl;
      exit(1);
    }
  } else if (time_unit_in == "hour") {
    if (time_unit_out == "day") {
      time_conversion_factor = 1.0 / 24.0;
    } else if (time_unit_out == "hour") {
      time_conversion_factor = 1;
    } else if (time_unit_out == "min") {
      time_conversion_factor = 60.0;
    } else if (time_unit_out == "sec") {
      time_conversion_factor = 60.0 * 60.0;
    } else {
      std::cout << "Error: Unknown time unit. Please use day, hour, min or sec."
                << std::endl;
      exit(1);
    }
  } else if (time_unit_in == "min") {
    if (time_unit_out == "day") {
      time_conversion_factor = 1.0 / (24.0 * 60.0);
    } else if (time_unit_out == "hour") {
      time_conversion_factor = 1.0 / 60.0;
    } else if (time_unit_out == "min") {
      time_conversion_factor = 1;
    } else if (time_unit_out == "sec") {
      time_conversion_factor = 60.0;
    } else {
      std::cout << "Error: Unknown time unit. Please use day, hour, min or sec."
                << std::endl;
      exit(1);
    }
  } else if (time_unit_in == "sec") {
    if (time_unit_out == "day") {
      time_conversion_factor = 1.0 / (24.0 * 60.0 * 60.0);
    } else if (time_unit_out == "hour") {
      time_conversion_factor = 1.0 / (60.0 * 60.0);
    } else if (time_unit_out == "min") {
      time_conversion_factor = 1.0 / 60.0;
    } else if (time_unit_out == "sec") {
      time_conversion_factor = 1;
    } else {
      std::cout << "Error: Unknown time unit. Please use day, hour, min or sec."
                << std::endl;
      exit(1);
    }
  }
  return time_conversion_factor;
}
