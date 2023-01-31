#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"

// Execute the sctipt with `root -l exp_fit.cpp` or open it in ROOT and type
// `.x exp_fit.cpp` in the command line

void exp_fit() {
  // ---------------------------------------------------------------------------
  // 1. Define the data points (see growth_data_combined.txt)
  // ---------------------------------------------------------------------------

  int n = 5;
  double *x = new double[n];
  double *y = new double[n];
  double *e = new double[n];

  x[0] = 7.0;
  x[1] = 14.0;
  x[2] = 23.0;
  x[3] = 29.0;
  x[4] = 34.0;

  // substract 7 from all x values
  for (int i = 0; i < n; i++) {
    x[i] = x[i] - 7;
    x[i] *= 24 * 60;  // convert to hours
  }

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
  graph_errors->GetXaxis()->SetTitle("time [min]");
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
  myc->SaveAs("exp_fit.pdf");

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
            << " [1/min]" << std::endl;

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
  std::cout << "\nThe doubling time is    : " << doubling_time / 24 / 60
            << " +- " << doubling_time_error / 24 / 60 << " [days]"
            << std::endl;
}
