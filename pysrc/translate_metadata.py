# This scipts loads the metadata from the metadata file and translates it to the
# latex notation of the paper.

import os
from absl import app
from absl import flags
import json

FLAGS = flags.FLAGS
flags.DEFINE_string(
    "filter", "bdm::SimParam", "Filter to apply to the metadata."
)
flags.DEFINE_string("inputfile", "metadata", "Input file.")
flags.DEFINE_string("outputfile", "metadata.tex", "Filename for output csv.")

translation_table = {
    "action_radius_factor": r"r_a",
    "adhesion_scale_parameter": r"c_a",
    "alpha_H_D_DOX": r"not used",
    "alpha_H_D_TRA": r"not used",
    "alpha_Q_D_N": r"a_{Q \rightarrow D}",
    "alpha_Q_SG2_N": r"c_{Q \rightarrow SG2}",
    "alpha_Q_SG2_TRA": r"\lambda_{Q \rightarrow SG2}",
    "alpha_SG2_D_DOX": r"c_{SG2 \rightarrow D}",
    "alpha_SG2_SG2_DOX": r"c_{SG2 \rightarrow SG2}",
    "apical_growth_gradient_weight": r"w_1",
    "apical_growth_old_weight": r"w_2",
    "apical_growth_random_weight": r"w_3",
    "apical_growth_speed": r"not used",
    "base_rate_H_D": r"r_{H \rightarrow D}",
    "cell_nuclear_radius": r"r_n",
    "cell_radius": r"r_p",
    "cell_radius_sigma": r"\sigma_{r_p}",
    "decay_rate_dox": r"\lambda_{d}",
    "decay_rate_nutrients": r"\lambda_{n}",
    "decay_rate_tra": r"\lambda_{t}",
    "decay_rate_vegf": r"\lambda_{v}",
    "default_vessel_length": r"not used",
    "diffusion_dox": r"D_{d}",
    "diffusion_nutrients": r"D_{n}",
    "diffusion_resolution_dox": r"(l/h_{d})",
    "diffusion_resolution_nutrients": r"(l/h_{n})",
    "diffusion_resolution_tra": r"(l/h_{t})",
    "diffusion_resolution_vegf": r"(l/h_{v})",
    "diffusion_tra": r"D_{t}",
    "diffusion_vegf": r"D_{v}",
    "dox_consumption_rate_tcell": r"\alpha_{d}",
    "dox_supply_rate_vessel": r"\beta_{d}",
    "duration_apoptosis": r"not used",
    "duration_cell_cycle": r"T_{SG2}",
    "duration_growth_phase": r"T_{G1}",
    "gamma_H_D_DOX": r"not used",
    "gamma_H_D_TRA": r"not used",
    "gamma_Q_D_N": r"gamma\_Q\_D\_N (not present in paper)",
    "gamma_SG2_D_DOX": r"gamma\_SG2\_D\_DOX (not present in paper)",
    "gamma_SG2_SG2_DOX": r"gamma\_SG2\_SG2\_DOX (not present in paper)",
    "hypoxic_threshold": r"u_n^H",
    "initial_concentration_dox": r"u_n(t=0)",
    "initial_concentration_nutrients": r"u_n(t=0)",
    "initial_concentration_tra": r"u_t(t=0)",
    "initial_concentration_vegf": r"u_v(t=0)",
    "k_H_D_DOX": r"not used",
    "k_H_D_TRA": r"not used",
    "k_Q_D_N": r"k_{Q \rightarrow D}",
    "k_SG2_D_DOX": r"not used",
    "k_SG2_SG2_DOX": r"not used",
    "max_speed": r"not used",
    "min_dist_to_bifurcation": r"d_{branch}",
    "min_dist_to_tip_cell": r"d_{tip}",
    "nutrient_consumption_rate_tcell": r"\alpha_{n}",
    "nutrient_supply_rate_vessel": r"\beta_{n}",
    "repulsive_scale_parameter": r"c_r",
    "secretion_rate_vegf": r"not used",
    "sprouting_probability": r"p_{s,rate}",
    "threshold_H_D_DOX": r"u_d^{H \rightarrow D}",
    "threshold_H_D_TRA": r"u_t^{H \rightarrow D}",
    "threshold_Q_D_N": r"u_n^{Q \rightarrow D}",
    "threshold_Q_SG2_N": r"u_n^{Q \rightarrow SG2}",
    "threshold_SG2_D_DOX": r"u_d^{SG2 \rightarrow D}",
    "threshold_SG2_SG2_DOX": r"u_d^{SG2 \rightarrow SG2}",
    "total_sim_time": r"T",
    "tra_consumption_rate_tcell": r"\alpha_{t}",
    "tra_supply_rate_vessel": r"\beta_{t}",
    "uptake_rate_glucose": r"not used",
    "vegf_consumption_rate_vessel": r"\beta_{v}",
    "vegf_grad_threshold_apical_growth": r"\nabla u_v^{thres}",
    "vegf_supply_rate_tcell": r"\alpha_{v}",
    "vegf_threshold_sprouting": r"u_v^{thres}",
    "viscosity": r"\eta (not present in paper)",
}


def parse_metadata(filename, filter="bdm::SimParam"):
    with open(filename, "r", newline="\n") as f:
        lines = f.readlines()
    # Condense the lines into a single string
    lines = "".join(lines)
    # Find the start of the metadata, e.g. first occurance of "{"
    start = lines.find("{")
    # Find the end of the metadata, e.g. last occurance of "}"
    end = lines.rfind("}")
    # Extract the metadata
    metadata = lines[start : end + 1]
    parameters = json.loads(metadata)
    parameters = parameters[filter]
    return parameters


def translate_dictionary(parameters):
    global translation_table
    keys_to_remove = ["not used"]
    translated = {}
    for key, value in parameters.items():
        if key in translation_table:
            if translation_table[key] in keys_to_remove:
                continue
            translated[translation_table[key]] = value
    return translated


def write_dictionary_to_file(parameters, filename):
    with open(filename, "w", newline="\n") as f:
        for key, value in parameters.items():
            f.write(r"${} = {}$,".format(key, value))
            f.write("\n")
    # Replace the last comma in the file with a period
    with open(filename, "rb+") as f:
        f.seek(-2, os.SEEK_END)
        f.truncate()
        f.write(b".")


def main(argv):
    # Check if the input file exists
    if FLAGS.inputfile is not None:
        filename = FLAGS.inputfile
        if not os.path.isfile(filename):
            print("Input file not found: " + filename)
            return
        else:
            print("Input file: " + filename)

    # Load the metadata
    metadata = parse_metadata(filename, FLAGS.filter)

    # Translate the metadata
    translated = translate_dictionary(metadata)

    # Write the metadata to a file
    write_dictionary_to_file(translated, FLAGS.outputfile)


if __name__ == "__main__":
    app.run(main)
