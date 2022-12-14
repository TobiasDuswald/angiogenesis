# This script parses the metadata of all simulation output files in a given
# folder. The metadata is saved as a csv file.  The metadata is parsed using the
# json module.  The metadata is stored in a dictionary, where the key is the
# name of the simulation parameter and the value is the value of the simulation
# parameter. We export the metadata as a pandas dataframe, where each row is a
# simulation output file and each column is a simulation parameter.
# usage:
# python parse_metadata.py --folder=output --filter=bdm::SimParam --simplify
# python parse_metadata.py --absfolder=/home/user/output --filter=bdm::SimParam

import pandas as pd
from absl import app
from absl import flags
import json
import os


FLAGS = flags.FLAGS
flags.DEFINE_string(
    "folder", "output", "Folder relative to project source to parse."
)
flags.DEFINE_string(
    "absfolder", None, "Absolute folder to parse. Overrides --folder."
)
flags.DEFINE_string(
    "filter", "bdm::SimParam", "Filter to apply to the metadata."
)
flags.DEFINE_string("filename", "metadata", "Filename for output csv.")
flags.DEFINE_boolean(
    "simplify", False, "Simplify the metadata, drop identical columns."
)


def simplify_metadata(df):
    """
    Drop columns that are identical for all rows.
    """
    # Get the columns that are identical for all rows
    identical_columns = df.apply(lambda x: x.nunique() == 1, axis=0)
    identical_columns = identical_columns[identical_columns == True]
    identical_columns = identical_columns.index
    # Drop the columns
    df = df.drop(columns=identical_columns)
    return df


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
    # Convert to a dataframe
    df = pd.DataFrame(parameters, index=[0])
    # If "output" is in the filename, then strip anything before it
    if "output" in filename:
        filename = filename.split("output")[1]
    # Add the filename
    df["filename"] = filename
    # Set the index to the filename
    df = df.set_index("filename")
    return df


def parse_all_files(files, filter="bdm::SimParam"):
    dfs = []
    for file in files:
        dfs.append(parse_metadata(file, filter=filter))
    df = pd.concat(dfs)
    return df


def search_metadata(folder, search_term="metadata"):
    """
    Find all metadata files in a given folder.
    """
    from pathlib import Path

    files = []
    for path in Path(folder).rglob("*metadata"):
        files.append(str(path.absolute()))
    return files


def main(argv):
    # Get the location of this script
    dir_path = os.path.dirname(os.path.realpath(__file__))

    # Get the location of the output folder
    if FLAGS.absfolder is not None:
        folder = FLAGS.absfolder
    else:
        folder = os.path.join(dir_path, "..", FLAGS.folder)
        folder = os.path.abspath(folder)
    print("Scanning folder: {}".format(folder))

    # Get the location of all metadata files
    files = search_metadata(folder)
    print("Found {} files ..".format(len(files)))

    # Parse all metadata files
    print("Parsing metadata ..")
    df = parse_all_files(files, filter=FLAGS.filter)

    # Simplify the metadata
    if FLAGS.simplify:
        print("Simplifying metadata ..")
        df = simplify_metadata(df)

    # Save the dataframe
    print("Saving metadata ..")
    save_path = os.path.join(
        dir_path, "..", "{}.csv".format(str(FLAGS.filename))
    )
    df.to_csv(save_path, index=True)


if __name__ == "__main__":
    app.run(main)
