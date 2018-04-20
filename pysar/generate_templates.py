#!/usr/bin/env python2

import os
import pandas as pd
import sys
import requests
import argparse

inps = None

def cmdLineParse(argv):
    global inps
    parser = argparse.ArgumentParser(description='Generate Processing Template Files',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=None)

    ##### Input
    infile = parser.add_argument_group('File to Generate', 'File to Generate')
    infile.add_argument('--file', dest='file', metavar='FILE', help='file to generate')
    infile.add_argument('--csv', dest='csv', metavar='FILE', help='CSV file containing template data')
    infile.add_argument('--output-dir', dest='output', metavar='FILE', help='directory to output template files to')

    inps = parser.parse_args(argv)
    return inps


#########################################################

def get_google_spreadsheet_as_string(url_id, output_type="csv"):
    response = requests.get('https://docs.google.com/spreadsheet/ccc?key=' + url_id + '&output=' + output_type)
    response.raise_for_status()
    return response.content


def write_file(content, output_file_location = None):
    name = "templateExample.csv"
    if output_file_location is not None:
        name = os.path.join(output_file_location, name)
    with open(name, 'wb') as f:
        f.write(content)
    return name


def get_google_spreadsheet_as_dataframe(url_id, output_file_location, output_type = "csv"):
    content = get_google_spreadsheet_as_string(url_id, output_type)
    loc = write_file(content, output_file_location)
    df = pd.read_csv(loc)
    return df


def get_spreadsheet_as_dataframe(file, output_file_location, output_type = "csv"):
    file_end = file[len(file)-4:len(file)]
    if file_end in ['.csv']:
        df = pd.read_csv(file)
    else:
        df = get_google_spreadsheet_as_dataframe(file, output_file_location, output_type)
    return df


##########################################################

def generate_template_file(names, subnames, column_vals):
    line_breaker = "\n" + "#"*20 + "\n\n"
    output_file = ""
    last_named_column = ""
    base = ""

    for i in range(len(names)):
        # Get set base name.
        if type(names[i]) != str:
            base = last_named_column
        else:
            # New Column
            if base != names[i] and i > 1 and "Unnamed:" not in names[i]:
                output_file += line_breaker
            base = names[i]
            last_named_column = base

        # Get subname
        if type(subnames[i]) == str:
            subname = "." + subnames[i]
        else:
            subname = ""

        # Get Value
        if type(column_vals[i]) ==  str:
            value = column_vals[i]
        else:
            continue

        # Need to create name in order to format string properly so that the "=" symbol is at the same location
        name = base + subname
        output_line = "{0:35} = {1:10}\n".format(name, value)
        output_file += output_line
    return output_file


def generate_template_files(df):
    names = list(df["Name"])
    subnames = list(df["Subname"])
    columns = list(df.columns)
    output_files = {}

    if inps.file is not None:
        file_base = inps.file
        output_files[file_base] = generate_template_file(names, subnames, list(df[inps.file]))
    else:
        for i, col_name in enumerate(columns[2:]):
            file_base = col_name
            output_files[file_base] = generate_template_file(names, subnames, list(df[col_name]))
    return output_files


def generate_and_save_template_files(df, output_location):
    # Create output directory if it doesn't exist
    if not os.path.isdir(output_location):
        os.mkdir(output_location)

    files_to_save = generate_template_files(df)

    for key, value in files_to_save.items():
        # Print contents for debugging purposes
        print(key)
        print()
        for i in value.split('\n'):
            print(i)

        with open(os.path.join(output_location, key + ".template"), "w") as f:
            f.write(value)


def generate_and_save_template_files_from_file(csv_file, output_location):
    df = pd.read_csv(csv_file)
    generate_and_save_template_files(df, output_location)


def generate_and_save_template_files_from_dataframe(df, output_location):
    generate_and_save_template_files(df, output_location)


# TODO: Properly name variables
# If output and input directories are declared, use them
if __name__ == "__main__":
    inps = cmdLineParse(sys.argv[1:])

    csv_file = "1Mvxf-O1NV-TJK9Ax7vWTvZ8q9jWx-GQD4y5WGgTOcMc"
    output_location = "/Users/Joshua/Desktop/output_dir"

    print(inps.output)

    if inps.csv:
        csv_file = inps.csv

    if inps.output:
        output_location = inps.output

    df = get_spreadsheet_as_dataframe(csv_file, output_location)
    generate_and_save_template_files_from_dataframe(df, output_location)