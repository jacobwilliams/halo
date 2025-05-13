"""
    Read a JSON file from the JPL database and convert it to the format required by Halo.

    See: https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html
"""

import os
import sys
import json

def convert_jpl_database(input_file: str, output_file: str):
    """Convert JPL database to Halo format.

    Args:
        input_file (str): JPL database file
        output_file (str): output file (required by Halo)
    """

    with open(input_file, 'r') as f:
        data_in = json.load(f)

    data_out = {}
    data_out['lstar'] = data_in['system']['lunit']
    data_out['tstar'] = data_in['system']['tunit']

    # "fields":["x","y","z","vx","vy","vz","jacobi","period","stability"]
    data_out['jc'] = []
    data_out['period'] = []
    data_out['x0'] = []
    data_out['z0'] = []
    data_out['ydot0'] = []
    for d in data_in['data']:
        # Convert to Halo format
        data_out['x0'].append(float(d[0]))
        data_out['z0'].append(float(d[2]))
        data_out['ydot0'].append(float(d[4]))
        data_out['jc'].append(float(d[6]))
        data_out['period'].append(float(d[7]))

    #sort based on period:
    sort_values = data_out['period']  # what to sort by
    sort_indexes = sorted(range(len(sort_values)), key=lambda i:sort_values[i])
    data_out['x0']      = [data_out['x0'][i] for i in sort_indexes]
    data_out['z0']      = [data_out['z0'][i] for i in sort_indexes]
    data_out['ydot0']   = [data_out['ydot0'][i] for i in sort_indexes]
    data_out['jc']      = [data_out['jc'][i] for i in sort_indexes]
    data_out['period']  = [data_out['period'][i] for i in sort_indexes]

    with open(output_file, 'w') as f:
        f.write(f'//converted from: {input_file}\n')
        f.write(json.dumps(data_out, indent=4))

if __name__ == '__main__':
    """test case"""

    if len(sys.argv) != 3:
        print("Usage: python convert_jpl_database.py <input_file> <output_file>")
        sys.exit(1)
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        if not os.path.exists(input_file):
            print(f"Input file {input_file} does not exist.")
            sys.exit(1)
        convert_jpl_database(input_file, output_file)
        print(f"Converted {input_file} to {output_file}.")