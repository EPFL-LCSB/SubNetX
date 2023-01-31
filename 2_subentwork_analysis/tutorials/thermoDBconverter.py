#!/usr/bin/env python

"""
Thermo Database Converter, Matlab to Python

.. module:: pyTFA
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

:date: June 13th, 2017
:author: pyTFA team

Usage: python thermoDBconverter.py source_database.mat destination.thermodb

You may also import matToPy in your script to convert a thermoDB from matlab to
a python-friendly equivalent of it.


"""


import sys
import os
import numpy
import re
import pickle
import zlib
from scipy.io import loadmat


def matToPy(db):
    """ Convert a Matlab-imported array to a Python Database

    :params numpy.ndarray db: The array of an imported .mat file

    :returns: A dictionnary containing the database data, put in a
    python-friendly format.
    :rtype: dict

    """

    def parse_cue_data(value):
        """ Convert a struct_cue string to a Python dictionnary

        Args:
           value (str): The string to parse

        Returns:
           dict. A dictionnary with the parsed data from the string.

        This is just a little tool used while parsing cue data.
        """

        res = {}
        cues = re.split(r';|\|', value)

        for cue in cues:
            cue = cue.split(':')

            if len(cue) != 2: # Invalid data, for example NULL, NA,...
                continue # Just skip it

            res[cue[0]] = int(cue[1])

        return res

    # Which unit are we using here ?
    units = db['thermo_units'][0,0][0]
    print("Units of the database : " + units)

    ###############
    # METABOLITES #
    ###############

    CompoundDB = db['compound'][0,0]

    # Number of metabolites in the database
    nb_compounds = len(CompoundDB['ID'][0,0])
    print(str(nb_compounds) + " metabolites in database")

    # Initialize the storage for our metabolites
    metabolites={}

    # Fill the data
    for i in range(nb_compounds):
        # Fill the data of the metabolite with the values from Matlab
        met = {}
        met['id'] = CompoundDB['ID'][0,0][i,0][0]
        met['mass_std'] = CompoundDB['Mass_std'][0,0][i,0]
        try:
            met['formula'] = CompoundDB['Formula'][0,0][i,0][0]
        except IndexError:
            met['formula'] = ''
        met['name'] = CompoundDB['metNames'][0,0][i,0][0]
        try:
            met['other_names'] = CompoundDB['AllNames'][0,0][i,0][0].split('|')
        except ValueError:
            met['other_names'] = []
        met['charge_std'] = CompoundDB['Charge_std'][0,0][i,0]

        # Make the cues' data a bit more friendly (as a dict)
        try:
            met['struct_cues'] = parse_cue_data(CompoundDB['struct_cues'][0,0][i,0][0])
        except IndexError:
            met['struct_cues'] = []

        # Continue the import...
        met['nH_std'] = CompoundDB['nH_std'][0,0][i,0]
        met['deltaGf_std'] = CompoundDB['deltaGf_std'][0,0][i,0]
        met['deltaGf_err'] = CompoundDB['deltaGf_err'][0,0][i,0]
        met['error'] = CompoundDB['error'][0,0][i,0][0]
        met['nH_std'] = CompoundDB['nH_std'][0,0][i,0]

        # Special treatment for pkA : convert it to a list of floating values
        met['pKa'] = CompoundDB['pKa'][0,0][i,0][0]
        # BUG: Somehow, the conversion fails for some values and we have sub-arrays...
        if type(met['pKa']) == numpy.ndarray:
            met['pKa'] = met['pKa'][0][0]

        # Create the list...
        met['pKa'] = met['pKa'].split('|')
        # Remove empty values
        met['pKa'] = [x for x in met['pKa'] if x not in [None, '']]
        # If there aren't any values, remove the list
        if met['pKa'][0] == 'NA':
            met['pKa'] = []
        else:
            # Convert the values to float (they're strings)
            met['pKa'] = [float(pKa) for pKa in met['pKa']]

        # Register the metabolite for easier access
        if met['id'] in metabolites:
            raise Exception("Duplicate ID : " + met['id'])
        metabolites[met['id']] = met

    ########
    # CUES #
    ########
    cues = {}

    CueDB = db['cue'][0,0]
    nb_cue = len(CueDB['ID'][0,0])

    for i in range(nb_cue):
        # Import the cue's data
        cue = {}
        cue['id'] = CueDB['ID'][0,0][i,0][0]
        cue['error'] = CueDB['Error'][0,0][i,0][0,0]
        cue['charge'] = CueDB['Charge'][0,0][i,0][0,0]
        cue['energy'] = CueDB['Energy'][0,0][i,0][0,0]
        # If the value is empty, set it to None
        try:
            cue['datfile'] = CueDB['Datfile'][0,0][i,0]
            cue['datfile'] = cue['datfile'][0] if cue['datfile'] else None
        except ValueError:
            cue['datfile'] = None
        cue['formula'] = CueDB['Formula'][0,0][i,0]
        cue['formula'] = cue['formula'][0] if cue['formula'] else None
        try:
            cue['small']   = True if CueDB['Small_Molecule'][0,0][i,0][0,0] == 1 else False
        except IndexError:
            cue['small'] = False
        cue['names']   = CueDB['AllNames'][0,0][i,0][0].split('|')

        # Register the cue in the database
        cues[cue['id']] = cue

    return (units, metabolites, cues)


# If the script was called directly (not imported)
if __name__ == "__main__":
    # Make sure we have all the parameters...
   
                    
    if len(sys.argv) < 3:
        print("Usage : python " + __file__ + " source_database.mat destination.thermodb")
        sys.exit(1)

    # Display a warning if the file already exists
    if os.path.isfile(sys.argv[2]):
        print("WARNING: The destination file already exists, it will be overwritten")

    # Load the Matlab file...
    print("Loading the database in memory...")

    try:
        ReactionDB = loadmat(sys.argv[1])
    except:
        print("Cannot load source file")
        sys.exit(2)

    # Guess the key of the array to the data of the database
    print("Guessing the name of the database...")
    name = None
    for key in ReactionDB:
        if key[:2] == "__": # Name begins with "__" (private attribute)
            continue
        if name == None:
            name = key
        else:
            print("Cannot guess the name of the database ! Found at least 2 candidates : ")
            print(name + " and " + key)
            sys.exit(3)

    if name == None:
        print("Cannot find name, no candidates found")
        sys.exit(3)

    print("Found name : " + name)
    ReactionDB = ReactionDB[name]

    # Convert all the matlab data to a more python-friendly version
    (units, metabolites, cues) = matToPy(ReactionDB)

    # Finally, save the data
    try:
        with open(sys.argv[2], 'wb') as file:
            # Convert it to string with Pickle
            data = pickle.dumps({
                "name":name,
                "units": units,
                "metabolites":metabolites,
                "cues":cues
            }, pickle.HIGHEST_PROTOCOL)

            # Compress it with zlib
            data = zlib.compress(data)

            # And write it to disk
            file.write(data)
    except:
        print('Could not save the data')
        sys.exit(4)

    sys.exit(0)
