#!/usr/bin/env python3 

from logging.config import dictConfig
import os
from inspect import signature
from random import randrange
import threading
import pymol
import sys
from tkinter.tix import CheckList
import traceback
import csv
import re
from pymol import cmd
import numpy as np


'''
Opens PyMOL with a given pdb-input file and creates
a selection of all the different rgroups, the backbone,
proteins and solvents.

After PyMOL has started up, there is also the possibility to call additional functions from the terminal

When adding new arguments to functions, be sure to update 'variable_dict' in 'read_input' as well. Its keys are used and kwargs.
Also check the expected input types and adjust them in the corresponding lists in 'check_inputs'.

Author: Rene Gall

'''

# List of aminoacids
AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# List of unnatural aminoacids (can be expanded for additional aa)
AA_UNNAT_LIST = ['CYSH', 'HIE', 'HISB']

# Default input_dict
INPUT_DICT_DEFAULT = {"file": None,
                      "ptp_file": None,
                      "representation_molecule": 'sticks', 
                      "representation_solvent": '', 
                      "representation_protein": 'cartoon', 
                      "representation_singleatom": 'nb_sphere', 
                      "peptide_atom_number": [], 
                      "combine_aa_to_protein": True,
                      "core_bound_to_rgroup":  True,
                      "verbose": True}


# Colors for PyMOL
HELIX_COLOR = '0xF00080' # magenta 
SHEET_COLOR = '0xFFFF00' # yellow
LOOP_COLOR = '0xFFFFFF' # white 
MOLECULE_COLOR = [[200, 200, 200], [50, 50, 50]] # grey
ADDITIONAL_SELECTION_COLOR = [[192, 250, 188], [113, 225, 105]] # green



def add_entry_to_dict(dict_add: dict, dict_key, dict_value) -> None:

    '''
    Usage
    -----------
    Add an element to a dictionary. If the key is already present, it will be added to the key as an additional entry,
    otherwise a new key is created and the entry added.

    Parameters
    -----------
    dict_add: dict
        The dictionary the value is added to

    dict_key:
        key to which the value is added.

    dict_value:
        value to be added
    
    '''

    # Check if key exists:
    if dict_key in dict_add.keys():
        dict_add[dict_key].append(dict_value)
        return
    
    # Key does not exist
    # Check if value is a list
    if isinstance(dict_value, list):
        dict_add[dict_key] = [dict_value]
        return
    
    # Value not a list
    dict_add[dict_key] = dict_value
    return




def atom_dict_to_list(atom_dict: dict,
                    skip_solvent: bool = True, skip_protein: bool = True) -> list:

    '''
    Usage
    ---------
    Takes the atom dict and turns it into a list. The list index corresponds to the atom and
    the entry to the molecule name

    Parameters:
    ---------
    atom_dict: dict
        dictionary containing the names and atom indices of the molecules
    
    skip_solvent: bool (Default: True)
        Skips all molecules with "solv_" in name.

    skip_protein: bool (Default: True)
        Skips all molecules with "PROTEIN_" in name.


    Returns
    ---------
    list of atom dictionary
    '''

    atom_molecule_list = [0]

    # Loop through molecules
    for molecule in atom_dict:

        # Check for solvent
        if skip_solvent and "solv_" in molecule:
            continue
            
        # Check for protein
        if skip_protein and "protein_" in molecule:
            continue

        # Loop through atoms
        for atom in atom_dict[molecule]:

            # Write atoms to list
            atom_molecule_list.append(molecule)
    
    return atom_molecule_list



def check_file(filepath: str) -> str:

    '''
    Usage
    -----------
    Checks if the filename is absolute or not and if the file exists. If it exists, the absolute path is returned.

    Parameters
    -----------
    filename: str
        Filename to check
    
    Returns
    -----------
    returns the absolut path to the file 
    '''

    # Check if filepath is absolute, else add current directory
    if not os.path.isabs(filepath):
        print('# No absolute path given. Using the current directory.')
        filepath = os.path.join(os.getcwd(), filepath)
    
    # Check if file exists
    if not os.path.isfile(filepath):
        print(f"# The given file '{filepath}' does not exist. Please check the filename and path and try again. Program terminated.")
        exit(1)

    # File exists, return full absolute path
    return filepath



def check_input(check_dict: dict) -> None:

    '''
    Usage
    -----------
    Checks if all inputs are of the correct type

    Parameters
    -----------
    check_dict: dict
        Dictionarry to check
    
    Returns
    -----------
    raises SystemExit if something is missing or wrong
    '''

    # Check if filename was given
    if 'file' not in check_dict.keys() or check_dict['file'] == '':
        print('# No filename was given. Program terminated.')

    
    # Check, if all entries have correct types. 
    bool_types = ['combine_aa_to_protein', 'core_bound_to_rgroup', 'rgroup_lowest_number', 'verbose']
    str_types = ['file', 'representation_molecule', 'representation_solvent', 'representation_protein', 'representation_singleatom']
    list_types = ['peptide_atom_number', 'coreatom_in_lower_atomindex']

    for entries in check_dict.keys():
        # Check for bools
        if entries in bool_types:
            if not isinstance(check_dict[entries], bool):
                print(f"# Bool expected for '{entries}'. '{type(check_dict[entries])}' gotten. Program terminated.")
                exit(1)

        # Check for strings
        elif entries in str_types:
            if not isinstance(check_dict[entries], str):
                print(f"# String expected for '{entries}'. '{type(check_dict[entries])}' gotten. Program terminated.")
                exit(1)
            

        # Check for lists
        elif entries in list_types:
            if not isinstance(check_dict[entries], list):
                print(f"# List expected for '{entries}'. '{type(check_dict[entries])}' gotten. Program terminated.")
                exit(1)
    
    return




def clean_up_dictionary(clean_up_dict: dict, default_dict: dict=None, clean_up_value = None,
                        set_bools: bool = True) -> None:

    '''
    Usage
    -----------
    Removes all keys with the clean up value and sets them to the default ones. If none is given, delete the elements instead.
    Also sets bools that were strings to real bools and lists to lists.

    Parameters
    -----------
    clean_up_dict: dict
        Dictionarry to clean up
    
    default_dict: dict (default: None)
        dictionary that contains the default values to be set.
    
    clean_up_value (default: None):
        Value which note a key to delete
    
    

    set_bools: bool (Default: True)
        Decides if the strings "True" and "False" should be converted to real bools.


    '''

    for entries in list(clean_up_dict.keys()):
        
        # Check for entries with noted value
        if clean_up_dict[entries] == clean_up_value:
            
            # If default dictionary is given, replace entries, else delete them
            if default_dict is not None:
                clean_up_dict[entries] = default_dict[entries]
                continue
            
            del clean_up_dict[entries]
            
            continue

        # Change bools
        if set_bools:
            change_bool = clean_up_dict[entries]

            # Check if Bool in string
            if change_bool == "True":
                clean_up_dict[entries] = True
                continue
            elif change_bool == "False":
                clean_up_dict[entries] =  False
                continue
    

    return


def combine_elements_to_list(list_to_change:list) -> list:

    '''
    Usage
    -----------
    Takes a list which has separated elements of a list (i. e. list=['x', 'y', '[1', '2', '3]']) and combines them into a list

    Parameters:
    -----------
    list_to_change: list
        list containing a split list
    
    
    Returns:
    -----------
    redone list
    '''


    new_list = [] # New list to be created
    list_in_list = [] # Save list of lists.
    list_found = False
    
    for element in list_to_change:

        # Check for list characters (i. e. '[' and ']')
        if element[0] == '[' and element[-1] == ']':
            list_in_list.append(element[1:-1])
            new_list.append(list_in_list)
            continue

        if element[0] == '[':
            
            # new list found. Add element between bracket and comma
            list_in_list.append(element[1:-1])
            list_found = True
            
            continue

        if element[-1] == ']':

            # End of list
            list_in_list.append(element[:-1])
            list_found = False

            # Add list to new list
            new_list.append(list_in_list)

            # Reset list in list
            list_in_list = []
            continue

        # Add element if currently in list (remove comma at)
        if list_found:
            list_in_list.append(element[:-1])
            continue

        # No list element
        new_list.append(element)

    return new_list


def create_connection_table(connection_table: list,
                    main_atom: int, connected_atoms: int) -> None:

    '''
    Usage
    -----------
    Add atoms to the connection table. The connections are saved in the index positions, the elements between are filled with 0

    Parameters
    -----------
    connection_table: list
        list to add connections to.

    main_atom: int
        is the atom, to which the connected atoms are connected. The number corresponds to the index in the list
    
    conneted_atoms: int
        The connected atoms.
    '''

    # Check, if main atom already in list
    if len(connection_table) > main_atom:
        connection_table[main_atom].append(connected_atoms)
    
    # Main atom not yet in list.
    else:
        # Add 0 to all unconnected indices
        while len(connection_table) != main_atom:
            connection_table.append(0)
        connection_table.append([connected_atoms])

    return



def distribute_colors(min_color: list, max_color: list,
                        length: int) -> list:

    '''
    Usage
    ----------
    Takes in two lists of rgb-values and evenly distributes the colors

    Parameters
    ----------
    min_color: list
        lower bound of the color ranges
    
    max_color:
        upper bound of the color ranges

    length: int
        amount of points to be distributed.


    returns
    -------
    list of list with each rgb-value
    '''  

    # Get distributed colors
    red = np.linspace(min_color[0], max_color[0], length)
    green = np.linspace(min_color[1], max_color[1], length)
    blue = np.linspace(min_color[2], max_color[2], length)

    # combine lists
    color_list = []
    for color in range(0, length):
        color_list.append([red[color], green[color], blue[color]])
    
    return color_list


def expand_connection_table(connection_table: list) -> list:

    '''
    Usage
    ----------
    Fills out all the connections in the connection table with every connection

    Parameters
    ----------
    connection_table: list
        connection_table to be expanded


    returns
    -------
    expanded connection table
    '''   

    # expanded table, first element to 0. List index corresponds to atomindex
    # Get maximal value from connection table to create new list
    max_value = max_value_nested_list(connection_table)
    connection_table_expanded = [[] for _ in range(0, max_value + 1)]

    # Loop over all atoms in the connection table
    for atoms, connections in enumerate(connection_table):
        
        # Skip if connections are 0
        if connections == 0:
            continue

        # Check if atom included in the connected atom
        for connected_atoms in connections:
            if not atoms in connection_table_expanded[connected_atoms]:
                connection_table_expanded[connected_atoms].append(atoms)
        
            # Add connection to current atom
            connection_table_expanded[atoms].append(connected_atoms)
    
    return connection_table_expanded




def get_connections(atom_index: str,
                    connection_table: list, atom_molecule_list: list, rgroup_connected_atoms: list = []) -> list:

    '''
    Usage
    ----------
    Looks through all connections of a given atom and returns all connected atoms except for core atoms

    Parameters
    ----------
    atom_index: str
        index of atom to be checked
    
    connection_table: list
        contains the connection table
    
    atom_molecule_list: list
        List of atom dict, where the atom indices are the list indices and the values the corresponding molecule

    rgroup_connected_atoms: list (Default: [])
        List of r-group connected atoms to be excluded.


    returns
    -------
    List of connections
    '''                

    # Get current molecule and only add connections within this molecule
    current_molecule = atom_molecule_list[atom_index]

    # Start with first connection
    connections = [atom_index]

    # Loop over all connections in list
    for atoms in connections:
        
        # Skip atoms without connections in table
        if connection_table[atoms] == 0:
            continue

        # Loop over connected atoms
        for connected_atoms in connection_table[atoms]:


            # Add to list if not already in and not a core atom and not from different molecule
            if connected_atoms not in connections and connected_atoms not in rgroup_connected_atoms and atom_molecule_list[connected_atoms] == current_molecule:
                connections.append(connected_atoms)
        
    return list(connections)



def get_rgroup_connections(atom_dict: dict, connection_table: list,
                    verbose: bool = True, **kwargs) -> list:

    '''
    Usage
    ----------
    Search for a core structure in a rgroup-core-system


    Parameters
    ----------
    atom_dict: dict
        dictionarry containing all molecule names as keys and corresponding atoms as lists
    
    connection_table: list
        contains the connection table
    
    verbose: bool (default True):
        Tells the user what the program is doing


    returns
    -------
    list of rgroup connected atoms
    '''                
    
    if verbose:
        print("# Looking for r-group connected atoms.")
    
    # Store atoms which are connected to multiple r-groups
    rgroup_connected_atoms = [] 

    # Convert dictionary to list, where the index is the atom and the entry the corresponding molecule
    atom_molecule_list = atom_dict_to_list(atom_dict) 
    
    # Loop through connection table
    for i, connected_atoms_list in enumerate(connection_table):

        # Skip lines with 0
        if connected_atoms_list == 0:
            continue
            

        # Get current molecule and stop if r-groups are finished (no more entries in the shortened atom list)
        if i > len(atom_molecule_list):
            break
        current_molecule = atom_molecule_list[i]

        # Check if connections belong to different molecule and add atom to r-group connecting atoms
        for connections in connected_atoms_list:
            if atom_molecule_list[connections] != current_molecule:
                rgroup_connected_atoms.append(i)
                break
    

    # Check if r-group connecting atoms were found
    if len(rgroup_connected_atoms):
        
        # r-group connecting atoms found. Going on to separate r-groups and cores
        if verbose:
            print(f"# The r-group connected atoms '{rgroup_connected_atoms}' were found. Going on to separate r-groups from cores.")

    # No r-group connecting atoms were found. Aborting separation.
    else:

        if verbose:
            print("# No r-group connected atoms could be found. No separation between r-group and potential cores possible. Continuing with program.")
        
    
    # Return core atoms
    return rgroup_connected_atoms

def read_ptp(filename):
    '''
    Usage
    ----------
    Reads information in the perturbed topology file to figure out the atom
    type codes of all atoms (in each of their perturbed states)

    Parameters
    ----------
    filename: str
        path to the ptp file

    returns
    -------
    atom_type: 
        numpy array of the atom type codes (for all atoms in all perturbed states)    
    dummy_iac:
        atom type code of the dummy atoms
    '''
    read = False

    data = []
    for line in open(filename, 'r'):
        if 'ATOM' in line: read = True
        if 'END' in line: read = False
        if '#' in line or not read: continue

        data.append(line.split())

    num_atoms = int(data[1][0])
    num_states = int(data[1][1])

    # Then go through the data line by line to find the atom types assigned to each state 
    atom_types = []

    for atom in data:
        if len(atom) < num_states *2 + 4:
            continue # don't read lines that are not atoms (potential comments, header, etc.)
        atom_types.append(atom[2:-2:2])

    # convert the data to numpy array
    atom_types = np.array(atom_types, dtype=int)
    # Take the dummy atom type code as the maximum number (gromos convention)
    dummy_iac = np.max(atom_types)
    
    return atom_types, dummy_iac

def assign_to_groups_from_ptp(atom_dict, atom_types, dummy_iac):
    '''
    Usage
    ----------
    Sets all values in the atom_dict properly based on the atom type information
    read from the ptp file    

    Parameters
    ----------
    atom_dict: dict
        dictionarry containing all molecule names as keys and corresponding atoms as lists

    returns
    -------
    atom_dict: dict
        updated atom_dict
    '''
     
    # Number of ligands = number of perturbed states
    num_ligs = len(atom_types[0])
    lig_keys = []
    # Clear all current ligand selections
    for i in range(1, num_ligs+1):
        atom_dict[f'l{i}'] = []
        lig_keys.append(f'l{i}')
    atom_dict['core'] = []
    
    for i, atom in enumerate(atom_types):
        if np.alltrue(atom != dummy_iac):
            atom_dict['core'].append(i+1)
        else:
            for j, lig_key in zip(atom, lig_keys):
                if j != dummy_iac:
                    atom_dict[lig_key].append(i)
    return atom_dict

def max_value_nested_list(check_list: list) -> int:

    '''
    Usage
    -----------
    Finds the highest integer in a nested list containing lists and integers.

    Parameters:
    -----------
    check_list: list
        list of lists to be searched through
    
    '''
    
    max_value = 0
    
    for index, values in enumerate(check_list):
        
        # Check if entry is integer or list
        if isinstance(values, int):
            if values > max_value:
                max_value = values
                continue
            continue
        
        # Get highest element in each list
        if max(values) > max_value:
            max_value = max(values)
            continue

        
    return max_value





def pymol_additional_function(atom_dict: dict, rgroup_connected_atoms_list: list, input_dict: dict) -> None:

    '''
    Usage
    ---------
    Lets the user enter a command from the terminal to execute additional functions in PyMOL

    Parameters:
    ---------
    atom_dict: dict
        dictionary containing the names and atom indices of the molecules


    rgroup_connected_atoms_list: list
        list with all the core atoms

    input_dict: dict
        dictionary with all user defined inputs

    '''
    
    print("# Everything loaded. Additional commands possible. For help, please enter 'help'.\n")
    
    # Dict with all possible functions to call, including their parameters
    possible_functions = {'change_atoms' : 'Input: atoms_to_change (list), selection_they_are_in (str), selection_to_move_to (str, optional)\n\t\t\t\tThe atomindices of the atoms which should be removed from the given selection can be entered.\n\t\t\t\tIf they should be moved to a different selection, the new selection can be added as well (all atoms go to the same new selection).',
                            'change_representation' : '\tInput: selection_name (str), representation (str)\n\t\t\t\tChanges the given selection to the given representation.',
                            'create': '\t\t\tInput: new_name (str), selection_name (str, optional), atom_index (str, optional) \n\t\t\t\tThe selection names and or atom indices may be given to create a new selection. Please enter selections befor atom indices',
                            'hide': "\t\t\tInput: selection_name (list) \n\t\t\t\tWill hide the given selection names. If all should be hidden enter 'all'",
                            'hide_id': '\t\t\tInputs: selection_name (list)\n\t\t\t\tWill hide the atom ID of the given selections',
                            'print_rgroup_connections': '\tPrints a list of found r-group connecting atoms.',
                            'remove_additional_selections': 'Input: None\n\t\t\t\tDeletes all additional sections while keeping changes made to the original selections (like atom swapping). However, a renamed selection will be deleted as well.',
                            'rename': '\t\t\tInputs: name_to_change (str), name_to_change_into (str)\n\t\t\tUsed, if a selection should be renamed. The first name is the name of the current selection, the second name the new one.',
                            'restart': '\t\t\tInput: None\n\t\t\t\the PyMOL-Session will be reset to the inital startup, meaning all changes to selection names, new selections and atom swaps will be unmade.',
                            'show_id':'\t\t\tInputs: selection_name (list)\n\t\t\t\tWill show the atom ID of the given selections.',
                            'show_protein_residues': "\tInputs: list_of_residuenumbers (list), representation (str, optional)\n\t\t\t\tThe given protein residues will be selected and shown (only enter the numbers, e.g. ALA46 -> 46).\n\t\t\t\tUseful for showing certain residues in cartoon mode. Without representation input, the default one or the one given for 'representation_molecule' is chosen.",
                            'show': "\t\tInput: selection_name (list), representation 'rep_xxx' (str, optional)\n\t\t\t\tThe given representations will be made visible with the given selection. This also works to change the current representation. If all selections should be changed, enter 'all'.\n\t\t\t\tThe representation has to be given last and must start with 'rep_'If no representation is given, the user defined representation or the default one will be taken.",
                            'sort': '\t\t\tInput: None\n\t\t\t\tWill sort all current selections in PyMOL in descending order.',
                            'quit': "\t\t\tInput: None\n\t\t\t\tWill terminate the PyMOL-session and end the program. Regular termination by closing PyMOL is also possible."
                            }
    
    # Stay in response mode
    while True:

        try:
            # Get user input
            command_input = input()
            

            # Check for help command
            if 'help' in command_input or command_input == "":
                print("\n#####################################################################################")
                print("Additional commands can now be eyntered.\nThe possible functions are:\n\n")
                for function, info in possible_functions.items():
                    print(function + "\t" + info + '\n')
                    
                print("\nThe individual inputs to the functions (if necessary) can be given without commas, except lists. They should be given as usual with [1, 2, ..., ]")
                print("#####################################################################################\n")
                continue
            
            
            # Extract inputs
            input_list = command_input.rsplit()
            

            # Get length
            input_length = len(input_list)

            # Get command
            command = input_list[0]


            # Check if inputfunction is in dict
            if command not in possible_functions.keys():
                print("# Function '" + input_list[0] + "' not known. Please check and try again or enter 'help' for a list of all functions.")
                continue
                


            # Check which function is called

            # CHANGE_ATOMS
            if command == 'change_atoms':
                
                # Get inputs back into list
                input_list = combine_elements_to_list(input_list)

                # Check number of inputs
                if len(input_list) < 3:
                    print("# Please enter the atom indices to change and the current selection")
                    continue
                
                if len(input_list) == 3:
                    pymol_change_atoms_in_selection([input_list[1]], input_list[2])
                else:
                    pymol_change_atoms_in_selection(input_list[1], input_list[2], input_list[3])


            # CHANGE_REPRESENTATION
            elif command == 'change_representation':
                
                # Check inputs
                if len(input_list) < 3:
                    print("# Please provide a selection and a representation.")
                    continue
                    

                pymol_change_representation(input_list[1], input_list[2])
            
            # CREATE
            elif command == 'create':

                # Get inputs back into list
                input_list = combine_elements_to_list(input_list)
                
                # Check inputs
                if len(input_list) < 3:
                    print("# Please provide at least the new selection and one atom or selection.")
                    continue

                pymol_create_new_selection(input_list[1], input_list[2:], atom_dict.keys())
            
            # HIDE_ID
            elif command == 'hide_id':
                
                # Check inputs
                if len(input_list) < 2:
                    print("# Please provide at least one selection.")
                    continue

                pymol_hide_id(input_list[1:])


            # HIDE
            elif command == 'hide':
                # Check inputs
                if len(input_list) < 2:
                    print("# Please provide at least one selection.")
                    continue
                
                pymol_hide_selections(input_list[1:])
            

            # PRINT_RGROUP_CONNECTIONS
            elif command == 'print_rgroup_connections':
                pymol_print_rgroup_connected_atoms(rgroup_connected_atoms_list)
                

            # REMOVE_ADDITIONAL_SELECTIONS
            elif command == 'remove_additional_selections':
                print("# All additional selections will be deleted and only the original ones and changes made to them will be kept. Continue (y/n)?")
                if input() == 'y':
                    pymol_remove_additional_selections(atom_dict)
            
            # RENAME
            elif command == 'rename':

                # Check if two names are given
                if input_length < 3:
                    print("# Please enter the old and the new name.")
                    continue

                elif input_length > 3:
                    print("# Please only enter one old and one new name.")
                    continue

                pymol_rename_selection(input_list[1], input_list[2])



            # RESTART
            elif command == 'restart':
                print("# All changes made during the PyMOL-Session will be undone. Continue (y/n)?")
                if input() == 'y':
                    pymol_restart(atom_dict, input_dict)



            # SHOW_ID
            elif command == 'show_id':
                
                # Check inputs
                if len(input_list) < 2:
                    print("# Please provide at least one selection.")
                    continue

                pymol_show_id(input_list[1:])


            # SHOW_PROTEIN_RESIDUE
            elif command == 'show_protein_residues':
                
                # Get inputs back into list
                input_list = combine_elements_to_list(input_list)
                
                # Check inputs
                if len(input_list) < 2:
                    print("# Please enter at least a list of residue numbers and optionally the preferred representation.")
                    continue
                
                # Check for optional representation
                if len(input_list) == 3:
                    pymol_show_protein_residues(input_list[1], input_list[2])
                else:
                    pymol_show_protein_residues(input_list[1], input_dict['representation_molecule'])
            
            
            # SHOW
            elif command == 'show':
                
                # Check inputs
                if len(input_list) < 2:
                    print("# Please provide at least one selection.")
                    continue
                
                # Check if representation is also given
                if input_list[-1][0:4] == 'rep_':
                    pymol_show_selections(input_list[1:-1], input_list[-1][4:])
                else:
                    pymol_show_selections(input_list[1:])

            
            # SORT
            elif command == 'sort':
                pymol_sort_selection()


            # QUIT
            elif command == "quit":
                print("# Program is terminating.")
                cmd.quit()
                return
            
            print("# Command executed. Waiting for next command.")

            

            
        except Exception:
            print("\n#####################################################################################")
            print("# Unfortunately something went wrong with the previous command. Please try again. The following message specifies the issue:")
            print("#####################################################################################\n")
            traceback.print_exception(*sys.exc_info())
            print("\n# Waiting for next command.")
    
    return



def pymol_change_atoms_in_selection(atom_index: list, current_selection: str, new_selection: str = None) -> None:
    '''
    Usage
    -----------
    Takes a list of atom indices and deletes them from the current selection. If a new selection is given, they are copied into the new one

    Parameters
    -----------
    atom_index: list
        list of atoms to be changed

    current_selections: str
        name of selection the atoms are currently in

    new_selection: str (Default: None)
        name of new selection, the old atoms should be moved into.

    '''

    selections = cmd.get_names('all')

    # Check if selection is available
    if current_selection not in selections:
        print(f"# The selection '{current_selection}' is not available. Please check and try again.")
        return
    
    # Check if new selection is given and it is available
    if new_selection is not None and new_selection not in selections:
        print(f"# The selection '{current_selection}' is not available. Please check and try again.")
        return

    # Hide old atoms
    command = pymol_command(atom_index, 'id or')
    command = current_selection + ' and (' + command + ')'
    cmd.show_as('', command)

    # Remake selection with new atoms
    command = pymol_command(atom_index, 'id or')
    command = current_selection + ' and not (' + command + ')'

    cmd.select(current_selection, command)


    # Check if second selection is given
    if new_selection is not None:
        

        # Add atom to new selection
        command = pymol_command(atom_index, 'id or')
        command = new_selection + ' or (' + command + ')'
        cmd.select(new_selection, command)
        
        cmd.show_as(INPUT_DICT_DEFAULT['representation_molecule'], new_selection)
        
        # Try to use old color, if not available use new one
        try:

            current_color = cmd.get_color_index(new_selection)
            cmd.color(current_color, new_selection)
        
        except:
            
            # Create new color from range
            new_color = distribute_colors(MOLECULE_COLOR[0], MOLECULE_COLOR[1], 1)
            
            cmd.set_color(new_selection + '_new_color', new_color)
            cmd.color(new_color, new_selection + '_new_color')

        cmd.util.cnc(new_selection)
    
    return


def pymol_change_representation(selection_name:str, representation: str) -> None:
    
    '''
    Usage
    -----------
    Changes the representation of a given selection

    Parameters
    -----------
    selection_name: str
        Name of selection to be changed

    representation: str
        representation to be used.
    
    '''

    # Check if selection exists
    if selection_name not in cmd.get_names('all'):
        print(f"# The selection '{selection_name}' is not available. Please check and try again.")
        return

    # Try to set representation
    try:
        cmd.show_as(representation, selection_name)
        return
    except:
        print(f"# The representation '{representation}' could not be set. Please check and try again.")
        return


def pymol_command(entry_list: list, command_keyword: str) -> str:
    '''
    Usage
    ---------
    Creates a command of the type 'key value key value ...', e.g. 'id 7 or id 8 or id 9...' with 'id or' being the keyword and the numbers being the values in the list

    Parameters:
    ---------
    entry_list
        list which contains the atoms, molecules, ... to be used in pymol

    command_keyword: str
        keyword which will be used in pymol. Consists of either only an identifier (e.g. id, resi, resn, ...) or of an identifier and a selector (e.g. or, and, ...)

    Returns
    ---------
    String with selection command
    '''

    # Split keyword into identifier and selector if possible
    command_list = command_keyword.rsplit()
    identifier = command_list[0]
    selector = ''
    if len(command_list) == 2:
        selector = command_list[1]
    
    
    # Start command with identifier
    command = identifier + ' '

    # Create command
    for index, value in enumerate(entry_list):

        # Add value to command
        command += str(value)

        # Add keyword if not already last value
        if index + 1 < len(entry_list):
            command = command + ' ' + selector + ' ' + identifier + ' '

    
    return command


def pymol_create_new_selection(new_name: str, selection_list:list, original_selections: list) -> None:

    '''
    Usage
    ----------
    Creates a new selection out of a given list with atom indices and selection names.

    Parameters
    ----------
    new_name: str
        Name of the new selection

    selection_list:
        list of the atoms to be selected. Can be selection names or atom indices.

    original_selections: list
        list of automatically created selections
    '''


    # Current selection
    current = cmd.get_names('all')

    # Check if new name already occupied
    while new_name in current:
        print(f"# The name {new_name} is already being used. Please select a different one (or enter 'q' to go back):")
        new_name = input()
        
        # User entered q
        if new_name == 'q':
            return

    # Length
    length = len(selection_list)

    # Create command
    command = ''

    for i, selection in enumerate(selection_list):

        # Check if selection
        if selection in current:
            command += str(selection)
        else:
            command += 'id ' + str(selection)

        # Check if not last entry
        if i < length -1:
            command += ' or '
    
    # Select
    cmd.select(new_name, command)

    # Create representation. Use the default one for molecules
    cmd.show_as(INPUT_DICT_DEFAULT['representation_molecule'], new_name)
    
    # Prevent highlighting of selection
    cmd.disable(new_name)

    # Get all additional selections
    selections = []
    for names in cmd.get_names('all')[1:]:
        if names not in original_selections and 'residue' not in names and 'solv' not in names:
            selections.append(names)
    
            
    # Set colors.
    color_list = distribute_colors(ADDITIONAL_SELECTION_COLOR[0], ADDITIONAL_SELECTION_COLOR[1], len(selections))
    for i, selection in enumerate(selections):
        cmd.set_color(selection + "_color", color_list[i])
        cmd.color(selection + "_color", selection)
        cmd.util.cnc(selection)

    return



def pymol_hide_id(selection_list: list) -> None:

    '''
    Usage
    ----------
    Hides the ID of atoms in the given selection


    Parameters
    ----------
    selection_list: list
        List of selections to hide ID
    
    '''  

    # Loop through inputs
    selection_names = cmd.get_names('all')
    for selection in selection_list:
        if selection in selection_names:
            cmd.hide('labels', selection)
            print(f"# Id for '{selection}' hidden.")

    return  


def pymol_hide_selections(selection_list: list) -> None:

    '''
    Usage
    ----------
    Hides the selections


    Parameters
    ----------
    selection_list: list
        List of selections to hide
    
    '''  

    # Loop through inputs
    selection_names = cmd.get_names('all')
    
    # Check, if all should be hidden
    if 'all' in selection_list:
        for selection in selection_names:
            cmd.show_as('', selection)
            return

    # Not all hidden
    for selection in selection_list:
        if selection in selection_names:
            cmd.show_as('', selection)
            print(f"# '{selection}' hidden.")
            continue
        
        # Selection not in list
        print(f"# '{selection}' is not available. Skipping.")
        continue

    return  


def pymol_make_selections(atom_dict: dict,
                    verbose: bool = True,
                    **kwargs) -> None:
 
    '''
    Usage
    ----------
    Makes selections according to the naming and atom indices in the atom dict

    Parameters
    ----------
    atom_dict:
        dictionary with molecule names as keys and atom indices as list-values.
    
    verbose: bool (Default: True)
        Tell user what program is doing
    
    '''

    # Sort lists in atom dict
    sort_dict(atom_dict)

    if verbose:
        print("# Creating ligand selections")

   
    for molecules in atom_dict:

        # Create selections
        sel_command = pymol_sel_command(atom_dict, molecules)
        cmd.select(molecules, sel_command)

        # Prevents selection from already being highlighted when starting
        cmd.disable(molecules)
    
    
    # Order all selections alphabetically without changing the PyMOL-Object or 'all'
    if verbose:
        print("# Sorting selections.")

    pymol_sort_selection()
    cmd.util.cnc("all")

    return




def pymol_make_visible(selection_names: list,
                       input_dict: dict) -> None:

    '''
    Usage
    ---------
    Displays all the selection with a user defined representation and colors them

    Parameters:
    ---------
    atom_dict: dict
        dictionary containing the names and atom indices of the molecules

    input_dict: dict
        Userdefined inputs

    '''

    if input_dict['verbose']:
        print("# Creating the representations of the selections.")

    # Create representations
    for selection in selection_names:
        
       
        representation, default = pymol_show_as(selection, **input_dict)

        # If representation not possible, change it in input_dict
        if representation is not None:
            input_dict[representation] = default

    
    # Color selections
    if input_dict['verbose']:
        print('# Coloring everything.')
    
    pymol_set_color()


def pymol_print_rgroup_connected_atoms(rgroup_connected_atoms_list: list) -> None:
    
    '''
    Usage
    ---------
    Prints the found r-group connecting atoms of the structure

    Parameters:
    ---------
    rgroup_connected_atoms_list: list
        list with r-group connectint atoms

    '''

    # Check if list is empty
    if len(rgroup_connected_atoms_list):
        if len(rgroup_connected_atoms_list) == 1:
            print(f"# The found r-group connecting atom is at the index '{rgroup_connected_atoms_list[0]}'.")
        else:
            print(f"# The found r-group connecting atoms have the following indices: '{rgroup_connected_atoms_list}'.")
        return
    
    # No atoms found
    print("# No r-group connecing atoms have been identified.")

    return


def pymol_remove_additional_selections(atom_dict: dict) -> None:

    '''
    Usage
    ----------
    Delets all additional selections in PyMol and keeps the ones from the atom dict


    Parameters
    ----------
    atom_dict: dict
        contains information about atoms and selections


    '''  

    # Delete selections
    for selection in cmd.get_names('all')[1:]:
        if selection not in atom_dict.keys() and selection != 'all' and selection != 'current_session':
            cmd.delete(selection)
            continue
        
        if selection in atom_dict.keys() and 'protein' not in selection and 'solv' not in selection:
            cmd.color(selection + '_color' , selection)
            cmd.util.cnc(selection)
            continue
    
    return




def pymol_rename_selection(old_name: str, new_name: str) -> None:
    
    '''
    Usage
    ----------
    Changes the selection name to a new one. It also changes the assigned color (if a random color is present) to the new name


    Parameters
    ----------
    
    old_name: str
        current name of selection
    
    new_name: str
        new name of selection

    '''  

    # Check if old name in selection
    selection_names = cmd.get_names("all")
    if old_name in selection_names:
        
        # Check if new name not already used
        if new_name not in selection_names:
            
            # Rename selection
            cmd.set_name(old_name, new_name)

            # Rename color
            current_color = cmd.get_color_index(old_name)
            current_color = list(cmd.get_color_tuple(current_color))
            cmd.set_color(new_name + '_color', current_color)
            cmd.color(new_name + '_color', new_name)
            cmd.util.cnc(new_name)

            return

        print(f"# The new name '{new_name}' is already used. Please choose a different one.")

        return
    
    # Name not available
    print(f"# The selection name '{old_name}' is not present. Please choose one of the following names:")
    for name in selection_names:
        print(name)
    
    return



def pymol_restart(atom_dict: dict, input_dict: dict) -> None:

    '''
    Usage
    ----------
    Delets all selections in PyMol and recreates them from the atom dict


    Parameters
    ----------
    atom_dict: dict
        contains information about atoms and selections

    input_dict: dict
        dictionary containing user defined inputs

    '''  

    # Delete selections in PyMol
    for selection in cmd.get_names('all')[1:]:
        if selection != 'all' and selection != 'current_session':
            cmd.delete(selection)

    # Make selection
    try:

        pymol_make_selections(atom_dict, **input_dict)
    
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR creating the selections.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        return

    # Make selection visible
    try:

        pymol_make_visible(atom_dict.keys(), input_dict)
    
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR making selections visible.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        return




def pymol_set_color() -> None:

    '''
    Usage
    ---------
    Colors all structures according to the global color definitions

    '''

    # Color helices, sheets, and loops
    cmd.color(HELIX_COLOR, 'ss h')
    cmd.color(SHEET_COLOR, 'ss s')
    cmd.color(LOOP_COLOR, 'ss l')

    # Color selection in range of color
    # Get selection names without first two and without proteins or solvents
    selections = []
    for names in cmd.get_names('all')[1:]:
        if 'protein' in names:
            continue
        if 'solv' in names:
            continue
        selections.append(names)

    # Get colors
    color_list = distribute_colors(MOLECULE_COLOR[0], MOLECULE_COLOR[1], len(selections))

    for i, names in enumerate(selections):

        # Set color in PyMOL
        cmd.set_color(names + '_color', color_list[i])
        cmd.color(names + '_color', names)

    cmd.util.cnc('all')



def pymol_sel_command(atom_dict: dict, current_molecule: str) -> str:

    '''
    Usage
    ---------
    Creates a command to select all atoms in one selection

    Parameters:
    ---------
    atom_dict: dict
        dictionary containing the names and atom indices of the molecules

    current_molecule: str
        current molecule that is selected

    Returns
    ---------
    String with selection command
    '''

    # Check if numbers are consecutive (easier command)
    if len(atom_dict[current_molecule]) == (int(max(atom_dict[current_molecule]) - int(min(atom_dict[current_molecule])) + 1 )):
        command = "id " + str(atom_dict[current_molecule][0]) + "-" + str(atom_dict[current_molecule][-1])
    else:
        command = pymol_command(atom_dict[current_molecule], 'id or')

    
    return command


def pymol_show_as(sel: str,
                representation_molecule: str = INPUT_DICT_DEFAULT['representation_molecule'], 
                representation_solvent: str = INPUT_DICT_DEFAULT['representation_solvent'], 
                representation_protein: str = INPUT_DICT_DEFAULT['representation_protein'], 
                representation_singleatom: str = INPUT_DICT_DEFAULT['representation_singleatom'],
                verbose: bool = True,
                **kwargs) -> str:
    '''
    Usage
    --------
    Takes a selection in pymol and gives it the correct representation.
    If the given representation could not be set, the default one is chosen and returned.


    Parameters
    ----------
    sel: str
        The selections in PyMOL to modify
    
    representation_moleucle: str, optional (default: "sticks")
        representation of molecule in PyMOL
    
    representation_solvent: str, optional (default: "line")
        representation of solvent molecules in PyMOL
    
    representation_protein: str, optional (default: "cartoon")
        represenation of protein in PyMOL

    representation_singleatom: str, optional (default: "nb_sphere")
        representation of single atoms (like Ions) in PyMOL

    verbose: bool (Default: True)
        Tell user what program is doing


    Returns
    ----------
    Default representation together with which representation it was (only if representation from user could not be set)
    
    'None' otherwise
    '''

        
    # Check for moleculetype and try to set the representation from user. If it is not available, use the default arguments
        
    # Solvents
    if 'solv' in sel:
        
        try:
            cmd.show_as(representation_solvent, sel)
            return None, None
            
        except:
            default = signature(pymol_make_visible).parameters['representation_solvent'].default
            print(f"# Preferred representation '{representation_solvent}' for '{sel}' not available. Setting to default value '{default}'.")
            cmd.show_as(default, sel)
            return 'representation_solvent', default
        

    # Proteins
    if 'protein' in sel:
        
        # Show secondary structure
        cmd.dss(sel)
        
        # Color helix and beta sheets
        cmd.color(HELIX_COLOR, "ss h")
        cmd.color(SHEET_COLOR, "ss s")
        cmd.color(LOOP_COLOR, "ss l")
        
        try:
            cmd.show_as(representation_protein, sel)
            return None, None
        
        except:
            default = signature(pymol_make_visible).parameters['representation_protein'].default
            print(f"# Preferred representation '{representation_protein}' for '{sel}' not available. Setting to default value '{default}'.")
            cmd.show_as(default, sel)
            return 'representation_protein', default
  

    # Normal molecules
    # Check for single atoms
    atom_count = cmd.count_atoms(sel)
    if (atom_count == 1):
    
        try:
            cmd.show_as(representation_singleatom, sel)
            return None, None

        except:
            default = signature(pymol_make_visible).parameters['representation_singleatom'].default
            print(f"# Preferred representation '{representation_singleatom}' for '{sel}' not available. Setting to default value '{default}'.")
            cmd.show_as(default, sel)
            return 'representation_singleatom', default
            
    # Multiple atoms

    try:
        cmd.show_as(representation_molecule, sel)
        return None, None

    except:
        default = signature(pymol_make_visible).parameters['representation_molecule'].default
        print(f"# Preferred representation '{representation_molecule}' for '{sel}' not available. Setting to default value '{default}'.")
        cmd.show_as(default, sel)
        return 'representation_molecule', default


    return


def pymol_show_id(selection_list: list) -> None:

    '''
    Usage
    ----------
    Shows the ID of atoms in the given selection


    Parameters
    ----------
    selection_list: list
        List of selections to show ID
    
    '''  

    # Loop through inputs
    selection_names = cmd.get_names('all')[1:]
    for selection in selection_list:
        if selection in selection_names:
            cmd.label(selection, 'ID')
            print(f"# Id for '{selection}' shown.")

    return  



def pymol_show_protein_residues(aa_residues: list, represenation: str) -> None:

    '''
    Usage
    ----------
    Selects the given aminoacid residues and changes their appearance to the given represenation


    Parameters
    ----------
    aa_residues: list
        List of residues to be selected
    
    representation: str
        PyMOL representation
    
    '''   
    
    # Put name together
    command = pymol_command(aa_residues, 'resi or')
    
    # Get selection name
    current_selections = cmd.get_names("all")[1:]
    x = 1
    selection_name = "residues_"
    while selection_name + str(x) in current_selections:
        x += 1
        
    selection_name += str(x)

    # Select and check if protein as selection is available. Otherwise residues with same number
    # from different selections might be selected
    protein_name = ''
    for names in current_selections:
        if 'protein' in names:
            protein_name = names
            break
    if protein_name != '':
        cmd.select(selection_name, protein_name + ' and (' + command + ')')
    else:
        cmd.select(selection_name, command)

    # Color side chains
    pymol.stored.colors = []
    cmd.iterate(selection_name, "stored.colors.append( (chain, resi, name, color))")
    res_colors = {}
    for chain, resi, name, color in pymol.stored.colors:
        if name == 'CA': # c-alpha atom
            res_colors[resi] = cmd.get_color_tuple(color)
    
    for residue in aa_residues:
        residue_color = res_colors[residue]
        cmd.set_color('res_color_' + residue, residue_color)
        cmd.color('res_color_' + residue,'resi ' + residue)
        cmd.util.cnc('resi ' + residue)

    # Try to use selection, if not available use default from 'make_visible'
    try:
        cmd.show_as(represenation, selection_name)
    except:
        cmd.show_as(INPUT_DICT_DEFAULT['representation_molecule'], selection_name)
    cmd.disable(selection_name)


    return


def pymol_show_selections(selection_list: list, representation: str = None) -> None:

    '''
    Usage
    ----------
    Changes the representation of given selections


    Parameters
    ----------
    selection_list: list
        List of residues to be selected
    
    representation: str
        PyMOL representation
    
    '''   

    selection_names = cmd.get_names('all')[1:]

    # Check if all should be changed
    if 'all' in selection_list:
        for selection in selection_names:
            pymol_show_as(selection)
        
        return
    
    # Else change all given selections
    for selection in selection_list:
        if selection in selection_names:
            try:
                cmd.show_as(representation, selection)
            except:
                pymol_show_as(selection)
        else:
            print(f"# Selection '{selection}' could not be found. Continuing with next selection.")
    
    return


def pymol_sort_selection() -> None:

    '''
    Usage
    -----------
    Sorts the selections in PyMOl starting from the 3rd selection. The first one is 'all' and the second one is the PyMOL-Object

    Parameters:
    -----------
    
    '''

    order_string = ''

    # Get names of all selections except first two
    for selection in cmd.get_names('all')[1:]:
        order_string = order_string + ' ' + selection
    
    cmd.order(order_string, 'yes')




def pymol_start(file: str, verbose: bool = True, **kwargs) -> None:

    '''
    Usage
    -----------
    Starts pymol with the correct file

    Parameters:
    -----------
    file: str
        Path to file
    
    verbose: bool (Default: True)
        Tell user what the program is doing
    '''
    

    file = check_file(file)
    if verbose:
        print("# Start loading in file '" + str(file) +"'.")
    
    # Open PyMOL
    pymol.finish_launching(file)    
    cmd.load(file)

    return




def read_argument_file(filename: str) -> list:

    '''
    Usage
    -----------
    Reads the inputs from an argumentfile and returns them as a list

    Parameters
    -----------
    filename: str
        Filename to extract arguments
    
    Returns
    -----------
    list of arguments
    '''

    # Check for file and get absolute path
    filename = check_file(filename)

    # Read in file
    arguments = []
    with open(filename) as inp:
        for line in inp:
            
            # If list in line, split line and combine lists, else just split line
            if '[' in line and ']' in line:
                line = combine_elements_to_list(line.rsplit())
            else:
                line = line.rsplit()             
            
            arguments.extend(line)
            


    print('# Arguments extracted.')
    
    return arguments



def read_pdb_file(file: str,
                peptide_atom_number: list = None, combine_aa_to_protein: bool = True,
                verbose: bool = True,
                **kwargs) -> dict:
    '''
    Usage
    ----------
    Read in a pdb file and separate all neccessary information.

    Parameters
    ----------
    file: str,
        absolute path to file

    peptide_atom_number: list, optional (default: None)
        If peptides are present, the first atoms of each peptide can be given in a list.
        This way they stay as peptides and are not combined into a protein
    
    combine_aa_to_protein: bool, optional (default: True)

    verbose: bool, optional (default: True):
        print out certain steps of what the program ist doing
    
    ---------
    returns:
        filled dictionaries
    ''' 

    # Stores names and atoms of each molecule. Names are used as keys, the atoms are a list
    atom_dict = {}
    
    # Stores connection table of all atoms. The index of the list is the main atom, the atoms that are connected to the
    # main atom are stored at that index as a list.
    connection_table = []

    # Check file
    file = check_file(file)

    # TER in pdb-file denotes new block
    ter_encountered = False
    combined_label = ''

    # Read in file
    try:
        if verbose:
            print(f"# Reading in file '{str(file)}'.")

        with open(file) as inp:
            for line in inp:

                # Format line
                line_list = line.rsplit()

                # Only consider lines with 'ATOM', 'CONECT' or 'TER'

                # ATOM found
                if line_list[0] == 'ATOM':
                    
                    molecule_name = line_list[3].lower()
                    
                    # Check if key already in dict and no new block was encountered
                    if molecule_name in atom_dict.keys() and not ter_encountered:
                        atom_dict[molecule_name].append(int(line_list[1]))
                        continue

                    # Key not yet in dict. Decide on keyword
                    # New molecule block
                    elif ter_encountered:

                        # Check if next block should be combined and get moleculename
                        combined_label = ter_block_label(list(atom_dict.keys()), molecule_name, combine_aa_to_protein, peptide_atom_number, int(line_list[1]))
                        
                        # No name given => Use name in pdb file
                        if combined_label == '':
                            atom_dict[molecule_name] = [int(line_list[1])]
                        
                        # Name given
                        else:
                            atom_dict[combined_label] = [int(line_list[1])]

                        # Set TER back to False
                        ter_encountered = False
                        continue

                
                    # No new block, but atom is in peptide list
                    elif len(peptide_atom_number) > 0 and (int(line_list[1]) in peptide_atom_number and (molecule_name in AA_LIST or molecule_name in AA_UNNAT_LIST)):
                        combined_label = ter_block_label(list(atom_dict.keys()), molecule_name, combine_aa_to_protein, peptide_atom_number, int(line_list[1]))
                        atom_dict[combined_label] = [int(line_list[1])]
                        continue
                    
                    #  combined block already found. Add atoms to given name
                    elif combined_label != '':
                        atom_dict[combined_label].append(int(line_list[1]))
                        continue
                    
                    # name of molecule changed between blocks. Use new name
                    else:
                        atom_dict[molecule_name] = [int(line_list[1])]
                        continue
                            

                elif line_list[0] == 'TER':
                    ter_encountered = True
                    combined_label = ''
                    continue

                
                # CONNECTION found
                elif line_list[0] == 'CONECT':

                    # Add connected atom to list at position of main atom
                    create_connection_table(connection_table, int(line_list[1]), int(line_list[2]))
                    continue

                # End of Frame reached (only first frame considered)
                elif line_list[0] == 'ENDMDL':
                    break
                    
        if verbose:
            print("# All inputs read.")


        # If only one protein, solvent or peptide. encountered, remove added '_1' from name
        remove_label_ending(atom_dict, ['protein_', 'solv_','peptide_'], ['protein', 'solv', 'peptide'])




    except Exception:
        print("#####################################################################################")
        print("\t\tERROR while reading in file")
        print("#####################################################################################")

        traceback.print_exception(*sys.exc_info())
        return -1

    return atom_dict, connection_table




def read_user_inputs(arg_list: list) -> dict:

    '''
    Usage
    -----------
    Read in inputs from command line and store them in the correct variables
    in an input_dictionarry

    Parameters:
    -----------
    arg_list: list
        all arguments from the command line

    Returns
    -----------
    Dictionarry with all inputs or 0 if 'help' was called
    '''

    # Create input_dictionarry with all possible options. The keys are the same as the function inputs
    # in order to pass the dictionarry as a kwarg
    variable_dict = {"file": "\t\t\t\t'str': Enter a filename or absolute path of the pdb-file. If no absolute path is given, the current directory is chosen.", 
                     "ptp_file":"\t'str': Enter the filename of the perturbed topology file",
                     "representation_molecule": "\t'str': Specify the preferred representation of the molecule. Default is 'sticks'.", 
                     "representation_solvent": "\t'str': Specify the preferred representation of the solvent. Default is '' (solvent hidden).", 
                     "representation_protein": "\t'str': Specify the preferred representation of the protein. Default is 'cartoon'.", 
                     "representation_singleatom": "\t'str': Specify the preferred representation of single atoms. Default is 'nb_sphere'.", 
                     "peptide_atom_number": "\t\t'list': Used if molecules in the pdb-file are peptides and should not be combined into a protein. The first atom of every peptide can be given in a list. Default is 'None' .", 
                     "combine_aa_to_protein": "\t\t'bool': Used, if aminoacids are present that should belong to a protein and be combined. Default is 'True'.",
                     "core_bound_to_rgroup": "\t\t'bool': Used, if multiple rgroups share a common core structure. The core is then labelled separately. If multiple corestructures are present, the program tries to assign all cores correctly. Default is 'True'.",
                     "verbose": "\t\t\t'bool': Used, if feedback of the program is wished as to what it is doing at the moment. Default is 'True'."}
    key_list = list(variable_dict.keys())

    
    # Check, if 'help' is in arg_list and display possible option
    if 'help' in arg_list:
        print("\n#####################################################################################")
        print("\nThis program is used to auto select all molecules and rgroups with core structures in a given pdb-file.\nThe possible inputs are:\n\n")
        for variable, info in variable_dict.items():
            print(f"\t@{variable} \t {info}\n")
        print("\nTo set the inputs, either enter '@input value' to set specific inputs, or enter them in order without '@'.")
        print("\nThe inputs may also be read in with a separate file with '@f' followed by the argument file.\n")
        print("\nWhen using the 'core_bound_to_rgroup'-option (or by not giving an argument for it), the program tries to separate a rgroup core and the corresponding rgroup.\nWhen more than one core-structur is present, i.e. some rgroups on a main core, and then additional rgroups on a sub-core, the result may look a bit off.\n In this case, additional functions can be called after PyMOL has loaded to put individual atom to the preferred selection.\n")
        print("\nAfter the program has fully loaded, additional commands can be called from the terminal to be executed in PyMOL. These can again be seen when entering help.")
        print("#####################################################################################\n")
        
        # Exit program
        return 0
    
    
    # Clear out help-information while keeping the keys.
    for item in variable_dict:
        variable_dict[item] = None


    
    # Check if inputs are given in order (no '@')
    if '@' not in str(arg_list):
        
        # Check if too many arguments given
        if len(arg_list) > len(key_list):
            print('# Too many arguments are given. Program terminated.')
            exit(1)


        # Fill in dict in order of keys
        for i, argument in enumerate(arg_list):
            variable_dict[key_list[i]] = argument
    
    
    # Specific entries are given
    else:
        
        # Check, if an argument file is given.
        if '@f' in arg_list:
            
            # Argumentfile not as only argument given. Exit program
            if arg_list[0] != '@f' or len(arg_list) != 2: 
                print('# "@f" found. Please put "@f" as the first argument and only give the corresponding filename.')
                exit(1)

            print('# Argument file found. Extracting arguments.')

            # Get arguments from file into list input
            arg_list = read_argument_file(arg_list[1])
    
        
        # Set all arguments in variable_dict
        length_list = len(arg_list)

        skip_value = False
        for i, elements in enumerate(arg_list):
            
            # Skip first element after key
            if skip_value:
                skip_value = False
                continue

            # Check if key is in variable_dict (Subtract '@' from argument)
            if elements[0] == '@' and elements[1:] in key_list:

                # Check if an element is given and add it to the dict. Then skip input value (Assumption: always one key ('@') and one value)
                if i + 1 < length_list and '@' not in str(arg_list[i + 1]):
                    variable_dict[elements[1:]] = arg_list[i + 1]
                    skip_value = True
                
                else:

                    # No input was given. Set default value for argument.
                    print(f"# No input was given for '{elements}'. Using the default value '{INPUT_DICT_DEFAULT[elements[1:]]}' and continuing.")
                    variable_dict[elements[1:]] = INPUT_DICT_DEFAULT[elements[1:]]
                    
                    # Do not skip over next element. Will be an argument
                    skip_value = False

            # Key in arguments, but not in dict
            elif elements[:1] == '@':
                print(f"# Argument '{elements}' unknown. Please check inputs and try again. Program terminated.")
                exit(1)
            
            # More than one argument given
            else:
                print(f"# More than one argument given for '{arg_list[i -2]}'. Please check inputs and try again.")
                exit(1)

    
    # Replace all empty elements with the default elements and change string bools to bools
    clean_up_dictionary(variable_dict, INPUT_DICT_DEFAULT)

    
    
    # Check if filename is given and all entries have the correct type
    if check_input(variable_dict):
        print('# Error with argument types. Please check printed messages and try again. Program terminated.')
        exit(1)

    return variable_dict



def remove_label_ending(label_dict: dict,
                        label_beginning: list, replace_label: list) -> None:


    '''
    Usage
    ----------
    Checks if a given labelname appears only once in a dict and renames the key


    Parameters
    ----------
    label_dict: dict
        dictionarry containing the labels as keys.
    
    label_beginning: list
        Defines how the labelname starts. Multiple entries at once possible

    replace_label: list
        replaces the key with the given label. The first entry in this list will replace the first entry from label_beginning


    '''     

    # Create dict to store entries. Keys are the label beginning, the first entry in the list is the found name,
    # the second one the name it should be replaced with and the third is used if multiple instances are found.
    remove_dict = {}
    for beginning, replace in list(zip(label_beginning, replace_label)):
        remove_dict[beginning] = ['', replace, False]


    for names in label_dict.keys():
        
        # Check for label beginning
        remove = [beginning for beginning in label_beginning if names.startswith(beginning)]

        # Check if no instance found
        if not remove:
            continue
        
        remove = remove[0]

        # Check if instances were already found
        if remove_dict[remove][2]:
            continue

        # Check if one instance was found
        if remove_dict[remove][0] != '':
            remove_dict[remove][2] = True
            continue
        
        # Set name to be deleted
        remove_dict[remove][0] = names

    
    # Rename labels with just one occurence
    for names in remove_dict.values():
        
        if names[0] != '' and not names[2] and names[0] in label_dict.keys():
            label_dict[names[1]] = label_dict.pop(names[0])


    return



def separate_rgroups_and_cores(atom_dict: dict, connection_table: list, rgroup_connecting_atoms_list: list, 
                            verbose: bool = True,
                            **kwarg) -> None:
    
    '''
    Usage
    ----------
    Separates the core atoms from the r-groups


    Parameters
    ----------
    atom_dict: dict
        dictionarry containing all molecule names as keys and corresponding atoms as lists
    
    connection_table: list
        contains the connection table

    rgroup_connecting_atoms_list: list
        contains the rgroup connecting atoms

    
    verbose: bool (default True):
        Tells the user what the program is doing


    '''     

    # Get atom dictionary as list to easily search through
    atom_molecule_list = atom_dict_to_list(atom_dict)   

    # Get full connection table
    connection_table_expanded = expand_connection_table(connection_table)

    # Dict to store selections per r-group connecting atom
    selection_dict = {} 

    # Get first molecule name
    current_molecule = ''

    if verbose:
        print('# Checking the connections of found r-group connected atoms.')

    # Loop through found core atoms
    for rgroup_connection in rgroup_connecting_atoms_list:

        # Check if next r-group connecting atom is already in same molecule and skip
        # because the whole fragmentations is already done.
        if atom_molecule_list[rgroup_connection] == current_molecule:
            continue
        else:
            current_molecule = atom_molecule_list[rgroup_connection]
        
        # Store used atoms
        used_atoms = [[0]]

        # Loop throug atoms in selection
        for atom in atom_dict[current_molecule]:

            # Atom already added or it is a r-group connecting atom
            if  atom in used_atoms[0] or atom in rgroup_connecting_atoms_list:
                continue
            
            # Get all connections of current atom excluding r-group connecting atoms
            selection_connection_list = get_connections(atom, connection_table_expanded, atom_molecule_list, rgroup_connecting_atoms_list)
            used_atoms[0].extend(selection_connection_list)
            
            # Add selection to dict
            add_entry_to_dict(selection_dict, rgroup_connection, selection_connection_list)
    

    if verbose:
        print('# Connections have been found. Setting new selection names.')

    # Decide on r-groups and cores    
    for rgroup_connection, selected_atoms in selection_dict.items():

        # Check if list has more than one entry
        if len(selected_atoms) == 1:
            continue

        # Delete entries in atom dict with r-group connected atoms and and r-groups
        del atom_dict[atom_molecule_list[rgroup_connection]]
        
        # Sort selected atoms by size
        selected_atoms.sort(key = len)

        rgroup_name = atom_molecule_list[rgroup_connection]

        # biggest list probably core.
        core_name = 'core_'+ rgroup_name
        atom_dict[core_name] = selected_atoms[-1]

        
        # Add r-group connected atoms to selection
        add_rgroup_connected_atoms = []
        add_single_atoms_to_rgroup_atoms = []
        
        if verbose:
            print('# All names set. Adding r-group connected atoms to new selections.')

        for rgroup_connection in rgroup_connecting_atoms_list:
            
            # Check if r-group connected atom belongs to same r-group
            if atom_molecule_list[rgroup_connection] == rgroup_name:
                add_rgroup_connected_atoms.append(rgroup_connection)
                add_single_atoms_to_rgroup_atoms.append(rgroup_connection)
        
        # Check if r-group connected atom have to be added, else skip
        if len(add_rgroup_connected_atoms):
            for rgroup_connection in add_rgroup_connected_atoms:
                
                # Add r-group connected if they are only bound to these.
                connections = connection_table_expanded[rgroup_connection]
                for connected_atom in connections:
                    if len(connection_table_expanded[connected_atom]) == 1 and connection_table_expanded[connected_atom][0] == rgroup_connection and atom_molecule_list[connected_atom] == rgroup_name:
                        add_single_atoms_to_rgroup_atoms.append(connected_atom)
                        for i, selection in enumerate(selected_atoms):
                            if len(selection) == 1 and connected_atom in selection:
                                del selected_atoms[i]
                                break
                            selection.remove(connected_atom)
                            break
        
        # Add found atoms to core and delete from subcores
        atom_dict[core_name].extend(add_single_atoms_to_rgroup_atoms)

        # Add remaining subcores to atom dict
        if len(selected_atoms) > 2:
            for subcore in range(1, len(selected_atoms)):
                subcore_name = 'sub_core_' + str(subcore + 1) + rgroup_name
                atom_dict[subcore_name] = selected_atoms[subcore]
        elif len(selected_atoms) == 2:
            atom_dict[rgroup_name] = selected_atoms[0]

    # Remove label ending from "CORE_" and "SUB_CORE_" if only one has been found
    remove_label_ending(atom_dict, ['core_', 'sub_core_'], ['core', 'sub_core'])
    
    
    return
    


def sort_dict(dict_to_sort: dict) -> None:
    
    '''
    Usage
    ----------
    Sorts the elements of a dictionary of lists in ascending order


    Parameters
    ----------
    dict_to_sort: dict
        dictionary of lists to be sorted
    
    '''        
    
    for sort_list in dict_to_sort.values():
        sort_list.sort()
    return




def ter_block_label(molecule_keys: list, molecule_name: str, 
                    combine_aa_to_protein: bool, peptide_atom_number: list,
                    atom_number: int) -> str:

    '''
    Usage
    -----------
    Looks at the input line and creates a new name if a protein, peptide, solvent, ... is encountered.
    If none of these are present, an empty string is returned.

    Parameters:
    -----------
    molecule_keys: list
        The molecule names currently saved in the atom dictionary.

    molecule_name: str
        label according to pdb file

    combine_aa_to_protein: bool
        if True, the encountered aminoacids will be put toghether as one protein
    
    peptide_atom_number: list
        Used, if some molecules are peptides. In this case, the first atom of the respective peptide is given.

    atom_number: int
        Current atom index

    Returns
    -----------
    returns label for molecule
    '''

    combined_label = ''
    # Check if aminoacid encountered
    if molecule_name.upper() in AA_LIST or molecule_name.upper() in AA_UNNAT_LIST:

        # Check if peptide
        if len(peptide_atom_number) > 0 and atom_number in peptide_atom_number:
            combined_label = 'peptide_'
        
        # Check if protein
        elif combine_aa_to_protein:
            combined_label = 'protein_'
    
    # Check if solvent
    elif molecule_name == 'solv':
        combined_label = 'solv_'
    
    # Check if name is empty
    if combined_label != '':
        x = 1
        # Check if name is in dictionary, else increase numbering
        while combined_label + str(x) in molecule_keys:
            x += 1
        
        combined_label = combined_label + str(x)
    
    # New block, but same label as already given
    elif molecule_name in molecule_keys:
        x = 2
        while molecule_name + '_' + str(x) in molecule_keys:
            x += 1
        combined_label = molecule_name + '_' + str(x)
    
    return combined_label






# Main program
def main():


    # Input dictionarry used to store all parameters for the program
    input_dict = {}
    
    # Check if  program has enough inputs
    
    if (len(sys.argv) < 2):
        print("#####################################################################################")
        print("\tNot enough inputs.\n\tPlease enter a filename or enter 'help' for more information.")
        print("\tStopping program.")
        print("#####################################################################################")
        exit(0)
    

    # Read in the inputs and create the input_dictionarry
    try:
        # Remove first entry from input (only program information)
        arguments = sys.argv[1:]

        input_dict = read_user_inputs(arguments)

        # If 'help' was in arguments, input_dict is 0 and program terminates.
        if not input_dict:
            exit(0)
        
        
    # Unknown Error
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR when getting input variables.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        exit(1)


    # Start PyMOL and load in file with separate thread
    try:
        
        pymol_thread = threading.Thread(target=pymol_start, kwargs=input_dict)
            
        pymol_thread.start()
        
        # Unknown Error
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR when starting up PyMOL.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        exit(1)
        

    # Read in input file and get atom numbering and connection table
    try:
        atom_dict, connection_table = read_pdb_file(**input_dict)
        
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR when reading in input file.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        sys.exit(1)

    if input_dict['ptp_file'] is not None:
        atom_types, dummy_iac = read_ptp(input_dict['ptp_file'])
        atom_dict = assign_to_groups_from_ptp(atom_dict, atom_types, dummy_iac)

    else:
        try:# Find r-group connecting atoms and separate  into r-group and cores
            rgroup_connecting_atoms_list = []

            # 'Check for corestructure' is on per default
            if input_dict['core_bound_to_rgroup'] == True :
                rgroup_connecting_atoms_list = get_rgroup_connections(atom_dict, connection_table, **input_dict)    
                
                # Continue with separation if r-group connecting atoms were found, else skip
                if len(rgroup_connecting_atoms_list):
                    separate_rgroups_and_cores(atom_dict, connection_table, rgroup_connecting_atoms_list, **input_dict)
                
        except Exception:
            print("#####################################################################################")
            print("\t\tERROR searching for core structures and separating them.")
            print("#####################################################################################")
            traceback.print_exception(*sys.exc_info())
            sys.exit(-1)

    # Wait for PyMOL to finish loading
    print("# Waiting for PyMOL to finish loading. This may take a few moments.")
    pymol_thread.join()

    # Make selection
    try:

        pymol_make_selections(atom_dict, **input_dict)
    
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR creating the selections.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        sys.exit(-1)

    # Make selection visible
    try:

        pymol_make_visible(atom_dict.keys(), input_dict)
        
    except Exception:
        print("#####################################################################################")
        print("\t\tERROR making selections visible.")
        print("#####################################################################################")
        traceback.print_exception(*sys.exc_info())
        sys.exit(-1)


    # Use new thread to receive additional inputs
    additional_commands = threading.Thread(target=pymol_additional_function, args=[atom_dict, rgroup_connecting_atoms_list, input_dict], )
            
    additional_commands.start()





if __name__ == "__main__":
    main()
