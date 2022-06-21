# HybridTopologyViewer

----------

## OVERVIEW

This program can be used to help view '.pdb' topology files in PyMOL. It is mainly focussed on viewing hybrid topologies, however it can also be used for regular '.pdb'-files.

The input to this program is a '.pdb'-file where in the 2nd column the atom index and in the 4th column the group label is written for 'ATOM' (in the first column) and the connection table is given with 'CONECT' (in the first column).

The output of the program will be in PyMOL. There all groups with the same label in the pdb-file will be combined into a PyMOL-Object with the label as the name.
Per default every protein will be shown as cartoon with alpha-helices colored in magenta, the beta-sheets in yellow and the loops in white. The solvents will be hidden. All other molecules will have a range of grey colors and be shown in sticks.
All of these default values can easily be changed by adjusting the global variables in the program.

----------

## STARTING THE PROGRAM

When starting the program from the command line, only the filename has to be given:

```
HybridTopologyViewer.py filename.pdb
```

This can be either done with an absolute or relative path. If the path is relative, the current directory is chosen to look for.

There are also additional commands to be given. These can be added with '@command argument'. In this case, the filename has to be given with '@file filename.pdb'. The commands are as follows:

```
@file: 				'str' 
		                Enter a filename or absolute path of the pdb-file. If no absolute path is given, the current directory is chosen. 

@representation_molecule: 	'str' 
		                Specify the preferred representation of the molecule. Default is 'sticks'.

@representation_solvent: 	'str' 
		                Specify the preferred representation of the solvent. Default is '' (solvent hidden). 

@representation_protein: 	'str'
		                Specify the preferred representation of the protein. Default is 'cartoon'.

@representation_singleatom: 	'str'
		                Specify the preferred representation of single atoms. Default is 'nb_sphere'. 

@peptide_atom_number: 		'list'
		                Used if molecules in the pdb-file are peptides and should not be combined into a protein. The first atom of every peptide can be given in a list. Default is 'None' .

@combine_aa_to_protein: 	'bool'
		                Used, if aminoacids are present that should belong to a protein and be combined. Default is 'True'.

@core_bound_to_rgroup: 		'bool'
		                Used, if multiple rgroups share a common core structure. The core is then labelled separately. If multiple corestructures are present, the program tries to assign all cores correctly. Default is 'True'.

@verbose:                       'bool'
		                Used, if feedback of the program is wished as to what it is doing at the moment. Default is 'True'.
```

A possible input to load in a file, hide the protein and cancel any feedback from the program could therefore be:

```
HybridTopologyViewer.py @file filename.pdb @representation_protein '' @verbose False
```

If multiple keywords have to be given for multiple files then these can be saved in a separate textfile which can then be given. Note that the filename has to be included within this textfile and no other arguments can be given:

```
HybridToplogyViewer.py @f textfile.arg
```

These keywords can also be viewed when typing 'help' after the program name:

```
HybridTopologyViewer.py help
```

If the keywords need to be adjusted on a regular base it is also possible to adjust them in the global dictionary names 'INPUT_DICT_DEFAULT'.

----------

## ADDITIONAL COMMANDS

After everything is loaded into PyMOL the program stays responsive for further commands. These can be entered into the terminal the program was started with (NOT into the PyMOL-terminal). The possible functions are:

```
change_atoms: 			Input: atoms_to_change (list), selection_they_are_in (str), selection_to_move_to (str, optional)
			        The atomindices of the atoms which should be removed from the given selection can be entered. If they should be moved to a different selection, the new selection can be added as well (all atoms go to the same new selection).

change_representation: 		Input: selection_name (str), representation (str)
			        Changes the given selection to the given representation.

create:				Input: new_name (str), selection_name (str, optional), atom_index (str, optional)
			        The selection names and or atom indices may be given to create a new selection. Please enter selections befor atom indices.

hide:				Input: selection_name (list)
			        Will hide the given selection names. If all should be hidden enter 'all'.

hide_id:			Inputs: selection_name (list)
			        Will hide the atom ID of the given selections.

print_rgroup_connections:	Inputs: None
			        Prints a list of found r-group connecting atoms.

remove_additional_selections:	Input: None
			        Deletes all additional sections while keeping changes made to the original selections (like atom swapping). However, a renamed selection will be deleted as well.

rename:				Inputs: name_to_change (str), name_to_change_into (str)
			        Used, if a selection should be renamed. The first name is the name of the current selection, the second name the new one.

restart:			Input: None
			        The PyMOL-Session will be reset to the inital startup, meaning all changes to selection names, new selections and atom swaps will be unmade.

show_id:			Inputs: selection_name (list)
			        Will show the atom ID of the given selections.
show_protein_residues:		Inputs: list_of_residuenumbers (list), representation (str, optional)
			        The given protein residues will be selected and shown (only enter the numbers, e.g. ALA46 -> 46).Useful for showing certain residues in cartoon mode. Without representation input, the default one or the one given for 'representation_molecule' is chosen.

show:				Input: selection_name (list), representation 'rep_xxx' (str, optional)
			        The given representations will be made visible with the given selection. This also works to change the current representation. If all selections should be changed, enter 'all'. The representation has to be given last and must start with 'rep_'If no representation is given, the user defined representation or the default one will be taken.

sort:				Input: None
			        Will sort all current selections in PyMOL in descending order.

quit:				Input: None
			        Will terminate the PyMOL-session and end the program. Regular termination by closing PyMOL is also possible.
```

A possible way to hide everything in PyMOL and to only show the core-structure and the first to r-groups with two commands could be:

```
hide all
show core r1 r2
```

where 'core' is the name of the core-object and 'r1' and 'r2' are the names of rgroup 1 and 2 respectively.

These commands can also be shown by typing 'help' into the terminal.
----------

WORKING OF PROGRAM

After reading everything in from the command line the program starts to read in the .pdb-file and in a first round groups together everything according to the labels in the file. This way, all solvents and proteins are already being separated and only the core-structure and the r-groups of the hybrid toplogy remain.
Also most of the r-groups will be separated in this step as well. In the last step the program tries to determine which part of the groups belongs to the core-structure. For this it first searches for an atom that is bound to multiple r-groups. Then it will split the selection at this point and regroups everything on both sides into separate groups.

After all the splitting and regrouping is done, the program starts to color in everything in PyMOL according to the default values and stays in response mode for further commands.
