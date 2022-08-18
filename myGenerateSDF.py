"""
Created on Sun Aug 14 18:05 2022

@author: SmirnygaTotoshka
"""

#! /usr/bin/env python
# -*- coding: utf-8 -*-
#generate sdf file containing peptide with input sequence
#reference file must be special,ask it to author, they should be with script
import argparse
import json
import os
import traceback
import multiprocessing
import time
import re
import sys
import copy
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem

alphabets = {"protein" : list("ACDEFGHIKLMNPQRSTVWY"),
             "DNA" : list("ATGC"),
             "RNA" : list("AUGC")}

#TODO fix charge aminoacids
def isValidNamedReference(referencePath, alphabet, isCharged):
    files = [os.path.splitext(f)[0] for f in os.listdir(referencePath)]
    charged_exceptions = [i.join("-charge") for i in list("DEHKR")]
    for a in alphabets[alphabet]:
        if isCharged:
            charge = a + "-charge"
            if a not in files and charge not in charged_exceptions:
                return False
        else:
            if a not in files:
                return False
    return True

def loadReference(referencePath):
    molecules = {}
    for file in os.listdir(referencePath):
        name = os.path.splitext(file)[0]
        molecules[name] = open(os.path.join(referencePath, file),"r",encoding="utf-8").readlines()
    return molecules

def isCorrectSequence(sequence, alphabet):
    seq = sequence.strip()
    for i in range(0,len(seq)):
        if seq[i] not in alphabet:
            return False
    return True

def parseReference(sequence, monomers, isCharged):
    atoms = []
    bonds = []
    counts = []
    #get all atoms and bonds of future peptide from reference file
    #reference files have to named after one-letter name amino acid and be "mol" format
    charged_aminoacids = list("DEHKR")
    for s in sequence:
        if isCharged and s in charged_aminoacids:
            s += "-charge"
        ref = monomers[s]
        f = False
        f1 = False
        tmp_atoms = []
        tmp_bonds = []
        for line in ref:
            if "M  V30 BEGIN ATOM" in line.strip():
                f = True
                continue
            if "M  V30 END ATOM" in line.strip():
                f = False
                atoms.append(tmp_atoms)
            if "M  V30 BEGIN BOND" in line.strip():
                f1 = True
                continue
            if "M  V30 END BOND" in line.strip():
                f1 = False
                bonds.append(tmp_bonds)
            if ("M  V30 COUNTS" in line.strip()):
                counts.append(line.strip())
            if (f):
                tmp_atoms.append(line.strip())
            if (f1):
                tmp_bonds.append(line.strip())
    return (atoms,bonds,counts)

#change parameter in line
#returning edited line
#line - input line
#index - paramenter place
#cur_val - new value of parameter
#isCoord - true for ccordinate, write with format 4 digits after point
def setParam(line,index,cur_val,isCoord = False):
    characteristics = re.split("\s+",line)
    out="M  V30 "
    for i in range(2,len(characteristics)):
        if i == index:
            if isCoord:
                out += "{:.4f}".format(cur_val)+ " "
            else:
                out += str(cur_val)+ " "
        else:
            out += characteristics[i] + " "
    return out

#get value of parameter from line
#index must be > 2
def getParam(line,index):
    if index < 2:
        raise IndexError("parameters starts from position 2")
    characteristics = re.split("\s+",line)
    return characteristics[index]

#need delete OH from COOH for build protein skeleton
#+need shift number of atoms after OH. This OH always have number 5 in reference files
def shiftAtoms(atoms):
    shift_ref_atom = copy.deepcopy(atoms)
    for i in range(0,len(shift_ref_atom)):
        flag = False
        j = 0
        while j < len(shift_ref_atom[i]):
            if (int(getParam(shift_ref_atom[i][j],2)) == 5 and not flag and i != len(shift_ref_atom)-1):
                flag = True
                del shift_ref_atom[i][j]
                j -= 1
            elif flag:
                n = int(getParam(shift_ref_atom[i][j],2)) - 1
                shift_ref_atom[i][j] = setParam(shift_ref_atom[i][j],2,n)
            j += 1
    return shift_ref_atom

#According to atom shifting, we shift bonds and bind first C in ref molecule and alpha N in next ref
def shiftBonds(bonds,shift_atoms):
    shifted_bond = copy.deepcopy(bonds)
    shift = 1
    for i in range(0,len(shifted_bond)):
        j = 0
        while j < len(shifted_bond[i]):
            if (int(getParam(shifted_bond[i][j],2)) == 4):
                if (i != len(shifted_bond)-1):
                    shift += len(shift_atoms[i])
                    shifted_bond[i][j] = setParam(shifted_bond[i][j],5,shift)
                j+=1
                continue
            first = int(getParam(shifted_bond[i][j],4))
            second = int(getParam(shifted_bond[i][j],5))
            if(i != len(shifted_bond) - 1):
                if (first >=5):
                    first -= 1
                if (second >=5):
                    second -= 1
            shifted_bond[i][j] = setParam(shifted_bond[i][j],4,first)
            shifted_bond[i][j] = setParam(shifted_bond[i][j],5,second)
            j += 1
    return shifted_bond


#increase x coordinate for nice view and renumber atoms. DEPRECATED
def niceView(shift_atoms):
    ready_atom = copy.deepcopy(shift_atoms)
    shift_num = len(shift_atoms[0])+1
    delta_x = 5
    for i in range(1,len(ready_atom)):
        j = 0
        while j < len(ready_atom[i]):
            ready_atom[i][j] = setParam(ready_atom[i][j],2,shift_num)
            cur_x = float(getParam(ready_atom[i][j],4))
            shift_x = cur_x + delta_x*i
            ready_atom[i][j] = setParam(ready_atom[i][j],4,shift_x,True)
            shift_num += 1
            j += 1
    return ready_atom

#renumber atoms in bonds
def renumberBonds(shift_bonds,ready_atoms):
    ready_bond = copy.deepcopy(shift_bonds)
    shift_num = len(ready_atoms[0])+1
    n = shift_num
    for i in range(1,len(ready_bond)):
        j = 0
        while j < len(ready_bond[i]):
            first = int(getParam(shift_bonds[i][j],4)) + shift_num-1
            second = int(getParam(shift_bonds[i][j],5)) + shift_num-1
            ready_bond[i][j] = setParam(ready_bond[i][j],4,first)
            if not (int(getParam(shift_bonds[i][j],2)) == 4 and i != len(shift_bonds)-1):
                ready_bond[i][j] = setParam(ready_bond[i][j],5,second)
            ready_bond[i][j] = setParam(ready_bond[i][j],2,n)
            n += 1
            j += 1
        shift_num += len(ready_atoms[i])
    return ready_bond

#generate_counts
def generateCounts(ready_atoms,ready_bonds):
    atoms = bonds = 0
    for i in range(0,len(ready_atoms)):
        atoms += len(ready_atoms[i])
    for i in range(0,len(ready_bonds)):
        bonds += len(ready_bonds[i])
    edit_count = "M  V30 COUNTS "
    edit_count += str(atoms) + " "+ str(bonds) +" 0 0 0\n"
    return edit_count

def writeStructure(sequence, reference, isCharged):
    ref_atoms, ref_bonds, ref_counts = parseReference(sequence,reference,isCharged)
    shift_atoms = shiftAtoms(ref_atoms)
    shift_bonds = shiftBonds(ref_bonds, shift_atoms)
    ready_atoms = niceView(shift_atoms)
    ready_bonds = renumberBonds(shift_bonds,ready_atoms)
    edit_count = generateCounts(ready_atoms,ready_bonds)
    structure = ""
    structure += sequence + "\n"
    structure += "Generated by the script of A. Smirnov, Dept. Bioinformatics, Pirogov RNRMU.\n\n"
    structure += "  0  0  0     0  0            999 V3000\n"
    structure += "M  V30 BEGIN CTAB\n"
    structure += edit_count
    structure += "M  V30 BEGIN ATOM\n"
    for mol in ready_atoms:
        for atom in mol:
            structure += atom.strip()+"\n"
    structure += "M  V30 END ATOM\n"
    structure += "M  V30 BEGIN BOND\n"
    for mol in ready_bonds:
        for bond in mol:
            structure += bond.strip()+"\n"
    structure += "M  V30 END BOND\n"
    structure += "M  V30 END CTAB\n"
    structure += "M  END\n"

    molecule = Chem.MolFromMolBlock(structure)
    AllChem.Compute2DCoords(molecule)
    optim_geom = Chem.MolToMolBlock(molecule,forceV3000 = True)
    return optim_geom
    #return structure

def write(proc, partition, output, sequence_column, isCharged, alphabet, reference, output_filename):
    out = os.path.join(output, output_filename + "_thread_" + str(proc) + ".sdf")
    log = os.path.join(output, output_filename + "_thread_" + str(proc) + "_log" + ".txt")
    fail = os.path.join(output, output_filename + "_thread_" + str(proc) + "_failed" + ".txt")
    try:
        with open(out, "w", encoding="utf-8") as o, open(log, "w", encoding="utf-8") as l, open(fail, "w", encoding="utf-8") as f:
            l.write("Start thread #" + str(proc) + " " + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())) + '\n' + "-------------------------------------------------\n")
            for i in partition.index:
                try:
                    l.write("Convert record #" + str(i) + " " + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())) + "\n")
                    if isCorrectSequence(partition.loc[i, sequence_column], alphabets[alphabet]):
                        o.write(writeStructure(partition.loc[i, sequence_column], reference,isCharged))
                        for c in partition.columns:
                            if c != sequence_column:
                                o.write(">  <" + c + ">\n")
                                o.write(str(partition.loc[i, c]) + "\n\n")
                        o.write("$$$$\n")
                        l.write("SUCCESS!\n")
                    else:
                        l.write("FAILED!\n")
                        f.write(str(i) + "\t" + partition.loc[i, sequence_column] + "\n")
                except Exception as er:
                    l.write("Exception on record #" + str(i) + "\t" + partition.loc[i, sequence_column] + "\n")
                    traceback.print_exc()
                    l.write(str(er) + "\n")
                finally:
                    l.write("-------------------------------------------------\n")
    except BaseException as e:
        print("Something went wrong in thread #" + str(proc)+"\n")
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':

    '''
        Path to config file
    '''
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Path to config file.")
    args = parser.parse_args()

    config = args.config

    if os.path.exists(config):
        try:
            with open(config,"r") as cfg:

                '''
                    Parsing arguments
                '''
                parameters = json.loads(''.join(cfg.readlines()))

                input = parameters["input"]
                output = parameters["output"]
                sequence_column = parameters["column"]
                referencePath = parameters["reference"]

                if "charged" in parameters.keys():
                    isCharged = parameters["charged"]
                else:
                    isCharged = False

                if "alphabet" in parameters.keys():
                    alphabet = parameters["alphabet"]
                else:
                    alphabet = "protein"

                if "threads" in parameters.keys():
                    number_threads = parameters["threads"]
                else:
                    number_threads = 1

                if "separator" in parameters.keys():
                    separator = parameters["separator"]
                else:
                    separator = ";"

                if "filename" in parameters.keys():
                    output_filename = parameters["filename"]
                else:
                    output_filename = os.path.splitext(os.path.basename(input))[0]

                '''
                      Validate arguments
                '''

                if not os.path.exists(input) or not os.path.isfile(input):
                    raise BaseException("Input file doesn`t exist or it isn`t file.")
                else:
                    table = pd.read_csv(input, sep=separator, header=0)
                    if sequence_column not in table.columns:
                        raise BaseException("Table doesn`t contain such column.")

                if not os.path.exists(output):
                    raise BaseException("Output directory doesn`t exist.")

                if type(isCharged) != bool:
                    raise BaseException("Which type amino acids residues it should use?")

                if alphabet not in alphabets.keys():
                    raise BaseException("Invalid alphabet. Allow only " + str(list(alphabets.keys())))

                if not os.path.exists(referencePath):
                    raise BaseException("Reference directory doesn`t exist.")
                elif not isValidNamedReference(referencePath, alphabet, isCharged):
                    raise BaseException("Correctly named reference files don`t found.")

                if number_threads < 1 or number_threads > 2 * multiprocessing.cpu_count():
                    raise BaseException("Too many threads. Please use from 1 to " + str(2 * multiprocessing.cpu_count()) + " threads.")

            reference = loadReference(referencePath)

            total = len(table.index)
            size_part = total // number_threads
            size_last_part = size_part + (total - size_part * number_threads)

            # procs - количество ядер
            # calc - количество операций на ядро

            processes = []

            # делим вычисления на количество ядер
            for proc, start in zip(range(number_threads), range(0, total, size_part)):
                if proc == number_threads - 1:
                    partition = table[start:start + size_last_part]
                else:
                    partition = table[start:start + size_part]

                p = multiprocessing.Process(target=write, args=(proc, partition, output, sequence_column, isCharged, alphabet, reference, output_filename))
                processes.append(p)
                p.start()

            # Ждем, пока все ядра
            # завершат свою работу.
            for p in processes:
                p.join()

            sdf = os.path.join(output, output_filename + ".sdf")
            total_log = os.path.join(output, output_filename + "_log" + ".txt")
            total_fail = os.path.join(output, output_filename + "_failed" + ".txt")
            with open(sdf, "w", encoding="utf-8") as final, open(total_log, "w", encoding="utf-8") as log, open(total_fail, "w",encoding="utf-8") as failed:
                for i in range(number_threads):
                    out_t = os.path.join(output, output_filename + "_thread_" + str(i) + ".sdf")
                    log_t = os.path.join(output, output_filename + "_thread_" + str(i) + "_log" + ".txt")
                    fail_t = os.path.join(output, output_filename + "_thread_" + str(i) + "_failed" + ".txt")
                    with open(out_t, "r", encoding="utf-8") as o, open(log_t, "r", encoding="utf-8") as l, open(fail_t, "r", encoding="utf-8") as f:
                        final.write("".join(o.readlines()))
                        log.write("".join(l.readlines()))
                        failed.write("".join(f.readlines()))
            print("Success")
        except BaseException as e:
            print("Something went wrong\n")
            traceback.print_exc()
            sys.exit(1)
        finally:
            end_time = time.time()
            print("--- %s seconds ---" % (end_time - start_time))
    else:
        print("Config doesn`t exist.")
        end_time = time.time()
        print("--- %s seconds ---" % (end_time - start_time))
        sys.exit(1)


