
# # Function Declaration # #
# the following is heavily adapted from Vincent F's DTERMen analysis code
# https://github.com/KeatingLab/dTERMen_design/blob/master/jupyter-notebook/1_Concat_scoring_method.ipynb
# NOTE: this can only score a single chain - if you have a dTERMen etab for multiple disjoint chains,
#       you will need to modify this


# ## Import Modules # #
# the following are packages that must be installed if you want to run this script
import glob
import regex as re
import argparse

# ## Dictionary defined for converting three letter amino acid code to one letter ## #

three_to_one_letter = {"ARG":"R",
                       "HIS":"H",
                       "LYS":"K",
                       "ASP":"D",
                       "GLU":"E",
                       "SER":"S",
                       "THR":"T",
                       "ASN":"N",
                       "GLN":"Q",
                       "CYS":"C",
                       "GLY":"G",
                       "PRO":"P",
                       "ALA":"A",
                       "VAL":"V",
                       "ILE":"I",
                       "LEU":"L",
                       "MET":"M",
                       "PHE":"F",
                       "TYR":"Y",
                       "TRP":"W"}

# ## Function/ Definitions ## #

# Reads energy table file and generate dictionaries that store the self and pair energy terms
# The dictionaries are then used by scoreSequence() to calculate the pseudoenergy of the sequence

# NOTE: this is adapted for parsing MST-dTERMen energy tables, which have a slightly updated format
#       self energies should be written as CHAIN,RESID RESNAME ENERGY (e.g. B,601 ALA 3.03801)
#       pair energies should be written as CHAIN1,RESID1 CHAIN2,RESID2 RESNAME1 RESNAME2 ENERGY

def loadDTERMenEnergyTable(etab_path):
    #import all lines from the energy table
    with open(etab_path, "r") as file:
        lines = [line.rstrip('\n') for line in file] #reads in all of the lines from the etab

    #initialize the structures to store designable positions and energy parameters
    self_table = dict() #stores the self energy parameters
    pair_table = dict() #stores the pair energy parameters

    #determine the indexing scheme
    first_index = int(re.split(",",re.split("\s+",lines[0])[0])[1])
    #NOTE: this is sort of a hack that takes advantage of the fact that
    #       dTERMen starts by outputting the lowest index. This value is used
    #       to readjust the following values - such that everything is zero indexed

    for l in lines:
        sp = re.split("\s+",l) #parses the line into a list of values

        if len(sp) == 3: #self energy term
            chain,index = re.split(",",sp[0]) #chain is not currently stored
            index = int(index) - first_index
            amino_acid = three_to_one_letter[sp[1]]
            self_energy = float(sp[2])

            if index not in self_table:
                self_table[index] = dict()
            self_table[index][amino_acid] = self_energy
        elif len(sp) == 5: #pair energy term
            chain,index_1 = re.split(",",sp[0])
            chain,index_2 = re.split(",",sp[1])
            index_1 = int(index_1)-first_index
            index_2 = int(index_2)-first_index
            amino_acid_1 = three_to_one_letter[sp[2]]
            amino_acid_2 = three_to_one_letter[sp[3]]
            pair_energy = float(sp[4])

            if int(index_1) >= int(index_2):
                print("Warning: etab formatting does not match expectation")

            pair_index = (index_1,index_2)
            pair_amino_acid = (amino_acid_1,amino_acid_2)

            if pair_index not in pair_table:
                pair_table[pair_index] = dict()

            pair_table[pair_index][pair_amino_acid] = pair_energy
        else:
            print("line does not match the expected .etab formatting")
    return (self_table,pair_table)

#Score sequence using energy function
#NOTE: I forgot that dTERMen doesn't include the zero-valued energies
#       easy fix is just try/catch addition
def scoreSequence(sequence,self_table,pair_table):
    if len(sequence) != len(self_table.keys()):
        assertion = "sequence length (" + str(len(sequence)) + ") does not match the provided energy table (" + str(len(self_table.keys())) + ")"
        raise AssertionError(assertion)

    energy_score = 0.0

    #sum the self energies
    for index in self_table.keys():
        try:
            energy_score += self_table[index][sequence[index]]
        except:
            continue

    #sum the pair energies
    for pair_index in pair_table.keys():
        try:
            energy_score += pair_table[pair_index][(sequence[pair_index[0]],sequence[pair_index[1]])]
        except:
            continue

    return(energy_score)

def readSequences(sequence_path):
    sequences = list()
    with open(sequence_path,'r') as file:
        sequences = [line.rstrip('\n') for line in file] #reads in all of the lines from the etab
    return sequences

def writeScores(output_path,sequences,self_table,pair_table):
    with open(output_path,"w") as outputFile:
        for sequence in sequences:
            score = scoreSequence(sequence,self_table,pair_table)
            outputFile.write(sequence+'\t'+str(score)+'\n')
    return


def main():
    #parse the arguments
    parser = argparse.ArgumentParser(description='Score sequences using an dTERMen energy table')
    parser.add_argument("--etab",help="the path specifying the energy table that will be used for scoring", required=True)
    parser.add_argument("--sequences",help="the path specifying the list of sequences to be scored. This assumes each line is a new sequence.", required=True)
    parser.add_argument("--output",help="the path specifying where the scores should be output to.", required=True)
    args = parser.parse_args()
    args_dict = vars(args)

    #generate a list of all of the sequences in the provided file
    sequences = readSequences(args_dict["sequences"])

    #generate dictionaries storing self and pair energies for fast access
    (self_table,pair_table) = loadDTERMenEnergyTable(args_dict["etab"])

    #score the provided sequences and generate an output file
    writeScores(args_dict["output"],sequences,self_table,pair_table)
    return

# finally, the main function is called if this is run in terminal
if __name__ == "__main__":
    main()
