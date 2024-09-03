# Scientific programming                     Prof. Piro Rosario Michael
# Year:                                      2023/2024
# Python project                             N.6 FM Index
# Python version                             3.11.9
# Author                                     Nossa Leonardo
# GitHub URL                                 https://github.com/LeonardoNossa/FM-index-prj

#!/usr/bin/env python3 

import requests 
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import time
import progressbar
import textwrap


class MissedOperation(Exception):
    pass

class ReadError(Exception):
    pass

class IntervalError(Exception):
    pass

class ElementNotFound(Exception):
    pass

class FM_index():

    def __init__(self):
        """
        Create the construct of the class.
        Here, there are all the call of the function inside the class
        """

        self.bar = progressbar.ProgressBar(max_value=9)
        self.bar.start()

        self.directory = self.dir_creation()
        self.bar.update(1)

        self.fasta, self.header = self.readFile()
        self.bar.update(2)

        self.fasta_sd = self.single_file_directory()
        self.sa = self.SuffixArray()
        self.bar.update(3)

        self.bwt = self.BWTransform()
        self.bar.update(4)

        self.LF_array = self.C_c()
        self.bar.update(5)

        self.occ = self.Occ_matrix()
        self.bar.update(6)

        self.pt_list = self.pattern_input()

        if self.pt_list:  
            self.search_result = self.final_result()
            self.bar.update(7)

            self.pattern_plot = self.DNA_plot()
            self.bar.update(8)

        self.save_all = self.save()
        self.bar.update(9)

    def dir_creation(self):
        """
        Create the source directory for all the results/structures

        Results:
            creation of the directory where put all the results obtained by the other functions
        """
        
        desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
        dir = os.path.join(desktop_path, "FM_index_cache")
        os.makedirs(dir, exist_ok=True)

        time.sleep(1)
        if os.path.isdir(dir):
            print(" Opening the directory, FM_index_cache")
        else:
            print(" Creation of the directory, FM_index_cache")

        return dir

    def readFile(self):
        """
        Reads and processes a FASTA sequence from either a local file or an online source.

        Returns:
            tuple: The sequence and header information.

        Raises:
            ValueError: If the interval values are not integers.
            IntervalError: If the specified interval is out of the sequence range.
            MissedOperation: If an invalid operation is provided.
            FileNotFoundError: If the online sequence retrieval fails.
            ReadError: If an invalid option for reading is given.
        """

        if len(sys.argv) < 4:
            raise ValueError("Insufficient arguments, must be at least 4.")

        req = sys.argv[1]
        interval = sys.argv[3]
        seq = ""
        header = None

        if req == "1":
            fasta = sys.argv[2]
            if not os.path.isfile(fasta):
                raise FileNotFoundError(f"The loaded FASTA file doesn't exist: path '{fasta}'")
            if interval == "full":
                with open(fasta, "r") as f:
                    for line in f:
                        if line.startswith(">"):
                            header = line.strip()
                        else:
                            seq += line.strip()
                seq += "$"
            elif interval == "personalized":
                try:
                    start = int(sys.argv[4])
                    stop = int(sys.argv[5])
                except ValueError:
                    raise ValueError("Bad interval value. Must be an integer number.")

                with open(fasta, "r") as f:
                    for line in f:
                        if line.startswith(">"):
                            header = line.strip()
                        else:
                            seq += line.strip()
                    if start > len(seq) or stop > len(seq):
                        raise IntervalError("Out of range. Check the start and/or end points.")
                seq = seq[start:stop] + "$"
            else:
                raise MissedOperation("Must be one of these two values: 'full' or 'personalized'.")

        elif req == "2":
            sequence_id = sys.argv[2]
            if interval == "full":
                start = None
                stop = None
            elif interval == "personalized":
                try:
                    start = int(sys.argv[4])
                    stop = int(sys.argv[5])
                except ValueError:
                    raise ValueError("Bad interval value. Must be an integer number.")
            else:
                raise MissedOperation("Must be one of these two values: 'full' or 'personalized'.")

            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                "db": "nuccore",
                "id": sequence_id,
                "seq_start": start,
                "seq_stop": stop,
                "rettype": "fasta",
                "retmode": "text"
            }

            response = requests.get(url, params=params)

            if response.status_code == 200:
                lines = response.text.strip().split("\n")
                header = lines[0]
                seq = "".join(lines[1:])
                seq += "$"
            else:
                raise FileNotFoundError(f"Error: {response.status_code}")
        else:
            raise ReadError("Invalid option. This argument must be '1' or '2' depending on the data source (local file or online).")

        time.sleep(2)
        print("Loading the FASTA file")

        return seq, header

    
    def single_file_directory(self):
        """
        Create the directory for the loaded sequence inside the source directory

        Results:
            the directory created inside the main directory with the correct name(if it doesn't exist)
            - all the structures/inìformation arew stored in this directory -
        """

        if self.header is not None:
            fasta_dir = os.path.join(self.directory, self.header[1:].replace(":", "_"))
        else:
            fasta_dir = os.path.join(self.directory, sys.argv[2])

        if sys.argv[3] == "personalized":
            start = sys.argv[4]
            stop = sys.argv[5]
            fasta_dir = os.path.join(self.directory, f"{self.header[1:].split(':')[0]}_{start}_{stop}")
    
        if not os.path.exists(fasta_dir):
            os.mkdir(fasta_dir)
        
        return fasta_dir
    
    def SuffixArray(self):
        """
        Create and return the Suffix Array (SA)

        Returns:
            the single tuple with the data but the final implementation of the structure is in the save function and it is a dataframe where the first column is the indices of FASTA and then second
            one is the array sequence (sub-sequence of the FASTA)
        """

        idx_SA = sorted(range(len(self.fasta)), key=lambda i: self.fasta[i:])
        sa = [(i, self.fasta[i:]) for i in idx_SA]

        time.sleep(3)
        print("Suffix array creation")
        
        return sa
    
    def BWTransform(self):
        """
        Create and return the Burrow Wheeler Transform (BWT), also called block-sorting compression

        Returns:
            the FASTA file is processed with an algorithm that restructures data in such a way that the transformed message is more compressible.
            Technically, it is a lexicographical reversible permutation of the characters of as string 
        """

        n = len(self.fasta)
        rotations = [self.fasta[i:] + self.fasta[:i] for i in range(n)]
        rotations_sorted = sorted(rotations)

        bwt = "".join(rotation[-1] for rotation in rotations_sorted)
        
        time.sleep(4)
        print("Doing the Burrow Wheeler Transform")

        return bwt
    
    def C_c(self):
        """
        Create and return the C[c] table

        Returns:
            the C[c], that is a table that, for each character c in the alphabet, contains the number of occurrences of lexically smaller characters in the text
        """
        
        bwt_sorted = sorted(self.bwt)
        C_data = {ch: bwt_sorted.index(ch) for ch in set(bwt_sorted)}

        data = {" ": list(C_data.values())}
        df_C = pd.DataFrame(data, index=[k for k in C_data.keys()])

        time.sleep(5)
        print("Creation of C[c]")

        return df_C

    def Occ_matrix(self):
        """
        Create and return the Occ(c,k) table

        Returns:
            the Occ(c,k) is a matrix, is the number of occurrences of character c in the prefix L[1..k]
        """

        C_string = sorted(set(self.bwt))
        Occ_matrix = pd.DataFrame(0, index=C_string, columns=range(len(self.bwt) + 1))

        for i in range(1, len(self.bwt) + 1):
            for l in C_string:
                if self.bwt[i - 1] == l:
                    Occ_matrix.loc[l, i] = Occ_matrix.loc[l, i - 1] + 1
                else:
                    Occ_matrix.loc[l, i] = Occ_matrix.loc[l, i - 1]
        Occ_matrix = Occ_matrix.drop(0, axis=1)

        time.sleep(6)
        print("Creation of Occ(c, k)")

        return Occ_matrix

    def pattern_input(self):
        """
        This feature captures and verifies pattern input based on the reading mode you choose.
        """

        if sys.argv[3] == "full":
            pattern = sys.argv[4:]  
        elif sys.argv[3] == "personalized":
            pattern = sys.argv[6:]    
        else:
            raise ValueError("The third argument must be full or personalized")

        for pt in pattern:
            if not isinstance(pt, str):
                raise ValueError("Pattern must be a string")
            if not pt.isalpha():
                raise ValueError("Patterns cannot contain numbers and symbols ")

        return pattern

    def backward_search_pattern(self):
        """
        Performs pattern search using backward search based on the FM-index.

        Returns:
            A list with:
            1) The index of the first occurrence of the pattern in the sequence.
            2) The list of positions of the pattern in the sequence.
            3) The total number of occurrences of the pattern in the sequence.
        """

        pattern_list = [pattern.upper() for pattern in self.pt_list]

        find1 = []
        find2 = []
        find3 = []

        for pattern in pattern_list:
            start = 0
            end = len(self.bwt) - 1

            for i in range(len(pattern) - 1, -1, -1):
                char = pattern[i]
                
                if char in self.LF_array.index:
                    occ_start = self.occ.loc[char, start] if start > 0 else 0
                    occ_end = self.occ.loc[char, end + 1] if end + 1 <= len(self.bwt) else self.occ.loc[char, len(self.bwt)]

                    start = self.LF_array.loc[char, " "] + occ_start
                    end = self.LF_array.loc[char, " "] + occ_end - 1

                    if start > end:
                        break
                else:
                    start = end + 1
                    break

            if start <= end:
                positions = [self.sa[i][0] for i in range(start, end + 1)]  
                find1.append({"First occurrence index": min(positions)})
                find2.append({"Pattern position": positions})
                find3.append(len(positions)) 
            else:
                find1.append({"First occurrence index": "Not found"})
                find2.append({"Pattern position": "Not found"})
                find3.append(0)

        return find1, find2, find3

    def final_result(self):
        """
        Create a method that stores the results of the backward search.

        Returns:
            A dataframe containing the results.
        """

        first_idx, all_idx, occurrences = self.backward_search_pattern()

        final_results = []
        for i in range(len(self.pt_list)):
            info = {
                "First occurrence index": first_idx[i].get("First occurrence index", "Not found"),
                "Pattern position": all_idx[i].get("Pattern position", "Not found"),
                "Total occurrences": occurrences[i]
            }
            final_results.append(info)
            
        data_results = pd.DataFrame(final_results)
        data_results.index = [pattern.upper() for pattern in self.pt_list]

        time.sleep(7)
        print("Searching your patterns")

        return data_results

    
    def DNA_plot(self):
        """
        Create a way for trace the pattern on the sequence 

        Returns:
            printing a plot that spot the pattern belong the sequence
        """

        pattern_dir = os.path.join(self.fasta_sd, "Pattern founded")
        os.makedirs(pattern_dir, exist_ok=True)
        
        pattern_list = [pattern.upper() for pattern in self.pt_list]
        seq = self.fasta[:-1]
        sequence_lines = textwrap.wrap(seq, 124)

        for pt in pattern_list:
            plt.figure(figsize=(12.85, 5), constrained_layout=True)
            if pt in self.search_result.index:
                positions = self.search_result.loc[pt, "Pattern position"]
                if positions != "Not found":
                    pattern_legend = []
                    for pos in positions:
                        if pt not in pattern_legend:
                            plt.axvline(x=pos, color="red", linestyle="dotted", linewidth=1, label="pattern")
                            pattern_legend.append(pt)
                        else:
                            plt.axvline(x=pos, color="red", linestyle="dotted", linewidth=1)
        
            base_color = {"A":"#d9d9d9", "T":"#bdbdbd", "C":"#969696", "G":"#737373"}
            base_legend = []
            for i, base in enumerate(seq):
                if base not in base_legend:
                    plt.bar(i, 1, color=base_color.get(base), label=base)
                    base_legend.append(base)
                else:
                    plt.bar(i, 1, color=base_color.get(base))
            
            displayed_lines = "\n".join(sequence_lines[:2] + ["...\n"] + sequence_lines[-2:]) if len(sequence_lines) > 3 else "\n".join(sequence_lines)
            if self.header:
                plt.title(self.header + "\n" + f"Searching the pattern: {pt}")
                if self.header[1:].startswith("NG"):
                    plt.ylabel("DNA sequence")
                elif self.header[1:].startswith("NM"):
                    plt.ylabel("mRNA sequence")
                elif self.header[1:].startswith("NR"):
                    plt.ylabel("non-coding RNA sequence")
            else:
                plt.title(sys.argv[2] + "\n" + f"Searching the pattern: {pt}")
                plt.ylabel("Sequence")

            plt.text(0.045, -0.02, displayed_lines, ha="left", va="top", family="monospace", transform=plt.gca().transAxes)
            plt.tight_layout()
            plt.xticks([])
            plt.yticks([])
            plt.ylim(0, 1.5)
            plt.legend()
            plt.savefig(os.path.join(pattern_dir, f"pattern {pt}.png"), dpi=300)
            plt.close()

        time.sleep(8)
        print("Representing the patterns")

    def save(self):
        """
        Saving all the results/infdormation in single files inside the source directory 
        """

        file_name = os.path.join(self.fasta_sd, "FASTA.txt")
        with open(file_name, "w") as f:
            if self.header:
                f.write(f"{self.header}\n{self.fasta[:-1]}")
            else:
                f.write(f"{sys.argv}\n{self.fasta[:-1]}")
        
        data = [{" ": tupla[1]} for tupla in self.sa]
        df = pd.DataFrame(data, index=[tupla[0] for tupla in self.sa])
        path_sa = os.path.join(self.fasta_sd, "SuffixArray.csv")
        df.to_csv(path_sa, index=True)

        bwt_file = os.path.join(self.fasta_sd, "BWT.txt")
        with open(bwt_file, "w") as f2:
            f2.write(f"{self.bwt}")
        
        path_LF_array = os.path.join(self.fasta_sd, "C[c].csv")
        self.LF_array.to_csv(path_LF_array, index=True)

        path_Occ = os.path.join(self.fasta_sd, "Occ(c,k).csv")
        self.occ.to_csv(path_Occ, index=True)

        if hasattr(self, 'search_result'):
            path_results = os.path.join(self.fasta_sd, "TableOfPattern.csv")
            if os.path.exists(path_results):
                self.search_result.to_csv(path_results, mode='a', header=False, index=True)
            else:
                self.search_result.to_csv(path_results, mode='w', header=True, index=True)
        
        time.sleep(9)
        print("Saving all the information")

fm_index = FM_index()
