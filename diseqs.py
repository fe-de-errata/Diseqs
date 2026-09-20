# Libreries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Global variables
bases = {"A": "T", "T": "A", "C": "G", "G": "C"}


# Funciones


# Main function, here you choose the analysis
def main(file):
    """Type the name of your data file"""
    data = open_file(file)
    return convert(data)


# Open file function (fasta)
def open_file(file=""):
    # Open fasta files
    if file.endswith((".fasta", ".fna")):
        data = {}
        with open(file) as f:
            for i in f:
                i = i.rstrip()
                try:
                    if ">" in i:
                        header = i.replace(">", "").split(" ")[0]
                        s = ""
                    else:
                        s += i
                        data[header] = s
                except UnboundLocalError:
                    print("Wrong file syntaxis")
                    sys.exit(1)
    # open csv or tsv files
    elif file.endswith(".csv"):
        data = pd.read_csv(file, engine="python")
    elif file.endswith(".tsv"):
        data = pd.read_csv(file, sep=r"\t", engine="python")
    return data


# Writes the output of your analysis
def make_graph(skews):
    fig = plt.figure(figsize=(13,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axhline(0, color='black', linewidth=1.2, linestyle='-')
    ax.axvline(0, color='black', linewidth=1.2, linestyle='-')
    ax.plot(range(len(skews)), skews, color = 'olive', lw = 1)
    return plt.savefig("diseqs_image.png")


# This returns a dict filtered by length, a dataframe filtered and sorted by the arguments given
def filt(data="", leng="", dataframe="", dropdup="", quality=""):
    try:
        if not data:
            raise ValueError("Not data given")
        if dataframe != "" or not isinstance(dataframe, str):
            raise ValueError("Invalid dataframe")
        if leng != "" or not isinstance(leng, int):
            ValueError("Leng must be int only")
        for i in [dropdup, quality]:
            if i != "" or not isinstance(i, bool):
                raise ValueError(f"{i} must be bool")
    except ValueError:
        sys.exit(1)
    # Returns the dict of the data given filtered by the length given
    if leng:
        LENG = {}
        for i in data.keys():
            if len(data[i]) >= leng:
                LENG.updat({i: data[i]})
    if dataframe:
        if quality:
            filtered = []
            not_filtered = []
            for qacc, identity, coberture, evalue, bitscore, qlen in dataframe[
                ["qacc", "pident", "length", "evalue", "bitscore", "qlen"]
            ].values:
                cobertura = (coberture / qlen) * 100
                res = ((identity >= 90) or (cobertura >= 75)) and (
                    (evalue <= 1e-5) or (bitscore >= 50)
                )
                if res == True:
                    filtered.append(qacc)
                else:
                    not_filtered.append(qacc)
            DF_FALSE = dataframe[dataframe["qacc"].isin(not_filtered)]
            DF_TRUE = dataframe[dataframe["qacc"].isin(filtered)]
    return [LENG, DF_FALSE, DF_TRUE]


# Makes a orginal dataframe with the data given, can merge only 2 dataframes
def df_maker(data="", dataframes="", merge="", column="", how=""):
    try:
        if not data:
            raise ValueError("Data not given")
        for i in [column, how]:
            if i != "" or not isinstance(i, str):
                raise ValueError("Invalid columnns to merge")
        if dataframes != "" or not isinstance(dataframes, list):
            raise ValueError("Need to give dataframes to work on")
        if merge != "" or not isinstance(merge, bool):
            raise ValueError("Merge must be True or False")
    except ValueError:
        sys.exit(1)
    # Make a dataframe of the data given
    DF_OG = pd.DataFrame(data)
    if merge:
        MERGE = dataframes[0].merge(dataframes[1], on=column, how=how)
    return [DF_OG, MERGE]


# Turns de data into a Genome class
def convert(data):
    conv = Genome(data)
    return conv


# Classes
# Genome class only admits dict, only gene_search admits list
class Genome:
    def __init__(self, data):
        self.data = data

    @property
    # Returns de length of each value in a dict
    def length(self):
        lengths = {}
        for i in self.data.keys():
            lengths.update({i: len(self.data[i])})
        return lengths

    @property
    # Counts the nucleatids of each value in a dict
    def nb_count(self):
        counts = {}
        for i in self.data.keys():
            nucl = []
            if "U" in self.data[i]:
                bases.update({"U": "T"})
                del bases["T"]
            for nb in bases.keys():
                nucl.append({nb: self.data[i].count(nb)})
            counts.update({i: nucl})
            if "T" not in bases.keys():
                bases.update({"T": "A"})
                del bases["U"]
        return counts

    @property
    # Counts the N of each value in a dict
    def N(self):
        N_count = []
        for i in self.data.keys():
            N_count.append((i, self.data[i].count("N")))
        return N_count

    # Convcerts each value of a dict into RNA, can convert to the complement and reverse secuence
    def RNA(self, complement="", reverse=""):
        try:
            if (complement != "" and not isinstance(complement, bool)) or (
                reverse != "" and not isinstance(reverse, bool)
            ):
                raise ValueError("Complement and reverse must be '', True or False")
        except ValueError:
            sys.exit(1)
        rna = {}
        for i in self.data.keys():
            try:
                if "U" in self.data[i]:
                    raise ValueError(f"The secuence {i} is already RNA")
            except ValueError:
                sys.exit(1)
            seqc = self.data[i].replace("T", "U")
            rna.update({i: seqc})
        if complement == True:
            temp_bases = dict(bases)
            temp_bases.update({"U": "A"})
            for i in rna.keys():
                seqc = ""
                for j in rna[i]:
                    if j != "A":
                        seqc = seqc + temp_bases[j]
                    else:
                        seqc = seqc + "U"
                if reverse == True:
                    rna.update({i: seqc[::-1]})
                else:
                    rna.update({i: seqc})
            result = rna
        else:
            result = rna
        return result

    # Counts the GC content of each value of a dict, can return the total of GC or de percentaje of GC, can give a numer to round
    def GC_content(self, an="", rd=""):
        try:
            if an not in ["", "%"]:
                raise ValueError("Invalid argument", an, ",only '' or '%'")
            if rd != "" and not isinstance(rd, int):
                raise ValueError("Invalid argument", rd, ",'rd' must be int")
        except ValueError:
            sys.exit(1)
        gc = {}  # "ID" : ["G", "C"]
        for i in self.data.keys():
            gc.update({i: [self.data[i].count("G"), self.data[i].count("C")]})
        if an == "%":
            for i in gc.keys():
                if rd == "":
                    final = gc.update(
                        {i: round(((gc[i][0] + gc[i][1]) / self.length[i]) * 100, 2)}
                    )
                else:
                    final = gc.update(
                        {i: round(((gc[i][0] + gc[i][1]) / self.length[i]) * 100, rd)}
                    )
            final = gc
        else:
            final = gc
        return final

    # Calculates GC_skew in order to create a graph
    def GC_skew(self, nwin=100):
        try:
            if len(self.data) != 1:
                raise ValueError("GC_skew only admits one secuence")
        except ValueError:
            sys.exit(1)
        GC_skews = {}
        for i in self.data.keys():
            L = len(self.data[i])
            v = max(L // nwin, 1)
            d = v
            N = int((L - v) / d) + 1
            D = 0
            positions, GCSK, GCc = [], [], []
            G, C = self.data[i].count("G"), self.data[i].count("C")
            GC_glob = ((G + C) / L) * 100
            for j in range(N):
                window = self.data[i][D : D + v]
                G, C = window.count("G"), window.count("C")
                if (G + C) > 0:
                    GC_skew = (G - C) / (G + C)
                else:
                    GC_skew = 0
                if len(window) > 0:
                    GC_content = ((G + C) / len(window)) * 100
                else:
                    GC_content = 0
                GC_content_norm = GC_content - GC_glob
                positions.append(D)
                GCSK.append(round(GC_skew, 3))
                GCc.append(round(GC_content_norm, 3))
                D += d
            GC_skews.update({i: [positions, GCSK, GCc]})
        return GC_skews

    # The only method that admits list, it searchs the list o dict of secuences in the data values
    def gene_search(self, items):
        try:
            if not items or not isinstance(items, (list, dict)):
                raise ValueError(
                    "Invalid arguments", items, "itemes must be list or dict"
                )
        except ValueError:
            sys.exit(1)
        mtchs = {}
        if type(items) == dict:
            for i in self.data.keys():
                for j in items.keys():
                    if items[j] in self.data[i]:
                        mtchs.update({i: (j, items[j])})
        elif type(items) == list:
            for i in self.data.keys():
                for j in items:
                    if j in self.data[i]:
                        mtchs.update({i: j})
        return mtchs


if __name__ == "__main__":
    main()
