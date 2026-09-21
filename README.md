# Diseqs

Diseqs is a basic library for bioinformatics python scripts, includes the estimate of GC Skew, GC content, dataframes maker, etc.

## How to install

```bash
git clone https://github.com/fe-de-errata/Diseqs.git
cd Diseqs
pip install -r requirements.txt
```

## How to use

To import all

```python
import diseqs as ds
```

or

To import a special function

```python
from diseqs import open_file
```

## Functions

### Summary

- main()
- open_file()
- make_graph()
- filt()
- df_maker()
- convert()

#### main(file)

The main function opens the file and converts into the Genome class and returns the GC content percentage with 3 decimals.

Example:

```python
from diseqs import main

print(main(file="human.fasta"))

>>> e.g. {"ID": 23.333}
```

Requirements:

- 1 argument: filename
- main functions only works with fasta files, look at classes.Genome to read more.

#### open_file(file='')

The open_file funcion opens fasta, csv and tsv files. Use the python engine to open the file and reads the end of the filename to open it.
This function returns a dict is the file given is fasta. If the file is a csv or tsv it returns a dataframe.

Example:

```python
from diseqs import open_file

genome = open_file(file="human_genome.fasta")

print(genome)

>>> {"ID": "ATCGATGGTCAGT"}
```

Requirements:

- 1 argument: filename
- Only admits fasta, csv or tsv files

#### make_graph(skew, image = 'new_image.png')

The make_graph function makes a 2 dimension linear, graph of the GC Skew of the GC content given as the argument "skew".
This functions creates a file, by defualt called 'new_image.png'.

Look that Genome.GC_skew returns a dict with:
{"ID": ["positions", "GC skews", "GC content"]}

Example:

```python
import diseqs as ds

#This example contemplates only 1 sequence
gc_skews = (ds.convert(ds.open_file("human.fasta"))).GC_skew 
for i in gc_skews.keys():
    id = i
make_graph(gc_skews[id][2])

>>> open image.png
```

Requirements:

- 2 argument required: the list of GC skew, and the name of the image ('new_image.png' by default)
- This functios uses the Genome class to get the GC skew. To read more look at classes.Genome.

#### filt(data="", leng="", quality="", id= "", column="")

This function is a customizable filter. Can give a length to filter with or/and arguments to make dataframes filtered.
Returns LENG if only leng and column arguments given. Else if quality given it returns a list: [DF_FALSE, DF_TRUE].

Example:

```python
import diseqs as ds

df = ds.open_file("genome.csv")
leng_sorted = ds.filt(data=df, leng=10, id="gene", column="sequence")

print (leng_sorted)

>>> gene : sequence
    "gagf" : "ATCGATGCAT"
```

Requierements:

- 5 arguments:
  - data='' (neccesary).
  - leng='' (optional): can sort a dict by a length given
  - quality='' (optional, True or False): activates a quality filter.
  - id='' (optional, str): give the column ID of the dataframe.
  - column=''(optional, str): give a column to sort by length.
- Only admits dataframes from open_file()

#### df_maker(data="", merge="", column="", how="")

This functions makes dataframes and merge them.
Returns a single dataframe if only data given, else if merge given it returns a merge dataframe.

```python
import diseqs as ds

genome = ds.open_file("genome.tsv") #with ID : sequences
comp_genome = ds.open_file("complement.tsv") #with ID : Organism
merge = ds.df_maker(data=[genome, comp_genome], merge=True, column="ID", how="left")

print(merge)

>>> ID  Sequence    Organism
    "AHE" "ATCGGCGATG   Homo sapiens sapiens
```

Requieremnts:

- 4 arguments:
  - data=''(neccesary): can be dict, the output of open_file() or a list with 2 dataframes.
  - merge=''(optional, True or False): indicates to merge dataframes given.
  - column=''(optional, str): gives a columns to merge with.
  - how=''(optinal, 'left' or 'right'): the side to merge the dataframes.
- If not merge needed, it can make dataframes by a dict given in data.
- Can only merge 2 dataframes.

#### convert(data)

This function converts a data given into Genome class.

Example:

```python
import diseqs as ds

data = ds.open_file("genome.fasta")
gen = ds.convert(data)

print(type(gen))

>>> type Genome
```

Requirements:

- Only admits fasta files.

## Classes

### Summary

- Genome

#### Genome

This class has methos and funcionts and only admits fasta files opened as dicts.

##### Methods

- Genome.length
    Returns the length of each value of a dict:
        {"ID": length}

- Genome.nb_count
    Returns the count of each nucleotid in each value of a dict:
        {"ID": [{"A": #}, {"T": #}, {"G": 3}, {"C": #}]}

- Genome.N
    Returns the count of N in each value of a dict:
        {"ID": #N}

##### Functions

- Genome.RNA(complement="", reverse="")
    Converts DNA sequence to RNA if not argument given.
    If complement=True: returns the RNA complementary.
    If reverse=True: returns the sequence reversed.

- Genome.GC_content(an="", rd="")
    Returns the GC content as list if no arguments given:
        {"ID": [#G, #C]}
    If an="%": returns the percentage of GC in sequence.
    If rd=#: rounds the percentege decimals the number given.

- Genome.GC_skew(nwin=100)
    Returns de GC skew and GC content as a dict:
        {"ID": [positions, GC skew, GC content]}
    nwin= modifys the window for GC skew count.

- Genome.gene_search(items)
    Returns a dict with matches:
        {"ID": [match1, match2]}
    Requieres a item to search in data.
    Admits list or dict for the itmes to search.

## System Requierements

- Python 3.13.5 (and above)
- requirements.txt libreries installed
