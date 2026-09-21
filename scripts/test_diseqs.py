from diseqs import main, open_file, filt, convert, Genome
import pytest
import pandas


def test_main():
    # Correct
    assert main("test_files/sequence.fasta") == {"GS874858.1": 37.537}

    # Wrong file name
    with pytest.raises(FileNotFoundError):
        main("test_files/sequencesecuence.fasta")

    # Correct class name
    assert type(main("test_files/sequencesequence.fasta")) == dict


def test_open_file():
    # Correct fasta open
    assert open_file("test_files/sequencesequence.fasta") == {
        "GS874858.1": "TAGTGTAACTGGGTTGACGTTCCATGTAGCAAATACGTCTCTAGCTTTAATTACCTTATTGTAATCATTGACAGTTCCTTTTGGAAGATTTATAGTTACTCTTCCAGAAGAGGTATTAATAGCGTATGATTTTCCCCATTCGGCTTTTAGTGTTTGACCAGATGAAGCATCATAAGTTTTCCAGGCACCGGCTGAATATGGAACATCTCCATCACCAAGCTCGTAATAAAGCTCATCAAAGTTTTCATTTATTTTTATACCACCTTTACGCAGGTAGTCACCGGTACCATCATCAACAACATTACCGATATTAATATTTTGTTTCATTATTGAGCCACCCC"
    }

    # Wrong fasta file written
    with pytest.raises(SystemExit):
        open_file("test_files/sequencewrong_written.fasta")

    # Correct csv open
    assert type(open_file("test_files/sequencesequence.csv")) == pandas.DataFrame

    # Correct tsv open
    assert type(open_file("test_files/sequencegenome.tsv")) == pandas.DataFrame


def test_filt():
    # Correct leng filter
    seq = open_file("test_files/sequencesequence.csv")
    assert filt(data=seq, leng=3, id="qacc", column="sacc") == {
        "Query_5958997": "PP987310.1"
    }

    # Correct quality filter
    assert type(filt(data=seq, quality=True)[1]) == pandas.DataFrame

    # Correct quality filter
    assert type(filt(data=seq, quality=True)[0]) == pandas.DataFrame

    # Correct filter working
    with pytest.raises(SystemExit):
        filt(data=seq, leng="hola")

    with pytest.raises(SystemExit):
        filt(data=seq, quality="hola")

    with pytest.raises(SystemExit):
        filt(data=seq, leng=3, id=3, column="sacc")


def test_convert():
    seq = open_file("test_files/sequencesequence.fasta")
    # converts makes genome class
    assert type(convert(seq)) == Genome

    # Count nb
    assert convert(seq).nb_count == {
        "GS874858.1": [{"A": 95}, {"T": 118}, {"C": 73}, {"G": 55}]
    }

    # Length
    assert convert(seq).length == {"GS874858.1": 341}

    # Correct GC count with 3 decimals
    assert convert(seq).GC_content(an="%", rd=3) == {"GS874858.1": 37.537}

    # Correct RNA reverse convert
    assert convert(seq).RNA(reverse=True) == {
        "GS874858.1": "UAGUGUAACUGGGUUGACGUUCCAUGUAGCAAAUACGUCUCUAGCUUUAAUUACCUUAUUGUAAUCAUUGACAGUUCCUUUUGGAAGAUUUAUAGUUACUCUUCCAGAAGAGGUAUUAAUAGCGUAUGAUUUUCCCCAUUCGGCUUUUAGUGUUUGACCAGAUGAAGCAUCAUAAGUUUUCCAGGCACCGGCUGAAUAUGGAACAUCUCCAUCACCAAGCUCGUAAUAAAGCUCAUCAAAGUUUUCAUUUAUUUUUAUACCACCUUUACGCAGGUAGUCACCGGUACCAUCAUCAACAACAUUACCGAUAUUAAUAUUUUGUUUCAUUAUUGAGCCACCCC"
    }
