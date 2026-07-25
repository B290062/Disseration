import pytest
import subprocess
import sys

#the path of the python pipeline that is being tested, this will need to be changed by the user.
code = "/home/s2089123/Dissertation/Pipeline.py"
#code was loosely adapted from this source https://pythontest.com/testing-argparse-apps/ it displays the use of subprocess and pytesting in action
#which was helpful in designing tests. 
#in this example, they imported main, however as my command line interface isn't in main its easier to use the file path for testing.

#the following test specifies a mode that doesn't exist in the command line and makes sure the program fails successfully
def test_wrong_mode():
    command = subprocess.run(["python3", code, "--sra", "SRP091444", "--mode", "somerandomwords", "--mask", "mm39_refseq.bed"],capture_output=True, text=True)
    assert command.returncode != 0

#example adapter sequences were sourced from illuminas website https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
#this test makes sure the program exists when the mode is rnaseq and only one adapter is provided.
def test_one_adapter_rnaseq():
    command = subprocess.run(["python3", code, "--sra", "SRP091444", "--mode", "rnaseq", "--mask", "mm39_refseq.bed", "--trim", "--adapter1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],capture_output=True, text=True)
    assert command.returncode != 0 

#tests that when trim is used, and no adapters are entered that the code exits
def test_no_adapters_groseq():
    command = subprocess.run(["python3", code, "--sra", "SRP064201", "--mode", "groseq", "--mask", "mm39_refseq.bed", "--trim"],capture_output=True, text=True)
    assert command.returncode != 0

#this test checks whether the code exists correctly when no SRA study number is provided in the command line
def test_no_SRA():
    command = subprocess.run(["python3", code, "--sra", "--mode", "rnaseq", "--mask", "mm39_refseq.bed"],capture_output=True, text=True)
    assert command.returncode != 0

def test_negative_window():
    command = subprocess.run(["python3", code, "--sra", "SRP091444", "--mode", "groseq", "--mask", "mm39_refseq.bed", "--window", "-9"],capture_output=True, text=True)
    assert command.returncode != 0

#this test checks if the code exits correctly when no bed file coordinates are provided in the command line"
def test_no_mask():
    command = subprocess.run(["python3", code, "--sra", "SRP091444", "--mode", "groseq", "--mask"],capture_output=True, text=True)
    assert command.returncode != 0

#this test makes sure that the help function is working in the command line to make sure the user is able to see how the pipeline functions before running it
def test_help_function():
    command = subprocess.run(["python3", code, "--help"],capture_output=True, text=True)
    assert command.returncode == 0



#this test also tests the full pipeline but for gro-seq it can be uncommented out when necessary.    
def test_full_pipeline_gro():
    command = subprocess.run(["python3", code, "--sra", "SRP064201", "--mode", "groseq", "--mask", "mm39_refseq.bed"])
    assert command.returncode == 0