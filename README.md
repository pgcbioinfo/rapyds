
## Restriction Site Associated DNA Python-Digested Simulation  (RApyDS)


**RApyDS** is a python script that performs in silico digestion as an aid for choosing a restriction enzyme for RADseq experiments. 

**RApyDS** can perform in silico digestion of a given genome/sequence or gc content. It can simulate single-enzyme or double-enzyme digestion on a single file or multiple FASTA files. After digestion, the program can match fragments to a given an annotation file and run alignment analysis to check for unique and repeat regions.

**RApyDS** provides a detailed report and visualization of the simulated digestion. 


### Documentation
- For the [full RApyDS program technical manual](docs/rapyds_manual.pdf)

### Requirements

- Python 3 or greater
- pip3 installed (run: `pip3 install -r requirements.txt`)
- BWA 0.7.12 (http://bio-bwa.sourceforge.net/)
- Firefox (for viewing the html files)
- Linux OS

### Usage

#### Arguments
-  **-h, --help**  show this help message and exit
-  **-gc [GC]**    input gc frequency. Value must be between 0 and 1
-  **-dna [DNA]**  input dna estimated length
-  **-i [I]**      directory containing the input files
-  **-pre [PRE]**  prefix of the input files (must match the file name of the sequence, annotation, and/or index files)
-  **-at [AT]**    what to look for in gene annotation file (ex. gene region, exon, intron, etc) (default: gene)
-  **-db [DB]**    resteriction enzyme dabatase file. Format per line: SbfI,CCTGCA|GG (default: database/re_db.txt)
-  **-re [RE]**    file of list of restriction enzyme to be tested
-  **-min [MIN]**  minimum fragment size (default: 200)
-  **-max [MAX]**  maximum fragment size (default: 300)
-  **-bp [BP]**    base pair read length for FASTQ generation (default: 100)
-  **-p [P]**      radseq protocol: use ddrad for double digestion (default: orig)
-  **-o [O]**      output file name (default: report)
-  **-t [T]**      number of processes (default 16)

Optional Flags:
-  **--bwaskip**   skip BWA indexing and alignment
-  **--clean**     clean files after running



#### Common Usage
- To use RApyDS with a given genome file (and/or corresponding annotation file)

``python rapyds.py -i <input_directory> -pre <input_prefix> [other arguments/flags]``

- To use RApyDS by generating a genome file given GC content/frequency

``python rapyds.py -gc <gc frequency from 0 to 1> -dna <length> [other args]``


#### Sample Run
Given an E.Coli FASTA file ``ecoli_seq.fasta`` with annotation file ``ecoli_seq.gff`` both located inside the directory ``ecoli``, the RADSeq protocol is DDRad

``python rapyds.py -i ecoli -pre ecoli_seq -p ddrad -re enzyme/dd_re.txt``


#### Output

The output files of the program is a zip file named after the ``-o`` argument or by default, the ``report.zip``.
Inside the archive are 3 html files containing:
- Overview (``index.html``) - summary of the results
- Electrophoresis (``gel.html``) - electrophoresis simulation of up to 5 restriction enzymes of choice
- Cut Site Distribution (``cutsite.html``) - zoomable graphic vector of cut site locations in the sequence

> It is advisable open the html files using Firefox. There is an issue with Google Chrome when opening local files. Only the overview file will work fine for any browser.

> Slowdown may be experienced when loading visualisations for enzymes with large number of cut sites.

See sample output [here](docs/examples/)

### Authors

**RADSeq Team** - Arielle Gabriel. Mark Mendoza. Danielle Pamulaklakin. Francis Tablizo.
Project for the 2018 Internship Program in Bioinformatics

**IMBUE Q1** - Jobeth Domingo. Hannah Mae Magno. Marc Jermaine Pontiveros. Maria Rejane Nepacina.