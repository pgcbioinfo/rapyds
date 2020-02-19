
## Restriction Site Associated DNA Python-Digested Simulation  (RApyDS)


**RApyDS** is a python script that performs in silico digestion as an aid for choosing a restriction enzyme for RADseq experiments. 

**RApyDS** can perform in silico digestion of a given genome/sequence or gc content. It can simulate single-enzyme or double-enzyme digestion on a single file or multiple FASTA files. After digestion, the program can match fragments to a given an annotation file and run alignment analysis to check for unique and repeat regions.

**RApyDS** provides a detailed report and visualization of the simulated digestion. 


### Documentation
- For the [full RApyDS program technical manual](docs/rapyds_manual.pdf)

### Requirements

- Python 3 or greater
- pip3
- BWA 0.7.12 (http://bio-bwa.sourceforge.net/)
- Firefox (for viewing the html files)
- Linux OS

### Usage

#### Installation
- Download the source from the Github repository and extract it. Or clone using ``git clone https://github.com/pgcbioinfo/rapyds.git``

- Enter the directory and run ``pip install -r requirement.txt`` or ``pip3 install -r requirements.txt`` whichever is applicable


#### Arguments
-  **-h, --help**  show this help message and exit
-  **-gc [GC]**    input GC frequency. Value must be between 0 and 1
-  **-dna [DNA]**  input estimated DNA length
-  **-i [I]**      directory containing the input files
-  **-pre [PRE]**  prefix of the input files (must match the file name of the sequence, annotation, and/or index files)
-  **-at [AT]**    target feature in annotation file (ex. gene region, exon, intron, etc) (default: gene) 
-  **-db [DB]**    resteiction enzyme dabatase file. Format per line: SbfI,CCTGCA|GG (default: database/re_db.txt)
-  **-re [RE]**    file containing list of restriction enzyme to be tested
-  **-min [MIN]**  minimum fragment size (default: 200)
-  **-max [MAX]**  maximum fragment size (default: 300)
-  **-bp [BP]**    base pair read length for mapping (default: 100)
-  **-p [P]**      RADSeq protocol: use ddrad for double digestion (default: orig)
-  **-o [O]**      output report file name (default: report)
-  **-t [T]**      number of processes (default 16)

Optional Flags:
-  **--skip_bwa**   skip BWA indexing and alignment
-  **--skip_graph**   skip cut site location histogram graphing
-  **--clean**     skip cleaning intermediate files after running



#### Common Usage
- To use RApyDS with a given genome file (and/or corresponding annotation file)

``python rapyds.py -i <input_directory> -pre <input_prefix> [other arguments/flags]``

- To use RApyDS by generating a genome file given a percent GC content

``python rapyds.py -gc <gc frequency from 0 to 1> -dna <length> [other args]``


#### Sample Run
Given an E.Coli FASTA file ``ecoli_seq.fasta`` with annotation file ``ecoli_seq.gff`` both located inside the directory ``ecoli``, to run the program using the origial RADSeq protocol, the command is shown below:

``python rapyds.py -i ecoli -pre ecoli_seq``

To run the same input using the ddRAD protocol and the list of RE pairings in enzyme/RE_list.txt:

``python rapyds.py -i ecoli -pre ecoli_seq -p ddrad -re enzyme/RE_list.txt``


#### Input Requirements
- User can only choose between having a directory input (using ``-i`` with ``-pre``or generating a sequence based on GC frequency (using ``-gc`` with ``-dna``

``python rapyds.py -i <input_dir> -pre <prefix> [other arguments]``
``python rapyds.py -gc <value from 0-1> -dna <int:length> [other arguments]``

-  The input directory must contain at least sequence file in the standard FASTA format having the extension of either .fasta, .fna, or .fa and the prefix or file name as specified in the ``-pre`` argument.

- The annotation and index files with filenames same as the prefix (ex. ``ecoli.fasta`` and ``ecoli.gff3``), must be inside the speicified input directory.

- The annotation file must be in the standard GFF3 format with the extension either be ``.gff``, ``.gff3``, or ``.gtf``. While the index files must have the extensions ``.amb``, ``.ann``, ``.bwt``, ``.pac``, and ``.sa``.

- User can only choose between ``orig`` (default) and ``ddrad`` as protocol (``-p``)

- Using the ddRAD protocol requires a ``-re`` argument or a list of restriction enzymes to be tested on.
``python rapyds.py -i <input_dir> -pre <prefix> -p ddrad -re <path/to/RE_list.txt>``

- In original RADSeq protocol, if no ``-re`` argument is given, by default the program uses all the restriction enzymes in the database

- Formats for the list of REs and the database are found in the next section.


#### Output

The output files of the program is a zip file named after the ``-o`` argument or by default, the ``report.zip``.
Inside the archive are 3 html files containing:
- Overview (``index.html``) - summary of the results
- Electrophoresis (``gel.html``) - electrophoresis simulation of up to 5 restriction enzymes of choice
- Cut Site Distribution (``cutsite.html``) - zoomable graphic vector of cut site locations in the sequence (optional: histogram of the cut site locations)

> It is advisable open the html files using Firefox. There is an issue with Google Chrome when opening local files. Only the overview file will work fine for any browser.

> Slowdown may be experienced when loading visualisations for enzymes with large number of cut sites.

See sample output [here](docs/examples/)

### Authors

**RADSeq Team** - Arielle Gabriel. Mark Mendoza. Danielle Pamulaklakin. Francis Tablizo.
Project for the 2018 Internship Program in Bioinformatics

**IMBUE Q1** - Jobeth Domingo. Hannah Mae Magno. Marc Jermaine Pontiveros. Maria Rejane Nepacina.