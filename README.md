POPBAM
======

POPBAM is a tool to perform evolutionary or population-based analyses of next-generation sequencing data. 
POPBAM takes a BAM file as its input and can compute many widely used evolutionary genetics measures in 
sliding windows across a genome.

INTRODUCTION
------------

This is the third beta release (0.3) of the POPBAM program. The source code has primarily been alpha tested 
on Debian and Red Hat based operating systems using the GNU C++ compiler. A Makefile is provided for the 
GNU C++ compiler with three compilation modes: release (default), debug, and profile.

OBTAINING THE SOURCE CODE
-------------------------

The POPBAM source code is available to download freely from GitHub at

[http://dgarriga.github.io/POPBAM/]

as either a gzipped tar file or a zipped file with the prefix

	dgarriga-POPBAM-<version>

where <version> is the current version number of POPBAM.

The only dependency of the POPBAM program is the zlib compression library and headers. 
Please insure these are installed on your system before attempting to compile POPBAM.

COMPILING THE SOURCE CODE
-------------------------

Once downloaded, you can extract the popbam/ directory using the command

	dgarriga-POPBAM-<version>.tar.gz

or 

	unzip dgarriga-POPBAM-<version>.zip

Move into the POPBAM source code directory by 

	cd popbam

Then you can simply type

	make

to build the source code.  If you have administrator privileges, you can automatically install 
the POPBAM executable into the /usr/local/bin directory by typing

	sudo make install

Lastly, you can clean the POPBAM directory by using

	make clean

GETTING HELP USING THE PROGRAM
------------------------------

Currently, the primary resource for helping users run POPBAM is a manpage.
If one has system administrator privileges, the popbam.1 file can be installed
into the system manpath or, alternatively, the POBPAM manpage can be viewed in
the current directory by typing

    man ./popbam.1

EXAMPLE DATA SET
----------------

To test the build of POPBAM, the user may download an example BAM file
from the web site:

[http://kimura.biology.rochester.edu/data/popbam/trial.bam]
	
This BAM files comprises a single read group each from nine lines of *Drosophila melanogaster*
from sub-Saharan Africa and one line from France. For the outgroup sequence, there is a single
*Drosophila mauritiana* read group. This example contains only reads that map to the X chromosome
of the *Drosophila melanogaster* reference genome (build 5.45).
