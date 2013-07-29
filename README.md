POPBAM
======

POPBAM is a tool to perform evolutionary or population-based analyses of next-generation sequencing data. 
POPBAM takes a BAM file as its input and can compute many widely used evolutionary genetics measures in 
sliding windows across a genome.

INTRODUCTION
------------

This is the third beta release (0.3) of the POPBAM program.
The source code has primarily been alpha tested on Debian and
RedHat based operating systems using the GNU C++ compiler.
A Makefile is provided for the GNU C++ compiler with three 
compilation modes: release (default), debug, and profile.

OBTAINING THE SOURCE CODE
-------------------------

The POPBAM source code is available to download freely from
SourceForge at

[http://sourceforge.net/projects/popbam/]

The source code is a gzipped tar file called

	popbam-<version>.tar.gz

where <version> is the current version number of POPBAM.

The only dependency of the POPBAM program is the zlib 
compression library and headers. Please insure these
are installed on your system before attempting to compile 
POPBAM.

COMPILING THE SOURCE CODE
-------------------------

Once downloaded, you can extract the popbam/ directory using
the command

	tar -xzf popbam-<version>.tar.gz

Move into the POPBAM source code directory by 

	cd popbam

Then you can simply type

	make

to build the source code.  If you have administrator privileges, you can
automatically install the POPBAM executable into the /usr/local/bin 
directory by typing

	sudo make install

Lastly, you can clean the POPBAM directory by using

	make clean

VIEWING THE MANPAGE
-------------------

If one has system administrator privileges, you can install the popbam.1 file
into a system manpath or, alternatively, the POBPAM manpage can be viewed in
the current directory by typing

    man ./popbam.1

