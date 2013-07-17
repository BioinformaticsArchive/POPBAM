# Makefile for the popbam program
CXX=              g++
CXXFLAGS=         -D_FILE_OFFSET_BITS=64 -std=c++0x
C_RELEASE_FLAGS=  -Wno-unused -Wno-sign-compare -Wno-write-strings -O2
C_DEBUG_FLAGS=    -Wall -ggdb -DDEBUG
C_PROFILE_FLAGS=  -Wall -O2 -pg
OBJS=              popbam.o pop_utils.o pop_sample.o pop_tree.o \
                   pop_snp.o pop_nucdiv.o pop_ld.o pop_sfs.o pop_diverge.o \
                   pop_haplo.o pop_svnversion.o getopt_pp.o gamma.o\
                   bam_aux.o bam.o bam_index.o faidx.o kstring.o bgzf.o \
                   sam_header.o bam_import.o razf.o sam.o bam_pileup.o
PROG=              popbam
LIBFLAGS=          -lz -lm
INSTALL_DIR=       /usr/local/bin

all: CXXFLAGS +=  $(C_RELEASE_FLAGS)
all: $(PROG)

pop_svnversion.cpp: force
	echo -n 'const char* svn_version(void) { const char* SVN_Version = "' \
									   > pop_svnversion.cpp
	svnversion -n .                   >> pop_svnversion.cpp
	echo '"; return SVN_Version; }'   >> pop_svnversion.cpp
force: ;


.PHONY: debug
debug: CXXFLAGS += $(C_DEBUG_FLAGS)
debug: $(PROG)

.PHONY: release
release: CXXFLAGS += $(C_RELEASE_FLAGS)
release: $(PROG)

.PHONY: profile
profile: CXXFLAGS += $(C_PROFILE_FLAGS)
profile: $(PROG)

$(PROG): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBFLAGS)

%.o : %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.cpp.o: pop_svnversion.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

.PHONY: install
install: $(PROG)
	install -m 755 $(PROG) $(INSTALL_DIR)

.PHONY: clean
clean:
	rm -rf $(PROG) $(OBJS) pop_svnversion.cpp

