VERSION=$(shell perl -I../perl/lib -MSanger::CGP::AlleleCount -e 'print Sanger::CGP::AlleleCount->VERSION;')

#Compiler
CC = gcc -O3 -DALLELECOUNTER_VERSION='"$(VERSION)"'

#compiler flags
# -g adds debug info to the executable file
# -Wall turns on most warnings from compiler
CFLAGS = -g -Wall

#Location of samtools/htslib libraries
HTSLOC?=$(HTSLIB)

HTSTMP?=./htslib_tmp
prefix?=/usr/local/

#Define locations of header files
OPTINC?=-I$(HTSLOC)/
INCLUDES= -Isrc/ $(OPTINC) -rdynamic

JOIN_INCLUDES= -I$(prefix)/include
CAT_LFLAGS= -L$(prefix)/lib

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS?= -L$(HTSTMP)

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:
LIBS=-lhts -lpthread -lz -lbz2 -llzma -lm -ldl

# define the C source files
SRCS=./src/bam_access.c

#Define test sources
TEST_SRC=$(wildcard ./tests/*_tests.c)
TESTS=$(patsubst %.c,%,$(TEST_SRC))

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
MD := mkdir -p

OBJS = $(SRCS:.c=.o)

#Build target executable
COUNTER_TARGET=./bin/alleleCounter

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean coverage test make_htslib_tmp remove_htslib_tmp

.NOTPARALLEL: test

all: clean make_bin make_htslib_tmp $(COUNTER_TARGET) test remove_htslib_tmp
	@echo  Binaries have been compiled.

$(COUNTER_TARGET): $(OBJS)
	$(CC) $(JOIN_INCLUDES) $(INCLUDES) $(CFLAGS) -o $(COUNTER_TARGET) $(OBJS) $(LFLAGS) $(CAT_LFLAGS) $(LIBS) ./src/alleleCounter.c

#Unit Tests
#test: $(COUNTER_TARGET)
##test: CFLAGS += $(JOIN_INCLUDES) $(INCLUDES) $(OBJS) $(LFLAGS) $(LIBS) $(CAT_LFLAGS)
#test: $(TESTS)
#	sh ./tests/runtests.sh

#Unit tests with coverage
#coverage: CFLAGS += --coverage
#coverage: test

make_bin:
	$(MD) ./bin

make_htslib_tmp:
	$(MD) $(HTSTMP)
	#Do some magic to ensure we compile ALLELECOUNT with the static libhts.a rather than libhts.so
	ln -fs $(HTSLOC)/libhts.a $(HTSTMP)/libhts.a

remove_htslib_tmp:
	@echo remove tmp hts location
	-rm -rf $(HTSTMP)

valgrind:
	VALGRIND="valgrind --log-file=/tmp/valgrind-%p.log" $(MAKE)


# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(JOIN_INCLUDES) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) ./src/*.o *~ $(COUNTER_TARGET) ./bin/* ./tests/tests_log $(TESTS) ./src/*.gcda ./src/*.gcov ./src/*.gcno *.gcda *.gcov *.gcno ./tests/*.gcda ./tests/*.gcov ./tests/*.gcno

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
