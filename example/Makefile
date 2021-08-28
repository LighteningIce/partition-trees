path_only = `dirname $(realpath $(lastword $(MAKEFILE_LIST)))`

all: call-heuristics-main	call-heuristics		select-heuristics   memory-heuristics	memory-heuristics-test	hetro-heuristics   best-main

LIB_PATH = ${path_only}/../lib
INC_PATH = ${path_only}/../include
BIN_PATH = ${path_only}/

CPP = g++
PEDANTIC_PARANOID_FREAK =       -O0 -Wshadow -Wcast-align \
				-Waggregate-return -Wmissing-prototypes -Wmissing-declarations \
				-Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations \
				-Wmissing-noreturn -Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions
REASONABLY_CAREFUL_DUDE =	-Wall
NO_PRAYER_FOR_THE_WICKED =	-w -O2 
WARNINGS = $(REASONABLY_CAREFUL_DUDE)
CFLAGS = $(WARNINGS) -m64 -g -DNOASSERT -std=c++14 -fopenmp
#CFLAGS = $(PEDANTIC_PARANOID_FREAK) -std=c++11 -g
INCLUDES = -I${INC_PATH}
DEFS = 
LDADD =
LIBS = ${LIB_PATH}/heuristics.a 


call-heuristics-main: main4.cpp 
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@
call-heuristics: main5.cpp 
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@
select-heuristics:main6.cpp 
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@
hetro-heuristics:main7.cpp
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@	

memory-heuristics:main_memory.cpp
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@
memory-heuristics-test:main_memory_test.cpp
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@

best-main:best_main.cpp
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@


split_real_tree:split_real_tree.cpp
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@
select_thin:select_thin.cpp
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@

#gen-pro-data:gen_pro_data.cpp
#	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@

clean:
	rm -f ${BIN_PATH}/call-heuristics-main  *.o *~ 
	rm -f ${BIN_PATH}/call-heuristics  *.o *~ 
	rm -f ${BIN_PATH}/select-heuristics  *.o *~ 
	rm -f ${BIN_PATH}/split_real_tree  *.o *~ 
	rm -f ${BIN_PATH}/select_thin  *.o *~ 
	rm -f ${BIN_PATH}/memory-heuristics *.o *~ 
	rm -f ${BIN_PATH}/memory-heuristics-test *.o *~ 
	rm -f ${BIN_PATH}/hetro-heuristics *.o *~
	rm -f ${BIN_PATH/}best-main *.o *~

.PHONY: clean all

.SUFFIXES:
.SECONDARY:
