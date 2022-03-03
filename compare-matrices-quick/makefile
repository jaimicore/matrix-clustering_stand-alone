################################################################
## Compile compare-matrices-quick
##
## Author: Sebastien JAEGER <jaeger@ciml.univ-mrs.fr>

#CC      = g++ 
CC      = gcc
CCFLAGS = -O3 -Wall
SRC     = $(wildcard *.cpp) 
CSRC    = $(wildcard *.c)
OBJS    = $(SRC:.cpp=.o) $(CSRC:.c=.o)
APP     = compare-matrices-quick

compile: $(OBJS)
	$(CC) $(OBJS) -o $(APP) -lm

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

%.o: %.c
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f *.o $(APP)

all: clean compile

################################################################
## Test compare-matrices-quick (this C version)

test: test_quick

LTH_W=3
LTH_NCOR2=0.7
MATRIX_DB=DemoCompMat.txt
QUERY_MATRICES=DemoCompMat.txt
RESULT=result
test_quick:
	time compare-matrices-quick -o ${RESULT}_quick.tab -file1 ${QUERY_MATRICES} -file2 ${MATRIX_DB} -lth_ncor2 ${LTH_NCOR2} -lth_w ${LTH_W}
	@echo "	-o ${RESULT}_quick.tab "
	text-to-html -i ${RESULT}_quick.tab -o ${RESULT}_quick.html
	@echo "	-o ${RESULT}_quick.html "
