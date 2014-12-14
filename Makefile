CC = mpicc
CFLAGS = -O2 -Wall -Werror -Wformat-security -Wignored-qualifiers -Winit-self -Wswitch-default -Wfloat-equal -Wshadow -Wpointer-arith -Wtype-limits -Wempty-body -Wlogical-op -Wstrict-prototypes -Wold-style-declaration -Wold-style-definition -Wmissing-parameter-type -Wmissing-field-initializers -Wnested-externs -Wno-pointer-sign -std=gnu99
LDFLAGS = -s -fopenmp -lm
CFILES = main.c common/common.c common/parsecfg.c common/matrix.c common/random.c
HFILES = common/common.h common/parsecfg.h common/matrix.h common/random.h
OBJECTS = $(CFILES:.c=.o)
TARGET = nmf

all: $(TARGET)
$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@
include depth.make
clean:
	rm $(TARGET) $(OBJECTS)
depth.make: $(CFILES) $(HFILES)
	gcc -MM $(CFILES) > depth.make
