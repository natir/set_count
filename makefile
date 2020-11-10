CC=g++
AR=gcc-ar

DEBUG ?= 0
ifeq ($(DEBUG), 1)
        CFLAGS+=-O0
        WARNS= -Wextra -Wpedantic -Wno-format -Wswitch-default -Wswitch-enum -Wfloat-equal -Wconversion -Wsign-conversion \
        -Wold-style-cast -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 \
        -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 \
        -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers \
        -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wnarrowing -Wsuggest-final-methods \
        -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized \
        -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion \
        -Wduplicated-cond -Werror=write-strings -Werror=return-type -Wredundant-decls \
        -Werror=volatile-register-var -Wsuggest-final-types
        DEBUG_SYMS=1
else
        CFLAGS+=-O3 -fno-fat-lto-objects -flto=jobserver -march=native -mtune=native -mcmodel=large
        LDFLAGS+=-fuse-linker-plugin
        WARNS=-Wfatal-errors
endif

ASSERTS ?= $(DEBUG)
ifeq ($(ASSERTS), 1)
	CFLAGS+=-DDEBUG
else
	CFLAGS+=-DNDEBUG
endif

DEBUG_SYMS ?= 1
ifeq ($(DEBUG_SYMS), 1)
	CFLAGS+=-g
	LDFLAGS+=-g
endif

ifeq ($(VERBOSE), 1)
	SHELL=sh -x
endif

WARNS+= -Wall
CFLAGS+=-std=c++14 -pipe -fopenmp ${WARNS}
LDFLAGS+=-lpthread -lz  thirdparty/MQF/build/src/libMQF.a thirdparty/MQF/build/ThirdParty/stxxl/lib/libstxxl.a thirdparty/Brisk/build/brisk/libbrisk_lib.a

INCS=-isystem thirdparty/BBHash -isystem thirdparty/kseqpp/include -isystem thirdparty/MQF/include -isystem thirdparty/Brisk/brisk -I include/set_count 
SUBMODULE_TOKEN=thirdparty/BBHash/makefile thirdparty/kseqpp/CMakeLists.txt thirdparty/MQF/CMakeLists.txt thirdparty/Brisk/CMakeLists.txt
SUBMODULE_REQ=thirdparty/kseqpp/include/kseq++/config.hpp thirdparty/MQF/build/src/libMQF.a thirdparty/MQF/build/ThirdParty/stxxl/lib/libstxxl.a thirdparty/Brisk/build/brisk/libbrisk_lib.a

CPPS = $(wildcard src/*.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:%.o=%.d)
SET_COUNT_OBJ = 

EXEC=set_count

all: $(EXEC) $(LIB)

bin: submodule_build
	@mkdir --parent ./bin

-include $(DEPS)

set_count: bin bin/set_count
bin/set_count: $(SET_COUNT_OBJ) src/set_count.o
	@echo "[LD] $@"
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

%.o: %.cpp $(SUBMODULE_TOKEN)
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) $(INCS) -MMD -o $@ -c $<

$(SUBMODULE_TOKEN):
	git submodule update --init

submodule_build: $(SUBMODULE_TOKEN) $(SUBMODULE_REQ)

thirdparty/kseqpp/include/kseq++/config.hpp: thirdparty/kseqpp/include//kseq++/config.hpp.in
	cd thirdparty/kseqpp/ && cmake . && make

thirdparty/MQF/build/src/libMQF.a: thirdparty/MQF/CMakeLists.txt
	mkdir -p thirdparty/MQF/build && cd thirdparty/MQF/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

thirdparty/MQF/build/ThirdParty/stxxl/lib/libstxxl.a: thirdparty/MQF/CMakeLists.txt
	mkdir -p thirdparty/MQF/build && cd thirdparty/MQF/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

thirdparty/Brisk/build/brisk/libbrisk_lib.a: thirdparty/Brisk/CMakeLists.txt
	mkdir -p thirdparty/Brisk/build && cd thirdparty/Brisk/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

clean:
	@echo "[clean]"
	@rm -rf bin $(LIB) $(OBJS) $(DEPS) $(SUBMODULE_REQ)

clean_all: clean
	@echo "[clean all]"
	@rm -rf $(SUBMODULE_REQ)

rebuild: clean
	@$(MAKE) -s all

rebuild_all: clean_all
	@$(MAKE) -s all

check_buildsys: $(SUBMODULE_TOKEN)
	$(CC) --version
	$(AR) --version
	@echo CFLAGS=$(CFLAGS)
	@echo LDFLAGS=$(LDFLAGS)

.PHONY: all $(EXEC) tests clean rebuild
