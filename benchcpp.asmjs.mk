EMCC=$(EMSCRIPTEN)/emcc --memory-init-file 0 -s PRECISE_F32=1 -DNDEBUG=1 -I./
HEAP=570425344

OBJECTS = \
run.bc \
base.bc

all: benchcpp.js

%.bc: %.cpp
	$(EMCC) -O2 -g -s ASM_JS=1 $< -o $@

benchcpp.js.bc: $(OBJECTS)
	$(LLVM)/llvm-link -o $@ $(OBJECTS)

benchcpp.js: benchcpp.js.bc
	$(EMCC) -O2 -g -s ASM_JS=1 -s TOTAL_MEMORY=$(HEAP) $< -o $@

clean:
	rm benchcpp.js benchcpp.js.bc $(OBJECTS)
