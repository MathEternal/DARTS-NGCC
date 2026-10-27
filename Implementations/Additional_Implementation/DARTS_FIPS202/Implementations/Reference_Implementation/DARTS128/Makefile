CC ?= /usr/bin/cc
#CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -march=native -mtune=native -std=c99

CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -march=native -mtune=native -std=gnu99
NISTFLAGS += -Wno-unused-result -O3

BIN_DIR = bin

.PHONY: all speed test clean

# Default target - KAT testing
all: $(BIN_DIR)/PQCgenKAT_sign

# Speed benchmark
speed: $(BIN_DIR)/test_speed

# All tests (performance and correctness)
test: $(BIN_DIR)/test_speed $(BIN_DIR)/test_Sign $(BIN_DIR)/test_mul

$(BIN_DIR):
	mkdir -p $@

# KAT generation
$(BIN_DIR)/PQCgenKAT_sign: PQCgenKAT_sign.c \
  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ PQCgenKAT_sign.c \
	  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
	  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c -lcrypto

# Speed benchmarking test
$(BIN_DIR)/test_speed: test/test_speed.c test/speed_print.c test/cpucycles.c \
  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ test/test_speed.c test/speed_print.c test/cpucycles.c \
	  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
	  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c -lcrypto

# Signature verification test
$(BIN_DIR)/test_Sign: test/test_Sign.c \
  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ test/test_Sign.c \
	  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
	  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c -lcrypto

# Polynomial multiplication test
$(BIN_DIR)/test_mul: test/test_mul.c \
  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ test/test_mul.c \
	  sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
	  polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c randombytes.c -lcrypto

# Shared library
shared: $(BIN_DIR)/libsign_ref.so

$(BIN_DIR)/libsign_ref.so: sign.c packing.c polyvec.c poly.c ntt.c reduce.c sampler.c \
  encoding.c polyfix.c polymat.c fft.c fixpoint.c fips202.c symmetric-shake.c | $(BIN_DIR)
	$(CC) -shared -fPIC $(CFLAGS) -o $@ sign.c packing.c polyvec.c poly.c ntt.c \
	  reduce.c sampler.c encoding.c polyfix.c polymat.c fft.c fixpoint.c \
	  fips202.c symmetric-shake.c -lcrypto

# Clean build artifacts
clean:
	-rm -rf $(BIN_DIR)
	-rm -f *.o test/*.o
	-rm -f PQCsignKAT_*.req PQCsignKAT_*.rsp