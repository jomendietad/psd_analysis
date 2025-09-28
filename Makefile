# C Compiler
CC = gcc

# Compiler Flags: -Wall (all warnings), -O2 (optimization), -pthread for threads
CFLAGS = -Wall -O2 -pthread

# --- LIBRARIES ---
# Libraries for long double version
LIBS_LD = -lfftw3l -lfftw3l_threads -lm

# Libraries for double version
LIBS_D = -lfftw3 -lfftw3_threads -lm

# Libraries for arbitrary-precision (MPFR) version
LIBS_MPFR = -lfftw3l -lfftw3l_threads -lm -lmpfr -lgmp

# --- PATHS ---
SRC_DIR = src
BIN_DIR = bin

# --- SOURCE FILES ---
SRCS_ANALYZER_LD = $(SRC_DIR)/main.c
SRCS_ANALYZER_D = $(SRC_DIR)/main_double.c
SRCS_ANALYZER_MPFR = $(SRC_DIR)/main_mpfr.c

# --- EXECUTABLE TARGETS ---
TARGET_ANALYZER_LD = $(BIN_DIR)/psd_analyzer_long_double
TARGET_ANALYZER_D = $(BIN_DIR)/psd_analyzer_double
TARGET_ANALYZER_MPFR = $(BIN_DIR)/psd_analyzer_mpfr

# Default rule: runs when 'make' is called
all: $(TARGET_ANALYZER_LD) $(TARGET_ANALYZER_D) $(TARGET_ANALYZER_MPFR)

# Rule to ensure the output directory exists
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# --- BUILD RULES ---

# Rule to build the main executable (long double)
$(TARGET_ANALYZER_LD): $(SRCS_ANALYZER_LD) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(SRCS_ANALYZER_LD) -o $(TARGET_ANALYZER_LD) $(LIBS_LD)

# Rule to build the double precision version
$(TARGET_ANALYZER_D): $(SRCS_ANALYZER_D) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(SRCS_ANALYZER_D) -o $(TARGET_ANALYZER_D) $(LIBS_D)

# Rule to build the arbitrary-precision (MPFR) version
$(TARGET_ANALYZER_MPFR): $(SRCS_ANALYZER_MPFR) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(SRCS_ANALYZER_MPFR) -o $(TARGET_ANALYZER_MPFR) $(LIBS_MPFR)

# Rule to clean generated files
clean:
	rm -rf $(BIN_DIR) data plots

.PHONY: all clean