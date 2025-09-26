# Nombre del compilador de C
CC = gcc

# Flags para el compilador: -Wall (todas las advertencias), -O2 (optimización)
CFLAGS = -Wall -O2

# Bibliotecas a enlazar (-l)
LIBS = -lfftw3l -lfftw3l_threads -lpthread -lm

# --- RUTAS ---
# Directorio de código fuente
SRC_DIR = src
# Directorio de salida para ejecutables
BIN_DIR = bin

# Archivo(s) de código fuente
SRCS_ANALYZER = $(SRC_DIR)/main.c
SRCS_PLOTTER = $(SRC_DIR)/plotter.c

# Nombre del archivo ejecutable de salida
TARGET_ANALYZER = $(BIN_DIR)/psd_analyzer
TARGET_PLOTTER = $(BIN_DIR)/plotter_c

# Regla por defecto: se ejecuta al escribir 'make'
all: $(TARGET_ANALYZER) $(TARGET_PLOTTER)

# Regla para asegurar que los directorios de salida existan
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# Regla para construir el ejecutable principal
$(TARGET_ANALYZER): $(SRCS_ANALYZER) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(SRCS_ANALYZER) -o $(TARGET_ANALYZER) $(LIBS)

# Regla para construir el ejecutable del plotter de Gnuplot
$(TARGET_PLOTTER): $(SRCS_PLOTTER) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(SRCS_PLOTTER) -o $(TARGET_PLOTTER)

# Regla para limpiar los archivos generados
clean:
	rm -rf $(BIN_DIR) data plots

.PHONY: all clean