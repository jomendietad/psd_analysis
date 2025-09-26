#!/bin/bash

# Función para manejar errores
handle_error() {
    echo "Error: Falló el comando en la línea $1"
    exit 1
}
trap 'handle_error $LINENO' ERR

# Cambiar al directorio del script (la raíz del proyecto)
# Esto hace que el script se pueda ejecutar desde cualquier lugar
cd "$(dirname "$0")"

# Argumento principal: run o clean
ACTION=${1:-run} # Por defecto es 'run'

if [ "$ACTION" = "run" ]; then
    echo "================================================="
    echo "=== MODO: CONSTRUIR Y EJECUTAR                ==="
    echo "================================================="

    # 1. Crear directorios necesarios
    echo "--> Creando directorios: bin/, data/, plots/..."
    mkdir -p bin data plots

    # 2. Compilar el código C usando el Makefile
    echo "--> Compilando código C con 'make'..."
    make all

    # 3. Ejecutar el analizador para generar los archivos .csv
    echo "--> Ejecutando analizador para generar datos..."
    ./bin/psd_analyzer

    # 4. Ejecutar el plotter de C/Gnuplot para mostrar las ventanas
    echo "--> Ejecutando plotter de C/Gnuplot (abrirá 3 ventanas)..."
    ./bin/plotter_c

    # 5. Ejecutar el plotter de Python para guardar las gráficas y mostrarlas
    echo "--> Ejecutando plotter de Python (guardará 3 .png y mostrará ventanas)..."
    python3 scripts/plotter.py

    echo "================================================="
    echo "=== PROCESO COMPLETADO                        ==="
    echo "================================================="
    echo "Los datos están en la carpeta 'data/'."
    echo "Las gráficas de Python se han guardado en la carpeta 'plots/'."

elif [ "$ACTION" = "clean" ]; then
    echo "================================================="
    echo "=== MODO: LIMPIAR PROYECTO                    ==="
    echo "================================================="
    
    # El Makefile ahora se encarga de toda la limpieza
    make clean
    
    echo "--> Directorios 'bin', 'data' y 'plots' eliminados."
    echo "================================================="
    echo "=== LIMPIEZA COMPLETADA                     ==="
    echo "================================================="

else
    echo "Argumento no válido: '$ACTION'"
    echo "Uso: ./build.sh [run|clean]"
    exit 1
fi