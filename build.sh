#!/bin/bash

# build.sh: Master script to compile, run, and plot.

# Function to display help
show_help() {
    echo "Usage: ./build.sh [command]"
    echo "Commands:"
    echo "  run      Compiles all C sources, runs the executables to generate data,"
    echo "           and then runs the Python script for plotting and analysis."
    echo "  compile  Only compiles the C source files into the bin/ directory."
    echo "  clean    Removes all generated directories (bin/, data/, plots/)."
    echo "  help     Shows this help message."
}

# Function to compile all C sources
compile_all() {
    echo "--- Compiling all C source files... ---"
    make all
    if [ $? -ne 0 ]; then
        echo "--- Makefile compilation failed. Aborting. ---"
        exit 1
    fi
    echo "--- Compilation successful. Executables are in bin/. ---"
}

# Function to run all executables
run_all() {
    echo "--- Running executables to generate data... ---"
    # Ensure data directory exists
    mkdir -p data

    # Run with resource timing
    echo "--- Running Standard Precision (double) version... ---"
    /usr/bin/time -v -o data/resource_usage_d.txt ./bin/psd_analyzer_double

    echo "--- Running High Precision (long double) version... ---"
    /usr/bin/time -v -o data/resource_usage_ld.txt ./bin/psd_analyzer_long_double
    
    echo "--- Running Arbitrary Precision (MPFR) version... ---"
    echo "--- NOTE: This may take a significantly longer time. ---"
    /usr/bin/time -v -o data/resource_usage_mpfr.txt ./bin/psd_analyzer_mpfr

    echo "--- All data generation complete. ---"
}

# Function to run the python plotter
run_plotter() {
    echo "--- Running Python script for analysis and plotting... ---"
    python3 scripts/plotter.py
    echo "--- Python script finished. ---"
}

# Function to clean the project
clean_project() {
    echo "--- Cleaning project directories... ---"
    make clean
    echo "--- Clean complete. ---"
}

# Main script logic
case "$1" in
    run)
        compile_all
        run_all
        run_plotter
        ;;
    compile)
        compile_all
        ;;
    clean)
        clean_project
        ;;
    help|*)
        show_help
        ;;
esac

exit 0