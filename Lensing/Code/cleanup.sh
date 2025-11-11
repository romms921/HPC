#!/bin/bash

# This script cleans the Run and Output directories.

cd ..
echo "Changing to the 'Run' directory..."
cd Simulations/Run
if [ $? -eq 0 ]; then
    echo "Clearing all files in 'Run'..."
    rm -f *
    cd ..
    echo "Returned to the parent directory."
else
    echo "Error: Could not change to the 'Run' directory."
    exit 1
fi

echo "Changing to the 'Output' directory..."
cd Output
if [ $? -eq 0 ]; then
    echo "Clearing all files in 'Output'..."
    rm -f *
    cd ..
    echo "Returned to the parent directory."
else
    echo "Error: Could not change to the 'Output' directory."
    exit 1
fi

echo "Cleanup script finished successfully."