#!/bin/bash
# Build script for C++ components of the audio fingerprinting service

# Create directories if they don't exist
mkdir -p src/cpp/fingerprinter
mkdir -p src/cpp/bindings
mkdir -p build

# Make sure all files exist (even if empty)
touch src/cpp/fingerprinter/peaks.cpp
touch src/cpp/fingerprinter/fingerprint.cpp

# Install required packages
echo "Installing required packages..."
pip install -r requirements.txt

# Build the C++ extension module
echo "Building C++ extension module..."
python setup.py build_ext --inplace

# Test the module
echo "Testing if the module was built successfully..."
python -c "import audiofp_core; print('Successfully imported audiofp_core')"

if [ $? -eq 0 ]; then
    echo "Build completed successfully!"
else
    echo "Build failed. Please check the error messages above."
    exit 1
fi