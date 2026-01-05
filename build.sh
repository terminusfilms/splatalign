#!/bin/bash
# Build script for SplatAlign standalone app

set -e

echo "================================"
echo "SplatAlign Build Script"
echo "================================"

# Check if we're in the right directory
if [ ! -f "splat_align.py" ]; then
    echo "Error: Run this script from the SplatAlign directory"
    exit 1
fi

# Create/activate virtual environment
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv .venv
fi

echo "Activating virtual environment..."
source .venv/bin/activate

# Install dependencies
echo "Installing dependencies..."
pip install -r requirements.txt
pip install pyinstaller

# Build
echo "Building standalone app..."
pyinstaller --clean splat_align.spec

echo ""
echo "================================"
echo "Build complete!"
echo "================================"
echo ""
echo "Output locations:"
echo "  macOS app:  dist/SplatAlign.app"
echo "  Executable: dist/SplatAlign"
echo ""
echo "To test:"
echo "  open dist/SplatAlign.app"
echo ""
echo "To distribute:"
echo "  - Zip the .app bundle for macOS users"
echo "  - For Windows, rebuild on a Windows machine"
