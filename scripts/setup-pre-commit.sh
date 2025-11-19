#!/bin/bash

# setup-pre-commit.sh
# Sets up pre-commit hooks for the project

echo "Setting up pre-commit hooks..."

# Check if pre-commit is installed
if ! command -v pre-commit &> /dev/null; then
    echo "pre-commit not found. Installing..."

    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        if command -v brew &> /dev/null; then
            brew install pre-commit
        else
            echo "Homebrew not found. Installing pre-commit with pip..."
            pip install pre-commit
        fi
    else
        # Linux/Other
        pip install pre-commit
    fi
fi

# Install the git hooks
pre-commit install

# Run against all files to check current state
echo ""
echo "Running pre-commit on all files..."
pre-commit run --all-files

echo ""
echo "Pre-commit hooks installed successfully!"
echo "They will run automatically on git commit."
echo ""
echo "To run manually: pre-commit run --all-files"
echo "To skip hooks: git commit --no-verify"
